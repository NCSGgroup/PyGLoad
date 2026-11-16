"""
Terrestrial water storage from GNSS, compared with GRACE, for one region (settings in config.py):

  1. preprocessing  : read the NGL series, offsets (steps file + detection), trajectory model, quality control,
                      monthly means, neighbour test, common reference period          -> qc.csv, dataset.npz
  2. inversion      : monthly regularised inversion with the chosen basis                 -> inversion.npz
  3. comparison     : regional mean EWH and mass, GNSS against GRACE mascons             -> series.csv, summary.txt

    python studies/GNSS_GRACE/download.py     # once, needs internet
    python studies/GNSS_GRACE/run.py
"""
import time

import numpy as np

import pygwater as pgw
from pygwater.config import FileTool
from pygwater.preprocessing import GNSSNetwork, GNSSDataset, read_ngl_steps, ngl
from pygwater.inversion.workflow import MonthlyInversion
from pygwater.observations.grace import Mascon, compare, remove_far_field, cm_to_cf
from pygwater.mesh.land import land_fraction

from config import CONFIG as C


def decyear(month):
    return int(month[:4]) + (int(month[5:7]) - 1) / 12


def reference_ring():
    """names of the far reference stations for the common mode (cme_reference as a dict), else []"""
    ref = C.get('cme_reference')
    if not isinstance(ref, dict):
        return []
    holdings = FileTool.get_input_dir('gnss/NGL') / 'DataHoldings.txt'
    return ngl.select_reference_ring(ngl.read_holdings(holdings), ref['center'], ref['min_km'], ref['max_km'],
                                     C['start'], C['end'], C['min_years'])


def preprocess(out):
    src = FileTool.get_input_dir('gnss/NGL')
    files = sorted((src / 'tenv3').glob('*.tenv3'))
    if not files:
        raise FileNotFoundError('no .tenv3 files in %s; run download.py first' % (src / 'tenv3'))
    steps = read_ngl_steps(src / 'steps.txt', as_events=True) if (src / 'steps.txt').exists() else None
    ring = reference_ring()
    net = GNSSNetwork.from_files(files)
    if ring:                                             # far reference stations for the common mode are kept too
        box, keep = C['box'], []
        for st in net.stations:
            if (box[0] <= st.lat <= box[1] and box[2] <= st.lon <= box[3]) or st.name in ring:
                keep.append(st)
            else:
                net.rejected[st.name] = 'outside the region'
        net.stations = keep
    else:
        net.select(box=C['box'])
    corrections = None
    frame = C.get('cm_to_cf')
    if frame not in (False, None, 'grace', 'ngl'):
        raise ValueError("cm_to_cf must be 'grace', 'ngl' or False")
    if C.get('load_models') or frame == 'ngl':
        corrections = ngl.load_corrections(src / 'loadpred', net.names(), C.get('load_models') or (),
                                           cmcf_folder=src / 'loadpred_cmcf' if frame == 'ngl' else None)
        parts = list(C.get('load_models') or ()) + (['CM->CF (NGL)'] if frame == 'ngl' else [])
        print('corrections (%s) for %d of %d stations' % ('+'.join(parts), len(corrections), len(net.names())))
    net.clean(steps=steps, period=(decyear(C['start']), decyear(C['end']) + 1 / 12), min_years=C['min_years'],
              corrections=corrections)
    net.daily() if C.get('step') == 'day' else net.monthly()
    if C['spatial_filter_km']:
        net.spatial_filter(radius_km=C['spatial_filter_km'])
    ds = net.to_dataset(components=C['components'], reference=C.get('gnss_reference', C['reference']))
    if frame == 'grace':
        # geocentre motion of the same GRACE solution as the comparison (degree 1 only: no regional signal)
        lln = pgw.LoveNumber().config(lmax=C['lmax'], method=pgw.LLN_Data.Wang).get_Love_number()
        g = Mascon.read(FileTool.get_input_dir(C['grace_file'])).anomalies(C['reference'])
        ds = cm_to_cf(ds, g, lln)
        del g
        print('CM -> CF with the geocentre motion of %s' % C['grace_file'])
    net.save_qc(out / 'qc.csv')
    ds.save(out / 'dataset.npz')
    print('stations: %d read, %d kept; %d %s' % (len(files), len(ds.names), len(ds.months),
                                                'days' if C.get('step') == 'day' else 'months'))
    return ds


def main():
    out = FileTool.get_output_dir('GNSS_GRACE/%s' % C['name'])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    ds = GNSSDataset.load(out / 'dataset.npz') if (out / 'dataset.npz').exists() and C.get('reuse', False) \
        else preprocess(out)
    print('preprocessing: %.0f s' % (time.time() - t0))

    lat_min, lat_max, lon_min, lon_max = C['box']
    r, b = C['resolution'], C['buffer']
    lat1 = np.arange(lat_min - b + r / 2, lat_max + b, r)
    lon1 = np.arange(lon_min - b + r / 2, lon_max + b, r)
    lln = pgw.LoveNumber().config(lmax=C['lmax'], method=pgw.LLN_Data.Wang).get_Love_number()
    lln.convert(pgw.Frame.CF)
    opts = dict(subgrid=r / 4) if C['kernel'] == 'cluster_disk' else {}
    # cells on land only (no station constrains the ocean, where the inversion would otherwise place water)
    land = land_fraction(lat1, lon1, r) >= C['min_land'] if C['min_land'] else None
    cells = pgw.GreenFunctionBasis.from_grid(lat1, lon1, r, mask=land, kernel=C['kernel'], **opts)
    R = pgw.regularization.curvature(cells.lat, cells.lon, r)
    grace_file = FileTool.get_input_dir(C['grace_file'])
    if C.get('far_field'):
        # water outside the cells, predicted from GRACE, is removed first: GNSS alone cannot tell it from water inside
        f, t0 = C['far_field'], time.time()
        gfar = Mascon.read(grace_file, box=(lat_min - f, lat_max + f, lon_min - f, lon_max + f)).anomalies(C['reference'])
        ds = remove_far_field(ds, gfar, lln, exclude=(lat_min - b, lat_max + b, lon_min - b, lon_max + b))
        print('far field: GRACE load beyond %g deg of the box (up to %g deg) removed (%.0f s)' % (b, f, time.time() - t0))
    if C.get('cme_reference'):
        ring = reference_ring()
        ds = ds.remove_common_mode(ring if ring else C['cme_reference'])
        print('common mode removed by stacking %d reference stations' % len(ds.meta['common_mode_stations']))
        if ring:
            ds = ds.select(box=C['box'])               # the far reference stations are not inverted
    t0 = time.time()
    inv = MonthlyInversion(cells, lln, ds, components=C['components'], R=R, common_mode=C.get('common_mode', ()))
    # daily: one alpha from every 15th day (the choice is stable, and GCV on every day would take long)
    epochs = inv.months(5)
    if C.get('window'):                                     # e.g. the days around an event
        epochs = [e for e in epochs if C['window'][0] <= e[:len(C['window'][0])] <= C['window'][1]]
    alpha = inv.choose_alpha(method=C['alpha_method'],
                             months=epochs[::15] if C.get('step') == 'day' and len(epochs) > 60 else epochs)
    res = inv.run(alpha, months=epochs, lb=0.0 if C['nonnegative'] else None)
    res.save(out / 'inversion.npz')
    print('inversion: %d cells, alpha %.2e, %d %s (%.0f s)' % (cells.n_param, alpha, len(res.months),
                                                             'days' if C.get('step') == 'day' else 'months',
                                                                 time.time() - t0))

    glat = np.arange(lat_min + 0.125, lat_max, 0.25)
    glon = np.arange(lon_min + 0.125, lon_max, 0.25)
    inside = land_fraction(glat, glon, 0.25) >= 0.5          # comparison region: land in the box
    gnss = res.regional_series(glat, glon, inside)
    mass = res.regional_mass(glat, glon, inside)
    rows = ['month,gnss_ewh_mm,gnss_mass_gt,grace_ewh_mm,n_stations']
    if grace_file.exists():
        g = Mascon.read(grace_file, box=(lat_min - 1, lat_max + 1, lon_min - 1, lon_max + 1)).anomalies(C['reference'])
        gseries = g.regional_series(glat, glon, inside)
        gs = dict(zip(g.months, gseries))
        lines = []
        for detrend in (False, True):
            # station trends were removed (tectonics), so GNSS has none; the detrended line compares like with like
            stats = compare(res.months, gnss, g.months, gseries, detrend=detrend)
            lines.append('GNSS vs GRACE over the land in %s, %d common months%s: r = %.2f, RMS = %.1f mm, NSE = %.2f, '
                         'amplitude ratio = %.2f' % (C['box'], stats['n'], ', both detrended' if detrend else '',
                                                     stats['r'], 1e3 * stats['rms'], stats['nse'],
                                                     stats['amplitude_ratio']))
        summary = '\n'.join(lines) if stats['n'] >= 3 else \
            'no comparison with GRACE: fewer than 3 common months (e.g. inside the 2017-2018 GRACE/GRACE-FO gap)'
    else:
        gs, summary = {}, 'no GRACE file at %s: comparison skipped' % grace_file
    for m, e, ms, n in zip(res.months, gnss, mass, res.n_stations):
        rows.append('%s,%.2f,%.2f,%s,%d' % (m, 1e3 * e, ms, '%.2f' % (1e3 * gs[m]) if m in gs else '', n))
    (out / 'series.csv').write_text('\n'.join(rows) + '\n')
    (out / 'summary.txt').write_text(summary + '\n')
    print(summary)
    if C.get('step') == 'day':
        lo, hi = C.get('gnss_reference', C['reference'])
        base = [i for i, e in enumerate(res.months) if lo <= e[:len(lo)] and e[:len(hi)] <= hi]
        m = mass - (np.mean(mass[base]) if base else 0.0)
        m3 = res.running_mean(3).regional_mass(glat, glon, inside)
        m3 = m3 - (np.mean(m3[base]) if base else 0.0)
        k, k3 = int(np.nanargmax(m)), int(np.nanargmax(m3))
        print('water stored over the land in the box, relative to %s to %s: peak %.1f Gt on %s (3-day mean: %.1f Gt on '
              '%s); scatter before %.1f Gt (3-day mean: %.1f Gt)'
              % (lo, hi, m[k], res.months[k], m3[k3], res.months[k3], np.std(m[base]) if base else np.nan,
                 np.std(m3[base]) if base else np.nan))
    print('results in %s' % out)


if __name__ == '__main__':
    main()

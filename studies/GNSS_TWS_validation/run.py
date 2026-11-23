"""
GNSS-TWS validation: how well does the elastic loading predicted from a monthly terrestrial water storage (TWS) grid
match the monthly vertical GNSS displacements at 4,807 stations (2003-2016)? Settings in config.py.

  1. GNSS        cleaned monthly vertical series per frame (igs20cf: IGS20 - cmcf - NTAL - NTOL; igs14: IGS14 -
                 NTAL - NTOL), read from processed/ or rebuilt from the NGL daily series (--gnss raw)
  2. products    each TWS grid as monthly 0.5-degree anomalies about its own mean over its months on the common
                 months, on the common land domain (zero elsewhere: ocean, Greenland, Antarctica); spherical-harmonic
                 analysis to degree 360                                                              -> cache/
  3. forward     vertical displacement at the stations: SH design matrix, PREM load Love numbers in the CF frame,
                 degree 0 included
  4. scores      per station and product, medians with station-bootstrap 5-95 % ranges, globally, without Japan,
                 per WMO region and by distance to the nearest basin-unit border                     -> output/<frame>/

    python studies/GNSS_TWS_validation/run.py                                  # CSR + JPL, both frames
    python studies/GNSS_TWS_validation/run.py --products CSR,JPL,GLWS4 --frames igs20cf
    python studies/GNSS_TWS_validation/run.py --tws MYTWS=/path/my_tws.nc      # your own monthly TWS grid
    python studies/GNSS_TWS_validation/run.py --gnss raw                       # also redo the GNSS cleaning

Scores, per station on the months both series have (both demeaned there; then the mean and trend of a mean + trend +
annual + semi-annual least-squares fit are removed from each series, because the cleaning removed the GNSS trend):
    r_total        Pearson r
    r_seas         r of the 12 calendar-month climatologies (every calendar month needs >= 2 values)
    r_nonseas      r of the deseasonalised series, model trend kept;  r_nonseas_dt: model trend removed too
    amp_ratio      annual amplitude model / GNSS;  phase_diff_d: day of the annual maximum, model - GNSS (+ = later)
    wrms_red       1 - WRMS(GNSS - model) / WRMS(GNSS), weights 1 / sigma^2 of the GNSS monthly means
                   (wrms_red_nodt: model trend kept)
"""
import argparse
import sys
import time
import warnings
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import netCDF4 as nc
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repository root, so that `src` is importable
from src import LoveNumber, LLN_Data, Frame, Displacement, SphericalHarmonicBasis
from src.preprocessing import GNSSNetwork, StationSeries, read_tenv3, read_ngl_steps
from src.preprocessing.series import _decyear_to_month
from src.utils.geomath import GeoMathKit

import config
from config import CONFIG as C, FRAMES, PRODUCTS

ALLM = ['%d-%02d' % (y, m) for y in range(2003, 2017) for m in range(1, 13)]   # 2003-01 .. 2016-12
LAT = 89.75 - 0.5 * np.arange(360)                                              # grid of the study: north -> south
LON = -179.75 + 0.5 * np.arange(720)
SCORES = ['r_total', 'r_seas', 'r_nonseas', 'r_nonseas_dt', 'amp_ratio', 'phase_diff_d', 'wrms_red', 'wrms_red_nodt']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# ---------------------------------------------------------------------------------------------------------------
# inputs of input/validation/
# ---------------------------------------------------------------------------------------------------------------
def domain():
    """common land domain (360 x 720 bool, rows north -> south, columns -180 -> 180) and the common months"""
    with nc.Dataset(config.validation_dir('land_domain_0.5deg.nc')) as f:
        lat, lon = f['lat'][:], f['lon'][:]
        assert np.allclose(lat, LAT) and np.allclose(lon, LON), 'land_domain_0.5deg.nc: unexpected grid'
        land = np.asarray(f['land'][:]).astype(bool)
    months = [l.strip() for l in open(config.validation_dir('common_months.txt')) if l.strip() and l[0] != '#']
    return land, months


def stations():
    """all candidate stations; the used ones have used == True"""
    S = pd.read_csv(config.validation_dir('stations.csv'), index_col=0)
    S['used'] = S['used'].astype(bool)
    return S


# ---------------------------------------------------------------------------------------------------------------
# 1. GNSS
# ---------------------------------------------------------------------------------------------------------------
def read_loadpred(path):
    """NGL loadpredictions file (ntal.cf/<site>, ntol.cf/<site>, cmcf/<site>): decimal year and up [m]"""
    t, u = [], []
    for line in open(path):
        p = line.split()
        if len(p) == 6 and p[0] != 'site':
            try:
                t.append(float(p[2])); u.append(float(p[5]))
            except ValueError:
                pass
    return np.array(t), np.array(u)


_STEPS = None


def _init_worker(steps):
    global _STEPS
    _STEPS = steps


def _clean_one(args):
    """daily NGL series of one station -> cleaned monthly vertical series and the monthly corrections [m]

    The trajectory model (bias, rate, NGL offsets, annual, semi-annual) is fitted by weighted least squares with
    iterative 5 x 1.4826 MAD outlier rejection; bias, rate and offsets are removed (seasonal terms and residual kept).
    Only that: the extra data selection, offset detection and quality control of GNSSNetwork.clean() are switched
    off. Monthly means need >= 10 days. Each correction (cmcf, NTAL, NTOL) is interpolated to the daily epochs and
    averaged over exactly the days of each monthly mean (= subtracting daily, then averaging)."""
    site, frame = args
    f = FRAMES[frame]
    try:
        s = read_tenv3(str(config.ngl_dir('tenv3/%s/%s.tenv3' % (f['tenv3'], site))))
    except Exception as e:
        return site, 'read error: %s' % e
    k = (s.t >= C['period'][0]) & (s.t < C['period'][1])
    if k.sum() < C['min_daily']:
        return site, 'fewer than %d daily solutions' % C['min_daily']
    s = StationSeries(s.name, s.lat, s.lon, s.t[k], {c: v[k] for c, v in s.xyz.items()},
                      {c: v[k] for c, v in s.sig.items()})
    net = GNSSNetwork([s])
    net.clean(steps=_STEPS, remove='trend+steps', outlier=C['outlier'], min_years=0.0, min_completeness=0.0,
              max_sigma=np.inf, max_residual=np.inf, detect_jumps=False,
              segment_options=dict(min_segment=0.0, min_between=0.0))
    net.monthly(min_days=C['min_days'], sigma='sem')
    if not net.stations:
        return site, 'rejected by clean(): %s' % net.rejected.get(s.name, '')
    s = net.stations[0]
    # month of each daily epoch, as GNSSNetwork.monthly() assigns it (see README, 'Known limitations')
    keys = np.array([_decyear_to_month(x) for x in s.t])
    months = list(s.monthly_t)

    def monthly_mean(path):
        if not path.exists():
            raise FileNotFoundError(path)
        t, u = read_loadpred(path)
        if len(t) < 100:
            raise ValueError('%s: fewer than 100 epochs' % path)
        ud = np.interp(s.t, t, u)
        return np.array([ud[keys == m].mean() for m in months])

    try:
        corr = {m: monthly_mean(config.ngl_dir('%s.cf/%s' % (m, site.lower()))) for m in C['load_models']}
        corr['cmcf'] = monthly_mean(config.ngl_dir('cmcf/%s' % site.lower())) if f['cmcf'] else None
    except (FileNotFoundError, ValueError) as e:
        return site, 'correction missing: %s' % e
    U0 = np.asarray(s.monthly_xyz['U'])
    load = sum(corr[m] for m in C['load_models'])
    U = U0 - load if corr['cmcf'] is None else U0 - corr['cmcf'] - load
    return site, dict(months=months, U=U, sU=np.asarray(s.monthly_sig['U']), U_cleaned=U0,
                      **{m: corr[m] for m in C['load_models']}, cmcf=corr['cmcf'])


def gnss_raw(frame, sites):
    """clean the NGL daily series of the stations (input/gnss/NGL/) -> {station: monthly series}"""
    t0 = time.time()
    steps = read_ngl_steps(config.ngl_dir('steps.txt'))
    with Pool(C['workers'], initializer=_init_worker, initargs=(steps,)) as p:
        res = dict(p.imap_unordered(_clean_one, [(s, frame) for s in sites], chunksize=8))
    out = {s: res[s] for s in sites if isinstance(res[s], dict)}
    bad = {s: r for s, r in res.items() if not isinstance(r, dict)}
    log('GNSS %s from the daily series: %d of %d stations, %.0f s' % (frame, len(out), len(sites), time.time() - t0))
    for s, r in list(bad.items())[:20]:
        log('  not processed: %s %s' % (s, r))
    return out


def write_series(path, frame, GN, S):
    """cleaned monthly series as NetCDF (station x month, metres) and as a long CSV (millimetres)"""
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(GN)
    var = ['dU', 'sigma_dU', 'dU_cleaned', 'ntal', 'ntol'] + (['cmcf'] if FRAMES[frame]['cmcf'] else [])
    key = dict(dU='U', sigma_dU='sU', dU_cleaned='U_cleaned', ntal='ntal', ntol='ntol', cmcf='cmcf')
    A = {v: np.full((len(names), len(ALLM)), np.nan) for v in var}
    for i, s in enumerate(names):
        j = [ALLM.index(m) for m in GN[s]['months']]
        for v in var:
            A[v][i, j] = GN[s][key[v]]
    desc = dict(dU='vertical displacement to compare with hydrological loading: dU_cleaned - ntal - ntol'
                   + (' - cmcf' if FRAMES[frame]['cmcf'] else ''),
                sigma_dU='uncertainty of the monthly mean: max(std, mean formal sigma) / sqrt(n days)',
                dU_cleaned='cleaned series (bias, rate and NGL offsets removed; seasonal terms and residual kept), '
                           'monthly mean of the daily values',
                ntal='ESMGFZ non-tidal atmospheric loading, CF, mean over the same days',
                ntol='ESMGFZ non-tidal oceanic loading, CF, mean over the same days',
                cmcf='NGL CM -> CF translation (up), mean over the same days')
    with nc.Dataset(path, 'w') as d:
        d.title = 'GNSS vertical monthly series, frame %s: %s' % (frame, FRAMES[frame]['label'])
        d.period = '2003-01 .. 2016-12; months with >= 10 daily solutions'
        d.createDimension('station', len(names)); d.createDimension('month', len(ALLM))
        vs = d.createVariable('station', str, ('station',)); vs[:] = np.array(names, dtype=object)
        vm = d.createVariable('month', str, ('month',)); vm[:] = np.array(ALLM, dtype=object)
        for c in ('lat', 'lon'):
            x = d.createVariable(c, 'f8', ('station',)); x[:] = S.loc[names, c].values
            x.units = 'degrees'; x.long_name = 'station coordinate used for the model (5 significant digits)'
        for v in var:
            x = d.createVariable(v, 'f8', ('station', 'month'), fill_value=np.nan)
            x[:] = A[v]; x.units = 'm'; x.long_name = desc[v]
    rows = []
    for i, s in enumerate(names):
        k = np.isfinite(A['dU'][i])
        df = pd.DataFrame({'station': s, 'month': np.array(ALLM)[k]})
        for v in var:
            df[v + '_mm'] = 1e3 * A[v][i, k]
        rows.append(df)
    pd.concat(rows).to_csv(path.with_suffix('.csv'), index=False, float_format='%.10g')


def read_series(path):
    """{station: dict(months, U, sU, ...)} from a NetCDF written by write_series"""
    out = {}
    with nc.Dataset(path) as d:
        names, months = list(d['station'][:]), list(d['month'][:])
        A = {v: np.ma.filled(d[v][:], np.nan) for v in d.variables if v not in ('station', 'month', 'lat', 'lon')}
    key = dict(dU='U', sigma_dU='sU', dU_cleaned='U_cleaned')
    for i, s in enumerate(names):
        k = np.isfinite(A['dU'][i])
        out[s] = dict(months=[months[j] for j in np.nonzero(k)[0]], **{key.get(v, v): A[v][i, k] for v in A})
    return out


def compare_series(new, ref):
    d, n_miss = [], 0
    for s, r in ref.items():
        if s not in new or list(new[s]['months']) != list(r['months']):
            n_miss += 1
            continue
        d.append(np.max(np.abs(np.asarray(new[s]['U']) - r['U'])))
    return '%d stations equal in months, %d differ in months or missing; max |dU diff| %.3g mm' % (
        len(d), n_miss, 1e3 * max(d) if d else np.nan)


def gnss(frame, S, source):
    """monthly series of the used stations in this frame"""
    proc = config.processed_dir(frame) / 'gnss_vertical_monthly.nc'
    sites = list(S.index[S.used])
    if source == 'processed':
        if not proc.exists():
            raise FileNotFoundError('%s missing: run with --gnss raw (needs input/gnss/NGL/)' % proc)
        GN = read_series(proc)
        log('GNSS %s: %d stations from %s' % (frame, len(GN), proc))
        return GN
    GN = gnss_raw(frame, sites)
    if proc.exists():
        target = config.output_dir(frame) / 'gnss_vertical_monthly_rebuilt.nc'
        log('rebuilt vs %s: %s' % (proc, compare_series(GN, read_series(proc))))
    else:
        target = proc
    write_series(target, frame, GN, S)
    log('written %s (and .csv)' % target)
    return GN


# ---------------------------------------------------------------------------------------------------------------
# 2. products: monthly anomalies on the land domain, spherical-harmonic coefficients
# ---------------------------------------------------------------------------------------------------------------
def _months_from_days(times, since, m0='2003-01', m1='2016-12'):
    import datetime as dt
    seen = {}
    for i, x in enumerate(times):
        d = since + dt.timedelta(days=float(x))
        m = '%d-%02d' % (d.year, d.month)
        if m0 <= m <= m1:
            seen.setdefault(m, i)                    # first field of a month
    return seen


def read_csr(land, path):
    """CSR RL06.3 mascons (0.25 deg, cm) -> 2 x 2 block means on the 0.5-degree grid [mm]"""
    import datetime as dt
    g = nc.Dataset(path)
    seen = _months_from_days(g['time'][:], dt.datetime(2002, 1, 1))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)  # valid_min/max attributes of lat/lon that netCDF4 cannot apply
        glat = np.array(g['lat'][:]); glon = np.array(g['lon'][:]); glon = np.where(glon > 180, glon - 360, glon)
    o = np.argsort(glon); oy = np.argsort(-glat)
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    for m, i in seen.items():
        w = np.ma.filled(g['lwe_thickness'][i], np.nan)[:, o][oy, :] * 10.0
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            R[ALLM.index(m)] = np.nanmean(w.reshape(360, 2, 720, 2), axis=(1, 3))[land]
    g.close()
    return R


def read_jpl(land, path):
    """JPL RL06.3M CRI mascons (0.5 deg, cm) [mm]"""
    import datetime as dt
    j = nc.Dataset(path)
    seen = _months_from_days(j['time'][:], dt.datetime(2002, 1, 1))
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    for m, i in seen.items():
        w = np.ma.filled(j['lwe_thickness'][i], np.nan) * 10.0
        R[ALLM.index(m)] = np.roll(w[::-1, :], 360, axis=1)[land]                   # 0..360 -> -180..180, N -> S
    j.close()
    return R


def read_glws(land, path):
    """GLWS 2.0 / 4.0 ensemble mean (0.5 deg, mm)"""
    R = np.full((len(ALLM), int(land.sum())), np.nan, np.float32)
    with nc.Dataset(path) as d:
        d.set_auto_mask(True)
        lat, lon = d['latitude'][:], d['longitude'][:]
        assert np.allclose(lat, LAT) and np.allclose(lon, LON)
        tt = d['time'][:]
        mo = ['%d-%02d' % (int(np.floor(v)), int(np.floor((v - np.floor(v)) * 12)) + 1) for v in tt]
        v = d['Total water storage anomalies']
        for i, m in enumerate(ALLM):
            if m in mo:
                R[i] = np.ma.filled(v[mo.index(m)], np.nan)[land]
    return R.astype(float)


def _gldas_month(args):
    """daily GLDAS-2.2 files of one month -> monthly mean (cells with >= 80 % of the days) -> 0.5 deg (cos-lat weighted
    mean of the 2 x 2 child cells with data) on the land domain"""
    import datetime as dt
    m, files, land = args
    if not files:
        return m, None
    y, mm = int(m[:4]), int(m[5:])
    nd = (dt.date(y + (mm == 12), mm % 12 + 1, 1) - dt.date(y, mm, 1)).days
    lat4 = -59.875 + 0.25 * np.arange(600)
    w4 = np.cos(np.deg2rad(lat4))[:, None] * np.ones((1, 1440))
    acc, cnt = np.zeros((600, 1440)), np.zeros((600, 1440))
    for f in files:
        with nc.Dataset(f) as d:
            x = np.ma.filled(d['TWS_tavg'][0], np.nan).astype(np.float64)
        x[x <= -9998] = np.nan
        ok = np.isfinite(x); acc[ok] += x[ok]; cnt[ok] += 1
    mean = np.where(cnt >= 0.8 * nd, acc / np.maximum(cnt, 1), np.nan)
    ok = np.isfinite(mean)
    num = np.where(ok, mean * w4, 0).reshape(300, 2, 720, 2).sum((1, 3))
    den = (ok * w4).reshape(300, 2, 720, 2).sum((1, 3))
    g05 = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)               # rows -59.75 .. 89.75
    full = np.full((360, 720), np.nan)
    full[:300] = g05[::-1]                                                          # rows 89.75 .. -59.75
    return m, full[land].astype(np.float32)


def read_gldas(land, folder, pattern):
    by = {}
    for f in sorted(Path(folder).glob(pattern)):
        d = f.name.split('.A')[1][:8]
        by.setdefault('%s-%s' % (d[:4], d[4:6]), []).append(str(f))
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    with Pool(C['workers']) as p:
        for m, o in p.imap(_gldas_month, [(m, by.get(m, []), land) for m in ALLM]):
            if o is not None:
                R[ALLM.index(m)] = o
    return R


UNITS = {'mm': 1.0, 'kg m-2': 1.0, 'kg/m2': 1.0, 'kg m**-2': 1.0, 'kg/m^2': 1.0, 'cm': 10.0, 'm': 1000.0}


def read_grid(land, path, var=None, units=None):
    """
    Any monthly TWS grid in one NetCDF file -> 0.5-degree land cells [mm]. Expected: a variable with dimensions
    (time, lat, lon) in any order of lat / lon values; a CF time axis ('days since ...', any calendar); units
    mm, kg m-2, cm or m (attribute 'units' or the option --tws-units); a regular grid of 0.5 deg, a finer grid that
    divides 0.5 deg (block means, e.g. 0.25 or 0.1 deg) or a coarser one that is a multiple of it (e.g. 1 deg).
    Absolute storage or anomalies: the study subtracts the product's own mean over its months either way. One field
    per month (the first field of a calendar month is used).
    """
    d = nc.Dataset(path)
    if var is None:
        cand = [v for v in d.variables.values() if v.ndim == 3]
        if len(cand) != 1:
            raise ValueError('%s: give the variable with --tws-var (3-D variables: %s)' % (path, [v.name for v in cand]))
        v = cand[0]
    else:
        v = d[var]
    dims = list(v.dimensions)
    tname = [n for n in dims if n.lower() in ('time', 't', 'month')]
    yname = [n for n in dims if n.lower() in ('lat', 'latitude', 'y')]
    xname = [n for n in dims if n.lower() in ('lon', 'longitude', 'x')]
    if not (tname and yname and xname):
        raise ValueError('%s: dimensions %s, expected time, lat, lon' % (path, dims))
    tv = d[tname[0]]
    dates = nc.num2date(tv[:], tv.units, getattr(tv, 'calendar', 'standard'))
    lat = np.asarray(d[yname[0]][:], float); lon = np.asarray(d[xname[0]][:], float)
    lon = (lon + 180) % 360 - 180
    u = (units or getattr(v, 'units', '')).strip()
    if u not in UNITS:
        raise ValueError("%s: units '%s' unknown; give --tws-units (one of %s)" % (path, u, ', '.join(UNITS)))
    res = abs(lat[1] - lat[0])
    order = [dims.index(tname[0]), dims.index(yname[0]), dims.index(xname[0])]
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    seen = {}
    for i, t in enumerate(dates):
        m = '%04d-%02d' % (t.year, t.month)
        if '2003-01' <= m <= '2016-12' and m not in seen:
            seen[m] = i
    for m, i in seen.items():
        sl = [slice(None)] * 3; sl[order[0]] = i
        w = np.ma.filled(np.ma.masked_invalid(v[tuple(sl)]).astype(float), np.nan)
        if order[1] > order[2]:
            w = w.T
        R[ALLM.index(m)] = _to_half_degree(w, lat, lon, res)[land] * UNITS[u]
    d.close()
    return R


def _to_half_degree(w, lat, lon, res):
    """(lat, lon) grid -> 360 x 720 grid of the study (NaN where the input has no cell)"""
    out = np.full((360, 720), np.nan)
    if abs(res - 0.5) < 1e-6:
        rows = np.rint((89.75 - lat) / 0.5).astype(int); cols = np.rint((lon + 179.75) / 0.5).astype(int)
        out[np.ix_(rows, cols)] = w
    elif res < 0.5 and abs(0.5 / res - round(0.5 / res)) < 1e-6:
        R, Cc = np.meshgrid(np.floor((90 - lat) / 0.5).astype(int), np.floor((lon + 180) / 0.5).astype(int) % 720,
                            indexing='ij')
        ok = np.isfinite(w) & (R >= 0) & (R < 360)
        acc = np.zeros((360, 720)); cnt = np.zeros((360, 720))
        np.add.at(acc, (R[ok], Cc[ok]), w[ok]); np.add.at(cnt, (R[ok], Cc[ok]), 1)
        out = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)               # mean of the child cells with data
    elif res > 0.5 and abs(res / 0.5 - round(res / 0.5)) < 1e-6:
        k = int(round(res / 0.5))
        r1 = np.rint((90 - (lat + res / 2)) / 0.5).astype(int); c1 = np.rint((lon - res / 2 + 180) / 0.5).astype(int)
        for a, i in enumerate(r1):
            if i < 0 or i + k > 360:
                continue
            for b, j in enumerate(c1):
                out[i:i + k, (j + np.arange(k)) % 720] = w[a, b]
    else:
        raise ValueError('grid spacing %g deg: only 0.5 deg, divisors or multiples of it' % res)
    return out


def product_source(name, files=None):
    """(reader kind, path, options) of a product; files = {name: path} from --tws / --product-file"""
    if files and name in files:
        p = PRODUCTS.get(name, dict(kind='grid'))
        return (p['kind'] if name in PRODUCTS else 'grid'), Path(files[name]).expanduser()
    if name not in PRODUCTS:
        raise ValueError('product %s: not one of %s; give its file with --tws %s=FILE' % (name, list(PRODUCTS), name))
    p = PRODUCTS[name]
    path = config.data_path(p['path'])
    if 'pattern' in p:
        found = sorted(path.glob(p['pattern'])) if path.is_dir() else []
        if not found:
            raise FileNotFoundError('%s: no %s in %s; see download.py' % (name, p['pattern'], path))
        path = path if p['kind'] == 'gldas' else found[-1]
    elif not path.exists():
        raise FileNotFoundError('%s: %s missing; see download.py' % (name, path))
    return p['kind'], path


def prepare(name, land, months, files=None, var=None, units=None):
    """monthly anomaly grid on the land domain [mm] and SH coefficients of the EWH [m], cached per input file"""
    kind, path = product_source(name, files)
    out = config.cache_dir() / ('shc_%s.npy' % name)
    pat = PRODUCTS.get(name, {}).get('pattern', '*')
    fl = sorted(path.glob(pat)) if path.is_dir() else [path]
    tag = '%s\n%s %d files, %d bytes, var %s, units %s\n' % (kind, path, len(fl), sum(f.stat().st_size for f in fl),
                                                              var, units)
    tagf = out.with_suffix('.source.txt')
    if out.exists() and tagf.exists() and tagf.read_text() == tag:
        log('%s: cached (%s)' % (name, out))
        return
    t0 = time.time()
    if kind == 'gldas':
        R = read_gldas(land, path, PRODUCTS[name]['pattern'])
    elif kind == 'grid':
        R = read_grid(land, path, var, units)
    else:
        R = dict(csr=read_csr, jpl=read_jpl, glws=read_glws)[kind](land, path)
    R = R[[ALLM.index(m) for m in months]]
    okm = np.isfinite(R).mean(1) > 0.5                                  # months with data
    if okm.sum() < 24:
        raise ValueError('%s: only %d of the %d common months have data' % (name, okm.sum(), len(months)))
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)                  # cells without any data stay NaN -> 0
        A = R - np.nanmean(R[okm], axis=0)[None, :]                     # anomaly about the product's own mean
    nanfrac = float(np.mean(~np.isfinite(A[okm])))
    A[~np.isfinite(A)] = 0.0
    A[~okm] = np.nan
    np.savez_compressed(config.cache_dir() / ('field_%s.npz' % name), months=np.array(months),
                        anom_mm=A.astype(np.float32), month_ok=okm)
    plat, plon = GeoMathKit.global_equal_distance_grid(0.5)
    assert np.allclose(plat, LAT[::-1]) and np.allclose(plon, LON)
    G = np.zeros((len(months), 360, 720))
    G[:, land] = np.nan_to_num(A) / 1000.0
    G = G[:, ::-1, :]                                                    # EWH [m], south -> north
    sb = SphericalHarmonicBasis(C['lmax'])
    n = C['sh_chunk']
    Cf = np.concatenate([sb.load_parameters(G[i:i + n], plat, plon) for i in range(0, len(months), n)], axis=0)
    Cf[~okm] = np.nan
    np.save(out, Cf)
    tagf.write_text(tag)
    log('%s: %d months (%d with data), land cells without data set to 0: %.4f, |anomaly| rms %.1f mm, %.0f s'
        % (name, len(months), okm.sum(), nanfrac, np.sqrt(np.nanmean(A[okm] ** 2)), time.time() - t0))


# ---------------------------------------------------------------------------------------------------------------
# 3. forward model
# ---------------------------------------------------------------------------------------------------------------
def forward(names, lat, lon):
    """vertical displacement [m] (n_station, n_month) per product at the stations"""
    lln = LoveNumber().config(lmax=C['lmax'], method=getattr(LLN_Data, C['love'])).get_Love_number()
    lln.convert(target=Frame.CF)
    sb = SphericalHarmonicBasis(C['lmax'])
    coef = {n: np.load(config.cache_dir() / ('shc_%s.npy' % n)) for n in names}
    U = {n: np.empty((len(lat), coef[n].shape[0])) for n in names}
    for i in range(0, len(lat), C['station_block']):
        b = slice(i, i + C['station_block'])
        A = sb.design_matrix(lat[b], lon[b], lln, Displacement.Vertical)
        for n in names:
            U[n][b] = A @ coef[n].T
    return U


# ---------------------------------------------------------------------------------------------------------------
# 4. scores
# ---------------------------------------------------------------------------------------------------------------
def tyear(months):
    return np.array([int(m[:4]) + (int(m[5:]) - 0.5) / 12 for m in months])


def harm(t, y):
    tc = t - 2010.0
    X = np.column_stack([np.ones_like(tc), tc, np.cos(2 * np.pi * tc), np.sin(2 * np.pi * tc),
                         np.cos(4 * np.pi * tc), np.sin(4 * np.pi * tc)])
    return np.linalg.lstsq(X, y, rcond=None)[0], X


def clim(months, y):
    mm = np.array([int(m[5:]) for m in months]); c = np.full(12, np.nan); n = np.zeros(12, int)
    for k in range(12):
        s = mm == k + 1; n[k] = s.sum()
        if s.any():
            c[k] = y[s].mean()
    return c, n, mm


def corr(a, b):
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def station_scores(months, g, sg, m):
    t = tyear(months); g = g - g.mean(); m = m - m.mean()
    bg, X = harm(t, g); bm, _ = harm(t, m)
    gd = g - X[:, :2] @ bg[:2]; md = m - X[:, :2] @ bm[:2]
    cg, n, mm = clim(months, gd); cm, _, _ = clim(months, md)
    cg0, _, _ = clim(months, g); cm0, _, _ = clim(months, m)
    w = 1 / np.maximum(sg, 1e-5) ** 2
    wr = lambda x: np.sqrt(np.sum(w * x ** 2) / np.sum(w))
    ag, am = np.hypot(bg[2], bg[3]), np.hypot(bm[2], bm[3])
    pg = np.arctan2(bg[3], bg[2]) / (2 * np.pi) * 365.25; pm = np.arctan2(bm[3], bm[2]) / (2 * np.pi) * 365.25
    dph = (pm - pg + 182.625) % 365.25 - 182.625
    return dict(n_months=len(months), r_total=corr(gd, md), r_seas=corr(cg, cm) if (n >= 2).all() else np.nan,
                r_nonseas=corr(g - cg0[mm - 1], m - cm0[mm - 1]), r_nonseas_dt=corr(gd - cg[mm - 1], md - cm[mm - 1]),
                amp_gnss_mm=1e3 * ag, amp_model_mm=1e3 * am, amp_ratio=am / ag if ag > 0 else np.nan,
                phase_diff_d=dph, wrms_gnss_mm=1e3 * wr(gd), wrms_red=1 - wr(gd - md) / wr(gd),
                wrms_red_nodt=1 - wr((g - X[:, :2] @ bg[:2]) - m) / wr(gd), trend_model_mm_yr=1e3 * bm[1])


def _rng(*key):
    """one bootstrap stream per (group, product, score): results do not depend on which products are run together"""
    return np.random.default_rng([C['seed'], zlib.crc32('|'.join(key).encode())])


def boot_median(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0
    b = np.median(x[rng.integers(0, len(x), (C['n_boot'], len(x)))], axis=1)
    return float(np.median(x)), float(np.percentile(b, 5)), float(np.percentile(b, 95)), len(x)


def boot_share(win, rng):
    win = np.asarray(win, float); win = win[np.isfinite(win)]
    if len(win) == 0:
        return np.nan, np.nan, np.nan, 0
    b = win[rng.integers(0, len(win), (C['n_boot'], len(win)))].mean(1)
    return float(win.mean()), float(np.percentile(b, 5)), float(np.percentile(b, 95)), len(win)


def beats(a, b, s):
    """1 where a is better than b on score s, 0 where worse, 0.5 for ties, NaN if missing"""
    if s == 'amp_ratio':
        a, b = -np.abs(np.log(a)), -np.abs(np.log(b))
    elif s == 'phase_diff_d':
        a, b = -np.abs(a), -np.abs(b)
    return np.where(np.isfinite(a) & np.isfinite(b), np.where(a > b, 1.0, np.where(a < b, 0.0, 0.5)), np.nan)


def per_station(names, S, GN, U, months):
    """long table of the scores: one row per station and product"""
    JM = set(months)
    rows = []
    for k, st in enumerate(S.index):
        g = GN[st]; gm = np.array(g['months'])
        for p in names:
            u = U[p][k]
            sel = [i for i, m in enumerate(gm) if m in JM and np.isfinite(u[months.index(m)])]
            mo = gm[sel]; mv = np.array([u[months.index(m)] for m in mo])
            r = station_scores(list(mo), np.asarray(g['U'])[sel], np.asarray(g['sU'])[sel], mv)
            r.update(station=st, product=p); rows.append(r)
    return pd.DataFrame(rows)


def groups_of(PS):
    PS = PS.copy()
    PS['japan'] = PS.lat.between(24, 46) & PS.lon.between(122, 146)     # GEONET: groups without Japan too
    g = [('global', PS), ('global_exJapan', PS[~PS.japan])]
    g += [('region:' + k, v) for k, v in PS.groupby('wmo_region')]
    g += [('distance:' + k, v) for k, v in PS.groupby('boundary_class')]
    g += [('distance_exJapan:' + k, v) for k, v in PS[~PS.japan].groupby('boundary_class')]
    g += [('footprint:' + ('in' if k else 'out(>40km)'), v) for k, v in PS.groupby('in_footprint')]
    return g


def summarise(PS, names, baseline=None):
    """medians with station-bootstrap 5-95 % ranges per group, product and score (+ share of stations where a product
    beats the baseline product)"""
    summ = []
    for grp, sub in groups_of(PS):
        for p in names:
            for s in SCORES:
                med, lo, hi, n = boot_median(sub['%s|%s' % (p, s)], _rng(grp, p, s))
                d = dict(group=grp, product=p, score=s, n=n, median=med, p05=lo, p95=hi)
                if s == 'amp_ratio':
                    d['median_abs_log'] = float(np.nanmedian(np.abs(np.log(sub['%s|%s' % (p, s)]))))
                if s == 'phase_diff_d':
                    d['median_abs'] = float(np.nanmedian(np.abs(sub['%s|%s' % (p, s)])))
                summ.append(d)
            if baseline and p != baseline:
                for s in SCORES:
                    sh, lo, hi, n = boot_share(beats(sub['%s|%s' % (p, s)].values, sub['%s|%s' % (baseline, s)].values,
                                                     s), _rng(grp, p, 'share_beats_' + baseline, s))
                    summ.append(dict(group=grp, product=p, score='share_beats_%s:%s' % (baseline, s), n=n, median=sh,
                                     p05=lo, p95=hi))
    return pd.DataFrame(summ)


def headline(SM, names, baseline=None):
    """short table: global and regional medians of the main scores"""
    G = [('global', 'global'), ('global_exJapan', 'global ex-Japan'), ('region:Europe', 'Europe'),
         ('region:North America, Central America and the Caribbean', 'N America'), ('region:Asia', 'Asia'),
         ('region:South America', 'S America'), ('region:Africa', 'Africa'),
         ('region:South-West Pacific', 'SW Pacific')]
    rows = []
    for grp, gl in G:
        for p in names:
            x = SM[(SM.group == grp) & (SM['product'] == p)].set_index('score')
            if x.empty:
                continue
            d = dict(group=gl, product=p, n=int(x.loc['wrms_red', 'n']))
            for s in ('wrms_red', 'r_seas', 'r_nonseas_dt', 'amp_ratio', 'phase_diff_d', 'r_total'):
                d[s] = x.loc[s, 'median']; d[s + '_p05'] = x.loc[s, 'p05']; d[s + '_p95'] = x.loc[s, 'p95']
            d['phase_median_abs_d'] = x.loc['phase_diff_d', 'median_abs']
            if baseline and p != baseline:
                d['share_beats_%s' % baseline] = x.loc['share_beats_%s:wrms_red' % baseline, 'median']
            rows.append(d)
    return pd.DataFrame(rows)


def score_frame(frame, names, S, GN, U, months, baseline=None):
    out = config.output_dir(frame)
    L = per_station(names, S, GN, U, months)
    W = L.pivot(index='station', columns='product', values=SCORES + ['n_months', 'amp_model_mm', 'trend_model_mm_yr'])
    W.columns = ['%s|%s' % (p, s) for s, p in W.columns]
    full = max(names, key=lambda p: L[L['product'] == p].n_months.sum())   # product on the most months
    geo = ['basin_unit', 'dist_unit_boundary_km', 'in_footprint', 'wmo_region', 'boundary_class']
    base = S[['lat', 'lon'] + geo].join(L[L['product'] == full].set_index('station')[['n_months', 'amp_gnss_mm',
                                                                                     'wrms_gnss_mm']])
    PS = base.join(W[sorted(W.columns, key=lambda c: (names.index(c.split('|')[0]), c))])
    PS.index.name = 'station'
    PS.to_csv(out / 'gnss_scores_per_station.csv', float_format='%.12g')
    SM = summarise(PS, names, baseline)
    SM.to_csv(out / 'gnss_summary.csv', index=False, float_format='%.6g')
    H = headline(SM, names, baseline)
    H.to_csv(out / 'headline.csv', index=False, float_format='%.4g')
    return PS, SM, H


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--data-root', help='PyGLoad data root (default: installation/pygload_paths.cfg or PYGLOAD_DATA)')
    ap.add_argument('--frames', default=','.join(C['frames']), help='igs20cf (main), igs14 (sensitivity)')
    ap.add_argument('--products', default=None,
                    help='comma-separated public products: %s (default CSR,JPL; none with --tws only)'
                         % ','.join(PRODUCTS))
    ap.add_argument('--tws', action='append', default=[], metavar='NAME=FILE',
                    help='your own monthly TWS grid (NetCDF, see read_grid); repeatable')
    ap.add_argument('--tws-var', default=None, help='variable of the --tws files (default: the only 3-D variable)')
    ap.add_argument('--tws-units', default=None, help='units of the --tws files if the attribute is missing or wrong')
    ap.add_argument('--product-file', action='append', default=[], metavar='NAME=FILE',
                    help='a public product kept elsewhere, e.g. JPL=/data/GRCTellus...nc (GLDAS22: the folder)')
    ap.add_argument('--gnss', choices=('processed', 'raw'), default=C['gnss'],
                    help="'processed': the cleaned series of processed/; 'raw': clean the NGL daily series again")
    ap.add_argument('--baseline', default=None, help='also give the share of stations where each product beats this one')
    ap.add_argument('--workers', type=int, default=C['workers'])
    a = ap.parse_args()
    if a.data_root:
        config.set_data_root(a.data_root)
    C['workers'] = a.workers
    t00 = time.time()
    files = dict(x.split('=', 1) for x in a.tws + a.product_file)
    own = [x.split('=', 1)[0] for x in a.tws]
    names = [p for p in (a.products.split(',') if a.products is not None else ([] if own else C['products'])) if p]
    names += [n for n in own if n not in names]
    frames = [f for f in a.frames.split(',') if f]
    for f in frames:
        if f not in FRAMES:
            raise SystemExit('unknown frame %s (known: %s)' % (f, ', '.join(FRAMES)))
    if a.baseline and a.baseline not in names:
        raise SystemExit('--baseline %s must be one of the products run (%s)' % (a.baseline, names))
    land, months = domain()
    S = stations()
    U_st = S[S.used]
    log('data root %s' % config.FileTool.get_data_dir())
    log('%d stations used, %d common months, %d land cells, products %s, frames %s'
        % (len(U_st), len(months), land.sum(), names, frames))
    GN = {f: gnss(f, S, a.gnss) for f in frames}
    for n in names:
        prepare(n, land, months, files, a.tws_var if n in own else None, a.tws_units if n in own else None)
    t0 = time.time()
    sta = [s for s in U_st.index if all(s in GN[f] for f in frames)]
    Ssub = U_st.loc[sta]
    U = forward(names, Ssub.lat.values, Ssub.lon.values)
    log('forward model: %d stations x %d products, %.0f s' % (len(Ssub), len(names), time.time() - t0))
    for f in frames:
        PS, SM, H = score_frame(f, names, Ssub, GN[f], U, months, a.baseline)
        log('frame %s (%s): scores in %s' % (f, FRAMES[f]['label'], config.output_dir(f)))
        x = H[H.group == 'global'][['product', 'n', 'wrms_red', 'r_seas', 'r_nonseas_dt', 'amp_ratio', 'phase_diff_d',
                                    'phase_median_abs_d']]
        print(x.to_string(index=False, float_format=lambda v: '%.3f' % v))
    log('total %.0f s' % (time.time() - t00))


if __name__ == '__main__':
    main()

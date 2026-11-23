"""
Checks of a GNSS_GRACE run (after run.py), written to <data root>/output/GNSS_GRACE/<name>/verify.txt:

  1. station check: annual vertical amplitude and phase of each station against the displacement predicted from GRACE
     (CSR mascons, 0.5-degree disks, Love numbers in the CF frame); a ratio near 1 and no lag means the GNSS seasonal
     signal is continental water
  2. Fu et al. (2015), if their digitised data are in <data root>/input/stations/: uplift 1 Apr -> 1 Oct at shared
     stations, and the water added 1 Oct -> 1 Apr over Washington and Oregon (GNSS inversion and GRACE)

    python studies/GNSS_GRACE/verify.py
"""
import numpy as np

import PyGWater as pgw
from PyGWater.config import FileTool
from PyGWater.preprocessing import GNSSDataset
from PyGWater.observations.grace import Mascon
from PyGWater.mesh.land import land_fraction
from PyGWater.utils.constants import EarthConstant

from config import CONFIG as C

OCT1, APR1 = 273 / 365.25, 90 / 365.25


def decyear(months):
    return np.array([int(m[:4]) + (int(m[5:7]) - 0.5) / 12 for m in months])


def harmonics(t):
    return np.c_[np.ones_like(t), t - 2014, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t), np.sin(4 * np.pi * t),
                 np.cos(4 * np.pi * t)]


def seasonal(c, t0):
    """seasonal part (annual + semi-annual) of fitted coefficients c at the phase t0 of the year"""
    return (c[2] * np.sin(2 * np.pi * t0) + c[3] * np.cos(2 * np.pi * t0) + c[4] * np.sin(4 * np.pi * t0)
            + c[5] * np.cos(4 * np.pi * t0))


def fit(t, y):
    ok = np.isfinite(y)
    return np.linalg.lstsq(harmonics(t[ok]), y[ok], rcond=None)[0] if ok.sum() >= 36 else None


def main():
    out = FileTool.get_output_dir('GNSS_GRACE/%s' % C['name'])
    ds = GNSSDataset.load(out / 'dataset.npz')
    inv = np.load(out / 'inversion.npz')
    lines = []

    def say(s):
        print(s)
        lines.append(s)

    lat_min, lat_max, lon_min, lon_max = C['box']
    g = Mascon.read(FileTool.get_input_dir(C['grace_file']), box=(lat_min - 11, lat_max + 11, lon_min - 15,
                                                                   lon_max + 15)).anomalies(C['reference'])
    la, lo = g.lat.reshape(-1, 2).mean(1), g.lon.reshape(-1, 2).mean(1)
    E = g.ewh.reshape(len(g.months), la.size, 2, lo.size, 2).mean((2, 4)).reshape(len(g.months), -1)
    cells_g = pgw.GreenFunctionBasis.from_grid(la, lo, 0.5, kernel='disk')
    lln = pgw.LoveNumber().config(lmax=C['lmax'], method=pgw.LLN_Data.Wang).get_Love_number()
    lln.convert(pgw.Frame.CF)

    # 1. station check
    pred = E @ cells_g.design_matrix(ds.lat, ds.lon, lln, pgw.Displacement.Vertical).T
    gi = {m: i for i, m in enumerate(g.months)}
    common = [m for m in ds.months if m in gi]
    t = decyear(common)
    obs = ds.values['U'][[ds.months.index(m) for m in common]]
    prd = pred[[gi[m] for m in common]]
    amp_o, amp_p, lag = [], [], []
    for k in range(len(ds.names)):
        co, cp = fit(t, obs[:, k]), fit(t, prd[:, k])
        if co is None:
            continue
        amp_o.append(np.hypot(co[2], co[3]))
        amp_p.append(np.hypot(cp[2], cp[3]))
        lag.append(np.angle(np.exp(1j * (np.arctan2(co[3], co[2]) - np.arctan2(cp[3], cp[2])))) / (2 * np.pi) * 365.25)
    amp_o, amp_p, lag = map(np.array, (amp_o, amp_p, lag))
    say('Station check (%d stations): annual vertical amplitude GNSS %.1f mm, GRACE-predicted %.1f mm (medians); '
        'ratio %.2f; lag %.0f days (IQR %.0f to %.0f)' % (len(amp_o), 1e3 * np.median(amp_o), 1e3 * np.median(amp_p),
                                                          np.median(amp_o / amp_p), *np.percentile(lag, [50, 25, 75])))

    # 2. Fu et al. (2015)
    stations = FileTool.get_input_dir('stations')
    if (stations / 'fu2015_seasonal_uplift_digitized.txt').exists():
        fu = np.loadtxt(stations / 'fu2015_seasonal_uplift_digitized.txt')
        mine, theirs = [], []
        for la_, lo_, u in fu:
            d = np.hypot(ds.lat - la_, (ds.lon - lo_) * np.cos(np.deg2rad(la_))) * 111
            k = int(np.argmin(d))
            c = fit(decyear(ds.months), ds.values['U'][:, k]) if d[k] < 3 else None
            if c is not None:
                mine.append(1e3 * (seasonal(c, OCT1) - seasonal(c, APR1)))
                theirs.append(u)
        mine, theirs = np.array(mine), np.array(theirs)
        say('Fu et al. (2015) Fig. 1, uplift 1 Apr -> 1 Oct at %d shared stations: r = %.2f, median %.1f mm here, '
            '%.1f mm digitised (colour bar capped at 12 mm)' % (len(mine), np.corrcoef(mine, theirs)[0, 1],
                                                                 np.median(mine), np.median(theirs)))
        # water added 1 Oct -> 1 Apr over Washington and Oregon land
        r, b = C['resolution'], C['buffer']
        lat1 = np.arange(lat_min - b + r / 2, lat_max + b, r)
        lon1 = np.arange(lon_min - b + r / 2, lon_max + b, r)
        land = land_fraction(lat1, lon1, r) >= C['min_land'] if C['min_land'] else None
        opts = dict(subgrid=r / 4) if C['kernel'] == 'cluster_disk' else {}
        cells = pgw.GreenFunctionBasis.from_grid(lat1, lon1, r, mask=land, kernel=C['kernel'], **opts)
        months = [str(m) for m in inv['months']]
        cx = np.linalg.lstsq(harmonics(decyear(months)), inv['x'], rcond=None)[0]
        dx = seasonal(cx, APR1) - seasonal(cx, OCT1)
        glat, glon = np.arange(42.125, 49, 0.25), np.arange(-125.875, -116.5, 0.25)
        lo2, la2 = np.meshgrid(glon, glat)
        m = land_fraction(glat, glon, 0.25) >= 0.5
        area = np.cos(np.deg2rad(la2[m])) * np.deg2rad(0.25) ** 2 * EarthConstant.radiusm ** 2
        gnss = cells.load_synthesis(dx[None, :], la2[m], lo2[m])[0]
        cg = np.linalg.lstsq(harmonics(decyear(g.months)), g.sample(la2[m], lo2[m]), rcond=None)[0]
        grace = seasonal(cg, APR1) - seasonal(cg, OCT1)
        gt = lambda e: np.sum(e * area) * EarthConstant.rhow / 1e12
        say('Water added 1 Oct -> 1 Apr, Washington and Oregon land: GNSS %.0f Gt (mean %.0f mm), GRACE %.0f Gt '
            '(mean %.0f mm); Fu et al. (2015) reproduced: about 105 Gt' % (gt(gnss), 1e3 * gnss.mean(), gt(grace),
                                                                             1e3 * grace.mean()))
    (out / 'verify.txt').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()

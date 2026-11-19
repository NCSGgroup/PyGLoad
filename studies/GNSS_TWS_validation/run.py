"""
GNSS-TWS validation: how well does the elastic loading predicted from a global terrestrial water storage (TWS)
product match the monthly vertical GNSS displacements? Settings in config.py.

  1. GNSS        : cleaned monthly vertical series of the stations used (trend and offsets removed, ESMGFZ NTAL + NTOL
                   subtracted), from the reference files or rebuilt from the NGL daily series     -> gnss_<source>.npz
  2. products    : each TWS product as monthly 0.5-degree anomalies on a common land domain (ocean, Greenland and
                   Antarctica = 0), on the common months; spherical-harmonic analysis (degree 360)  -> cache
  3. forward     : vertical displacement at the stations, SH design matrix with PREM load Love numbers in the CF frame
                                                                                                    -> cache
  4. scores      : per station and product (correlations, annual amplitude ratio and phase, WRMS reduction), medians
                   with station-bootstrap ranges, globally, per region and by distance to basin-unit borders
                                                                                                    -> output

    python studies/GNSS_TWS_validation/download.py     # once, needs internet
    python studies/GNSS_TWS_validation/run.py
    python studies/GNSS_TWS_validation/verify.py       # compare with the reference scores

Scores (both series demeaned on the common months, then mean + trend of a mean + trend + annual + semi-annual LS fit
removed, because the GNSS trend was removed when cleaning):
    r_total        Pearson r
    r_seas         r of the 12 calendar-month climatologies (each calendar month needs >= 2 values)
    r_nonseas      r of the deseasonalised series, model trend kept;  r_nonseas_dt: model trend removed too
    amp_ratio      annual amplitude model / GNSS;  phase_diff_d: day of the annual maximum, model - GNSS
    wrms_red       1 - WRMS(GNSS - model) / WRMS(GNSS), weights 1 / sigma^2 of the GNSS monthly means
                   (wrms_red_nodt: model trend kept)
"""
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import netCDF4 as nc
import pandas as pd

from pygwater import LoveNumber, LLN_Data, Frame, Displacement, SphericalHarmonicBasis
from pygwater.preprocessing import GNSSNetwork, StationSeries, read_tenv3, read_ngl_steps
from pygwater.preprocessing.series import _decyear_to_month
from pygwater.utils.geomath import GeoMathKit

from config import CONFIG as C, PRODUCTS, reference_dir, product_path as _configured_path, ngl_dir, cache_dir, output_dir

ALLM = ['%d-%02d' % (y, m) for y in range(2003, 2017) for m in range(1, 13)]   # 2003-01 .. 2016-12
LAT = 89.75 - 0.5 * np.arange(360)                                              # grid of the study: north -> south
LON = -179.75 + 0.5 * np.arange(720)
SCORES = ['r_total', 'r_seas', 'r_nonseas', 'r_nonseas_dt', 'amp_ratio', 'phase_diff_d', 'wrms_red', 'wrms_red_nodt']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# ---------------------------------------------------------------------------------------------------------------
# reference files
# ---------------------------------------------------------------------------------------------------------------
def domain():
    """land mask (360 x 720, north -> south) and the common months"""
    with nc.Dataset(reference_dir('domain/land_domain_0.5deg.nc')) as f:
        land = f['land'][:].astype(bool)
    months = [l.strip() for l in open(reference_dir('domain/common_months.txt')) if l.strip() and l[0] != '#']
    return land, months


def stations():
    return pd.read_csv(reference_dir('stations/stations.csv'), index_col=0)


# ---------------------------------------------------------------------------------------------------------------
# 1. GNSS
# ---------------------------------------------------------------------------------------------------------------
def gnss_reference():
    """{station: dict(months, U, sU)} from the cleaned monthly series of the reference files"""
    with nc.Dataset(reference_dir('gnss_monthly/gnss_vertical_monthly_2003_2016.nc')) as f:
        names, months = list(f['station'][:]), list(f['month'][:])
        U, sU = np.ma.filled(f['dU'][:], np.nan), np.ma.filled(f['sigma_dU'][:], np.nan)
    out = {}
    for i, s in enumerate(names):
        k = np.isfinite(U[i])
        out[s] = dict(months=[months[j] for j in np.nonzero(k)[0]], U=U[i, k], sU=sU[i, k])
    return out


def _read_loadpred(path):
    """NGL loadpredictions/<model>.cf/<site> file: decimal year and up [m]"""
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


def _clean_one(args):
    """daily NGL series of one station -> cleaned monthly vertical series minus NTAL + NTOL"""
    site, variant = args
    try:
        s = read_tenv3(str(ngl_dir('tenv3_IGS14/%s.tenv3' % site)))
    except Exception as e:
        return site, 'read error: %s' % e
    k = (s.t >= C['period'][0]) & (s.t < C['period'][1])
    if k.sum() < C['min_daily']:
        return site, 'fewer than %d daily solutions' % C['min_daily']
    s = StationSeries(s.name, s.lat, s.lon, s.t[k], {c: v[k] for c, v in s.xyz.items()},
                      {c: v[k] for c, v in s.sig.items()})
    if s.t.max() - s.t.min() < C['min_span']:
        return site, 'span below %g years' % C['min_span']
    net = GNSSNetwork([s])
    if variant == 'rebuild':
        # the settings of the reference series: trajectory model with NGL offsets and outlier rejection only; the
        # additional data selection and quality control of GNSSNetwork.clean are switched off
        net.clean(steps=_STEPS, remove='trend+steps', outlier=C['outlier'], min_years=0.0, min_completeness=0.0,
                  max_sigma=np.inf, max_residual=np.inf, detect_jumps=False,
                  segment_options=dict(min_segment=0.0, min_between=0.0))
        net.monthly(min_days=C['min_days'], sigma='sem')
    else:
        # GNSSNetwork.clean with its default data selection, offset detection and quality control
        net.clean(steps=_STEPS, remove='trend+steps', outlier=C['outlier'])
        net.monthly(min_days=C['min_days'])
    if not net.stations:
        return site, 'rejected by clean(): %s' % net.rejected.get(s.name, '')
    s = net.stations[0]
    load = 0.0
    for model in C['load_models']:
        f = ngl_dir('loadpred_cf/%s/%s' % (model, site.lower()))
        if not f.exists():
            return site, 'no %s file' % model
        tl, ul = _read_loadpred(f)
        if len(tl) < 100:
            return site, '%s file short' % model
        load = load + np.interp(s.t, tl, ul)
    keys = np.array([_decyear_to_month(x) for x in s.t])
    L = np.array([load[keys == m].mean() for m in s.monthly_t])        # mean over the days of each monthly mean
    return site, dict(months=list(s.monthly_t), U=np.asarray(s.monthly_xyz['U']) - L,
                      sU=np.asarray(s.monthly_sig['U']), L=L)


def _init_worker(steps):
    global _STEPS
    _STEPS = steps


def gnss_rebuild(variant, sites):
    """rebuild the cleaned monthly series of the stations from the NGL daily data (download.py --gnss)"""
    steps = read_ngl_steps(ngl_dir('steps.txt'))
    with Pool(C['workers'], initializer=_init_worker, initargs=(steps,)) as p:
        res = dict(p.imap_unordered(_clean_one, [(s, variant) for s in sites], chunksize=8))
    out = {s: r for s, r in res.items() if isinstance(r, dict)}
    bad = {s: r for s, r in res.items() if not isinstance(r, dict)}
    log('GNSS %s: %d of %d stations' % (variant, len(out), len(sites)))
    if bad:
        reasons = pd.Series(bad).str.replace(r':.*', '', regex=True).value_counts()
        log('  not rebuilt: %s' % reasons.to_dict())
    return out


def compare_gnss(new, ref):
    """differences between rebuilt and reference series, on the station-months both have"""
    d, extra, missing = [], 0, 0
    for s, r in ref.items():
        if s not in new:
            missing += len(r['months'])
            continue
        a = dict(zip(r['months'], r['U'])); b = dict(zip(new[s]['months'], new[s]['U']))
        common = sorted(set(a) & set(b))
        d += [b[m] - a[m] for m in common]
        missing += len(set(a) - set(b)); extra += len(set(b) - set(a))
    d = np.abs(np.array(d))
    return 'station-months in both %d (max |diff| %.3g mm, median %.3g mm), only in reference %d, only rebuilt %d' % (
        len(d), 1e3 * d.max() if len(d) else np.nan, 1e3 * np.median(d) if len(d) else np.nan, missing, extra)


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
    import datetime as dt
    g = nc.Dataset(path)
    seen = _months_from_days(g['time'][:], dt.datetime(2002, 1, 1))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)  # valid_min/max attributes of lat/lon that netCDF4 cannot apply
        glat = np.array(g['lat'][:]); glon = np.array(g['lon'][:]); glon = np.where(glon > 180, glon - 360, glon)
    o = np.argsort(glon); oy = np.argsort(-glat)
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    for m, i in seen.items():
        w = np.ma.filled(g['lwe_thickness'][i], np.nan)[:, o][oy, :] * 10.0          # cm -> mm, 0.25 deg
        with np.errstate(invalid='ignore'):
            R[ALLM.index(m)] = np.nanmean(w.reshape(360, 2, 720, 2), axis=(1, 3))[land]   # 2 x 2 block mean
    g.close()
    return R


def read_jpl(land, path):
    import datetime as dt
    j = nc.Dataset(path)
    seen = _months_from_days(j['time'][:], dt.datetime(2002, 1, 1))
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    for m, i in seen.items():
        w = np.ma.filled(j['lwe_thickness'][i], np.nan) * 10.0                       # cm -> mm
        R[ALLM.index(m)] = np.roll(w[::-1, :], 360, axis=1)[land]                   # 0..360 -> -180..180, N -> S
    j.close()
    return R


def read_glws(land, path):
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
    for f in sorted(folder.glob(pattern)):
        d = f.name.split('.A')[1][:8]
        by.setdefault('%s-%s' % (d[:4], d[4:6]), []).append(str(f))
    R = np.full((len(ALLM), int(land.sum())), np.nan)
    with Pool(C['workers']) as p:
        for m, o in p.imap(_gldas_month, [(m, by.get(m, []), land) for m in ALLM]):
            if o is not None:
                R[ALLM.index(m)] = o
    return R


def product_path(name):
    p = PRODUCTS[name]
    path = _configured_path(name)
    if 'pattern' in p:
        files = sorted(path.glob(p['pattern'])) if path.is_dir() else []
        if not files:
            raise FileNotFoundError('%s: no %s in %s; see download.py' % (name, p['pattern'], path))
        return path if p['kind'] == 'gldas' else files[-1]
    if not path.exists():
        raise FileNotFoundError('%s: %s missing; run download.py' % (name, path))
    return path


def prepare(name, land, months):
    """monthly anomaly grid on the land domain [mm] and SH coefficients of the EWH [m], cached"""
    out = cache_dir() / ('shc_%s.npy' % name)
    p, path = PRODUCTS[name], product_path(name)
    files = sorted(path.glob(p['pattern'])) if path.is_dir() else [path]
    source = '%s\n%d files, %d bytes\n' % (path, len(files), sum(f.stat().st_size for f in files))
    tag = out.with_suffix('.source.txt')
    if out.exists() and tag.exists() and tag.read_text() == source:
        return                                                          # cached for the same input file(s)
    t0 = time.time()
    R = {'csr': read_csr, 'jpl': read_jpl, 'glws': read_glws}[p['kind']](land, path) if p['kind'] != 'gldas' \
        else read_gldas(land, path, p['pattern'])
    R = R[[ALLM.index(m) for m in months]]
    okm = np.isfinite(R).mean(1) > 0.5                                  # months with data
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)                  # cells without any data stay NaN -> 0
        A = R - np.nanmean(R[okm], axis=0)[None, :]                          # anomaly about the product's own mean
    nanfrac = float(np.mean(~np.isfinite(A[okm])))
    A[~np.isfinite(A)] = 0.0
    A[~okm] = np.nan
    np.savez_compressed(cache_dir() / ('field_%s.npz' % name), months=np.array(months), anom_mm=A.astype(np.float32),
                        month_ok=okm)
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
    tag.write_text(source)
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
    coef = {n: np.load(cache_dir() / ('shc_%s.npy' % n)) for n in names}
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


def boot_median(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0
    b = np.median(x[rng.integers(0, len(x), (C['n_boot'], len(x)))], axis=1)
    return float(np.median(x)), float(np.percentile(b, 5)), float(np.percentile(b, 95)), len(x)


def scores(names, S, GN, U, months, suffix=''):
    rng = np.random.default_rng(C['seed'])
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
    L = pd.DataFrame(rows)
    W = L.pivot(index='station', columns='product', values=SCORES + ['amp_model_mm', 'trend_model_mm_yr'])
    W.columns = ['%s|%s' % (p, s) for s, p in W.columns]
    geo = ['dist_unit_boundary_km', 'in_footprint', 'wmo_region', 'boundary_class']
    full = 'JPL' if 'JPL' in names else names[0]           # JPL defines the common months: GNSS figures on all of them
    base = S[['lat', 'lon'] + geo].join(L[L['product'] == full].set_index('station')[['n_months', 'amp_gnss_mm',
                                                                                     'wrms_gnss_mm']])
    PS = base.join(W[sorted(W.columns, key=lambda c: (names.index(c.split('|')[0]), c))])
    out = output_dir()
    PS.to_csv(out / ('gnss_scores_per_station%s.csv' % suffix), float_format='%.4g')
    PS['japan'] = PS.lat.between(24, 46) & PS.lon.between(122, 146)     # GEONET: extra groups without Japan
    groups = [('global', PS)] + [('region:' + k, v) for k, v in PS.groupby('wmo_region')]
    groups += [('distance:' + k, v) for k, v in PS.groupby('boundary_class')]
    groups += [('distance_exJapan:' + k, v) for k, v in PS[~PS.japan].groupby('boundary_class')]
    groups += [('global_exJapan', PS[~PS.japan])]
    groups += [('footprint:' + ('in' if k else 'out(>40km)'), v) for k, v in PS.groupby('in_footprint')]
    summ = []
    for grp, sub in groups:
        for p in names:
            for s in SCORES:
                med, lo, hi, n = boot_median(sub['%s|%s' % (p, s)], rng)
                d = dict(group=grp, product=p, score=s, n=n, median=med, p05=lo, p95=hi)
                if s == 'amp_ratio':
                    d['median_abs_log'] = float(np.nanmedian(np.abs(np.log(sub['%s|%s' % (p, s)]))))
                if s == 'phase_diff_d':
                    d['median_abs'] = float(np.nanmedian(np.abs(sub['%s|%s' % (p, s)])))
                summ.append(d)
    SM = pd.DataFrame(summ)
    SM.to_csv(out / ('gnss_summary%s.csv' % suffix), index=False, float_format='%.4g')
    grp_cols = [('global', 'global'), ('global_exJapan', 'global ex Japan'), ('region:Europe', 'Europe'),
                ('region:North America, Central America and the Caribbean', 'N America'), ('region:Asia', 'Asia'),
                ('region:South America', 'S America')]
    rows = []
    for p in names:
        d = dict(product=p)
        for grp, gl in grp_cols:
            x = SM[(SM.group == grp) & (SM['product'] == p)].set_index('score')
            if x.empty:
                continue
            d['%s n' % gl] = int(x.loc['wrms_red', 'n'])
            d['%s WRMS red %%' % gl] = round(100 * x.loc['wrms_red', 'median'], 1)
            d['%s amp err factor' % gl] = round(float(np.exp(x.loc['amp_ratio', 'median_abs_log'])), 2)
            d['%s r_total' % gl] = round(x.loc['r_total', 'median'], 2)
        rows.append(d)
    T = pd.DataFrame(rows)
    T.to_csv(out / ('table_gnss%s.csv' % suffix), index=False)
    return T


if __name__ == '__main__':
    t00 = time.time()
    names = C['products']
    land, months = domain()
    S = stations()
    S = S[S.used.astype(bool)]
    log('%d stations, %d common months, %d land cells, products %s' % (len(S), len(months), land.sum(), names))
    ref = gnss_reference()
    if C['gnss_source'] == 'reference':
        GN, suffix = ref, ''
    else:
        variant = 'rebuild' if '--default-cleaning' not in sys.argv else 'default'
        GN = gnss_rebuild(variant, list(S.index))
        log('rebuilt vs reference: ' + compare_gnss(GN, ref))
        suffix = '_rebuild' if variant == 'rebuild' else '_default_cleaning'
        np.savez(output_dir() / ('gnss%s.npz' % suffix), stations=np.array(list(GN)),
                 data=np.array([GN[k] for k in GN], dtype=object))
        S = S.loc[[s for s in S.index if s in GN]]
    for n in names:
        prepare(n, land, months)
    t0 = time.time()
    U = forward(names, S.lat.values, S.lon.values)
    log('forward model: %d stations x %d products, %.0f s' % (len(S), len(names), time.time() - t0))
    T = scores(names, S, GN, U, months, suffix)
    log('scores written to %s' % output_dir())
    print(T[['product'] + [c for c in T.columns if c.startswith('global ')]].to_string(index=False))
    log('total %.0f s' % (time.time() - t00))

"""
End to end: synthetic daily GNSS files from a known seasonal load -> preprocessing -> monthly inversion -> comparison
with a synthetic GRACE mascon file holding the same load.
"""
import numpy as np
import pytest

import src as pgw
from src.preprocessing import GNSSNetwork
from src.inversion.workflow import MonthlyInversion
from src.observations.grace import Mascon, compare

RES = 0.5
LAT1 = np.arange(43 + RES / 2, 49, RES)
LON1 = np.arange(3 + RES / 2, 10, RES)
MON = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def pattern(lat, lon):
    """spatial pattern of the load [m EWH per unit of the time function]"""
    return 0.15 * np.exp(-((lat - 46) ** 2 + (lon - 6.5) ** 2) / 2.0)


def time_function(t):
    return np.sin(2 * np.pi * t) + 0.3 * np.sin(2 * np.pi * t / 4)       # annual + a 4-year wave


def write_network(folder, cells, lln, n=40, seed=0):
    rng = np.random.default_rng(seed)
    slat, slon = rng.uniform(44, 48, n), rng.uniform(4, 9, n)
    gain = cells.design_matrix(slat, slon, lln, pgw.Displacement.Vertical) @ pattern(cells.lat, cells.lon)
    days = np.arange(np.datetime64('2010-01-01'), np.datetime64('2018-01-01'))
    t = np.array([int(str(d)[:4]) + ((d - np.datetime64(str(d)[:4] + '-01-01')).astype(int) + 0.5) /
                  (366 if int(str(d)[:4]) % 4 == 0 else 365) for d in days])
    head = 'site YYMMMDD yyyy.yyyy __MJD week d reflon _e0(m) __east(m) ____n0(m) _north(m) u0(m) ____up(m) ...\n'
    files = []
    for k in range(n):
        up = gain[k] * time_function(t) + rng.normal(0, 0.003) * (t - 2014) + rng.normal(0, 2e-3, t.size)
        if k == 3:
            up += np.where(t > 2013.7, 0.015, 0.0)                         # an unreported offset
        lines = [head] + ['S%03d %02d%s%02d %.4f 0 0 0 0 0 0.0 0 0.0 0 %.5f 0 0.001 0.001 0.003 0 0 0 %.4f %.4f 0\n'
                          % (k, int(str(d)[2:4]), MON[int(str(d)[5:7]) - 1], int(str(d)[8:10]), tt, u,
                             slat[k], slon[k]) for d, tt, u in zip(days, t, up)]
        f = folder / ('S%03d.tenv3' % k)
        f.write_text(''.join(lines))
        files.append(f)
    return files


def write_grace(path, months):
    import netCDF4 as nc
    lat = np.arange(-89.875, 90, 0.25)
    lon = np.arange(0.125, 360, 0.25)
    lo2, la2 = np.meshgrid((lon + 180) % 360 - 180, lat)
    with nc.Dataset(path, 'w') as f:
        f.createDimension('time', len(months)), f.createDimension('lat', lat.size), f.createDimension('lon', lon.size)
        t = f.createVariable('time', 'f8', ('time',)); t.units = 'days since 2002-01-01 00:00:00'
        f.createVariable('lat', 'f4', ('lat',))[:] = lat
        f.createVariable('lon', 'f4', ('lon',))[:] = lon
        v = f.createVariable('lwe_thickness', 'f4', ('time', 'lat', 'lon')); v.units = 'cm'
        P = pattern(la2, lo2)
        for i, m in enumerate(months):
            tm = int(m[:4]) + (int(m[5:7]) - 0.5) / 12
            t[i] = (np.datetime64(m + '-15') - np.datetime64('2002-01-01')).astype(int)
            v[i] = 100 * P * time_function(tm)


@pytest.fixture(scope='module')
def lln():
    return pgw.LoveNumber().config(lmax=2000, method=pgw.LLN_Data.PREM).get_Love_number()


def test_gnss_to_grace(tmp_path, lln):
    lo2, la2 = np.meshgrid(LON1, LAT1)
    cells = pgw.GreenFunctionBasis.from_grid(LAT1, LON1, RES, kernel='point_mass')
    files = write_network(tmp_path, cells, lln)

    net = GNSSNetwork.from_files(files).clean(detect_jumps=True).monthly().spatial_filter(radius_km=300)
    assert net.stations[3].qc['n_detected_jumps'] == 1 and len(net.stations) == 40
    ds = net.to_dataset(components=('U',), reference=('2011-01', '2016-12'))

    R = pgw.regularization.curvature(cells.lat, cells.lon, RES)
    inv = MonthlyInversion(cells, lln, ds, components=('U',), R=R)
    alpha = inv.choose_alpha(method='gcv')
    res = inv.run(alpha)
    assert len(res.months) == 96 and (res.n_stations == 40).all()

    # regional mean over the station area against the truth
    glat, glon = np.arange(44.125, 48, 0.25), np.arange(4.125, 9, 0.25)
    gl2, gb2 = np.meshgrid(glon, glat)
    w = np.cos(np.deg2rad(gb2))
    tm = np.array([int(m[:4]) + (int(m[5:7]) - 0.5) / 12 for m in res.months])
    truth = np.sum(w * pattern(gb2, gl2)) / w.sum() * time_function(tm)
    truth = truth - truth[[i for i, m in enumerate(res.months) if '2011-01' <= m <= '2016-12']].mean()
    est = res.regional_series(glat, glon)
    assert np.corrcoef(est, truth)[0, 1] > 0.97

    # the same comparison through a GRACE mascon file
    write_grace(tmp_path / 'mascon.nc', res.months)
    g = Mascon.read(tmp_path / 'mascon.nc').anomalies(ds.reference)
    stats = compare(res.months, est, g.months, g.regional_series(glat, glon))
    assert stats['n'] == 96 and stats['r'] > 0.97 and 0.6 < stats['amplitude_ratio'] < 1.4
    mass = res.regional_mass(glat, glon)
    assert mass.shape == (96,) and np.isfinite(mass).all()


def test_far_field_ignores_the_load_inside_and_predicts_the_load_outside(lln):
    from src.observations.grace import far_field
    lat = np.arange(30.125, 60, 0.25)
    lon = np.arange(-20 + 0.125, 30, 0.25)
    lo2, la2 = np.meshgrid(lon, lat)
    inside = (la2 > 40) & (la2 < 50) & (lo2 > 0) & (lo2 < 10)
    ewh = np.stack([np.where(inside, 0.2, 0.0), np.where(inside, 0.0, 0.05)])
    m = Mascon(['2020-01', '2020-02'], lat, lon, ewh)
    months, p = far_field(m, np.array([45.0]), np.array([5.0]), lln, exclude=(40, 50, 0, 10))
    assert abs(p[0, 0]) < 1e-9                                     # load only inside: nothing predicted
    cells = pgw.GreenFunctionBasis.from_grid(lat.reshape(-1, 2).mean(1), lon.reshape(-1, 2).mean(1), 0.5, kernel='disk')
    direct = cells.design_matrix(np.array([45.0]), np.array([5.0]), lln, pgw.Displacement.Vertical) @ \
        ewh[1].reshape(lat.size // 2, 2, lon.size // 2, 2).mean((1, 3)).ravel()
    assert p[1, 0] < 0 and abs(p[1, 0] - direct[0]) < 0.05 * abs(direct[0])   # outside load: as the full forward


def test_remove_far_field_keeps_every_gnss_month(lln):
    """a long mascon gap (as GRACE/GRACE-FO 2017-2018) is bridged by a seasonal fit, so no GNSS month is lost"""
    from src.observations.grace import remove_far_field
    from src.preprocessing import GNSSDataset
    lat = np.arange(30.125, 60, 0.25)
    lon = np.arange(-20 + 0.125, 30, 0.25)
    gm = ['%d-%02d' % (y, m) for y in range(2014, 2021) for m in range(1, 13)]
    keep = [m for m in gm if not ('2017-07' <= m <= '2018-05')]
    tg = np.array([int(m[:4]) + (int(m[5:7]) - 0.5) / 12 for m in keep])
    ewh = 0.1 * np.sin(2 * np.pi * tg)[:, None, None] * np.ones((1, lat.size, lon.size))
    m = Mascon(keep, lat, lon, ewh)
    ds = GNSSDataset(['A'], np.array([45.0]), np.array([5.0]), gm, {'U': np.zeros((len(gm), 1))},
                     {'U': np.full((len(gm), 1), 1e-3)}, ('U',), ('2014-01', '2020-12'))
    out = remove_far_field(ds, m, lln, exclude=(40, 50, 0, 10))
    assert np.isfinite(out.values['U']).all()
    t = np.array([int(x[:4]) + (int(x[5:7]) - 0.5) / 12 for x in gm])
    assert np.corrcoef(out.values['U'][:, 0], -np.sin(2 * np.pi * t))[0, 1] < -0.9   # loading: subsidence when wet


def test_geocentre_is_a_rigid_translation_and_cm_to_cf_removes_it(lln):
    """CM minus CF displacement of a global load = projection of one translation T, at any point and component"""
    import copy
    from src.observations.grace import geocentre, unit_vectors, cm_to_cf
    from src.preprocessing import GNSSDataset
    lat = np.arange(-89, 90, 2.0)
    lon = np.arange(-179, 180, 2.0)
    months = ['2015-%02d' % m for m in range(1, 13)]
    t = 2015 + (np.arange(12) + 0.5) / 12
    shape = np.cos(np.radians(lat - 40))[:, None] * (np.abs(lon + 100) < 40)[None, :]  # a continental-scale load
    ewh = 0.2 * np.sin(2 * np.pi * t)[:, None, None] * shape[None]
    g = Mascon(months, lat, lon, ewh)
    m, T = geocentre(g, lln, resolution=2.0, n_points=200)
    assert np.abs(T).max() > 1e-4                                   # a load off the pole moves the centre of mass
    plat, plon = np.array([45.0, -30.0, 10.0]), np.array([-120.0, 140.0, 20.0])
    cells = pgw.GreenFunctionBasis.from_grid(lat, lon, 2.0, kernel='disk')
    u = unit_vectors(plat, plon)
    for c, obs in (('U', pgw.Displacement.Vertical), ('N', pgw.Displacement.North), ('E', pgw.Displacement.East)):
        d = 0
        for frame, sign in ((pgw.Frame.CM, 1), (pgw.Frame.CF, -1)):
            l_ = copy.deepcopy(lln)
            l_.convert(frame)
            d = d + sign * (ewh.reshape(12, -1) @ cells.design_matrix(plat, plon, l_, obs).T)
        assert np.allclose(T @ u[c].T, d, atol=1e-6)
    # a GNSS series in CM = the CF series + the projection: cm_to_cf brings it back
    cf = 1e-3 * np.cos(2 * np.pi * t)[:, None] * np.ones((1, 3))
    ds = GNSSDataset(['A', 'B', 'C'], plat, plon, months, {'U': cf + T @ u['U'].T}, {'U': np.full((12, 3), 1e-3)},
                     ('U',), ('2015-01', '2015-12'))
    out = cm_to_cf(ds, g, lln, detrend=False, resolution=2.0, seasonal=False)
    assert np.allclose(out.values['U'], cf, atol=1e-7)


def test_common_mode_keeps_a_network_wide_horizontal_signal_out_of_the_load(lln):
    """a uniform north offset (not loading) is estimated as common mode, and the load is recovered"""
    from src.preprocessing import GNSSDataset
    cells = pgw.GreenFunctionBasis.from_grid(LAT1, LON1, RES, kernel='disk')
    rng = np.random.default_rng(3)
    slat, slon = rng.uniform(44.2, 47.8, 60), rng.uniform(4.2, 8.8, 60)
    months = ['2015-%02d' % m for m in range(1, 13)]
    x = np.array([pattern(cells.lat, cells.lon) * np.sin(2 * np.pi * (k + 0.5) / 12) for k in range(12)])
    obs = {c: x @ cells.design_matrix(slat, slon, lln, d).T for c, d in
           (('U', pgw.Displacement.Vertical), ('N', pgw.Displacement.North), ('E', pgw.Displacement.East))}
    offset = 3e-3 * np.cos(2 * np.pi * (np.arange(12) + 0.5) / 12)              # 3 mm, not loading
    obs['N'] = obs['N'] + offset[:, None]
    ds = GNSSDataset([str(i) for i in range(60)], slat, slon, months, obs,
                     {c: np.full((12, 60), 5e-4) for c in obs}, ('U', 'N', 'E'), (months[0], months[-1]))
    R = pgw.regularization.curvature(cells.lat, cells.lon, RES)
    err = {}
    for cm in ((), ('N', 'E')):
        inv = MonthlyInversion(cells, lln, ds, R=R, common_mode=cm)
        res = inv.run(0.1)
        err[cm] = np.sqrt(np.mean((res.x - x) ** 2)) / np.sqrt(np.mean(x ** 2))
        if cm:
            assert np.allclose(res.meta['common_mode']['N'], offset, atol=5e-4)
            assert np.abs(res.meta['common_mode']['E']).max() < 5e-4
    assert err[('N', 'E')] < 0.5 * err[()]


def test_compare_detrended_ignores_a_trend_missing_from_the_estimate():
    months = ['%d-%02d' % (y, m) for y in range(2010, 2016) for m in range(1, 13)]
    t = 2010 + (np.arange(72) + 0.5) / 12
    seasonal = 0.1 * np.sin(2 * np.pi * t)
    grace = seasonal - 0.02 * (t - 2013)                 # drought trend
    raw, det = compare(months, seasonal, months, grace), compare(months, seasonal, months, grace, detrend=True)
    assert raw['r'] < 0.97 and det['r'] > 0.999 and abs(det['amplitude_ratio'] - 1) < 0.01


def test_daily_inversion_matches_the_truth_and_averages_to_months(tmp_path, lln):
    """daily epochs through the same chain; compare() averages them to months"""
    cells = pgw.GreenFunctionBasis.from_grid(LAT1, LON1, RES, kernel='point_mass')
    files = write_network(tmp_path, cells, lln, n=30, seed=1)
    net = GNSSNetwork.from_files(files).clean().daily()
    ds = net.to_dataset(components=('U',), reference=('2011-01', '2016-12'))
    assert len(ds.months[0]) == 10 and len(ds.months) == 2922                  # 'YYYY-MM-DD', 2010-2017
    assert 1.9e-3 < np.nanmedian(ds.sigma['U']) <= 3.0e-3                      # daily noise 2 mm, formal 3 mm
    R = pgw.regularization.curvature(cells.lat, cells.lon, RES)
    inv = MonthlyInversion(cells, lln, ds, components=('U',), R=R)
    days = [d for d in inv.months(20) if '2015' <= d < '2016']
    alpha = inv.choose_alpha(months=days[::10])
    res = inv.run(alpha, months=days)
    glat, glon = np.arange(44.125, 48, 0.25), np.arange(4.125, 9, 0.25)
    gl2, gb2 = np.meshgrid(glon, glat)
    w = np.cos(np.deg2rad(gb2))
    from src.preprocessing.series import epoch_to_decyear
    td = np.array([epoch_to_decyear(d) for d in res.months])
    truth = np.sum(w * pattern(gb2, gl2)) / w.sum() * time_function(td)
    est = res.regional_series(glat, glon)
    assert np.corrcoef(est, truth)[0, 1] > 0.8                                 # daily noise 2 mm
    months = sorted({d[:7] for d in res.months})
    tm = np.array([epoch_to_decyear(m) for m in months])
    st = compare(res.months, est - est.mean(), months, np.sum(w * pattern(gb2, gl2)) / w.sum() * time_function(tm))
    assert st['n'] == 12 and st['r'] > 0.98


def test_remove_common_mode_by_stacking():
    from src.preprocessing import GNSSDataset
    rng = np.random.default_rng(5)
    lat, lon = rng.uniform(27, 33, 40), rng.uniform(-99, -91, 40)
    days = ['2017-08-%02d' % d for d in range(1, 31)]
    cme = 3e-3 * rng.standard_normal(30)
    signal = np.where((lat < 30.5) & (lon > -96), -0.02, 0.0)[None, :] * (np.arange(30) >= 25)[:, None]
    ds = GNSSDataset([str(i) for i in range(40)], lat, lon, days, {'U': signal + cme[:, None]},
                     {'U': np.full((30, 40), 3e-3)}, ('U',), (days[0], days[23]))
    out = ds.remove_common_mode((31, 33, -99, -97) if ((lat > 31) & (lon < -97)).sum() >= 3 else (30.5, 33, -99, -91))
    assert np.allclose(out.meta['common_mode']['U'], cme, atol=1e-12)
    assert np.allclose(out.values['U'], signal, atol=1e-12)

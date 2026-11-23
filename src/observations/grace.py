"""
GRACE/GRACE-FO mascon solutions for comparison with GNSS-inverted loads.

    from src.observations.grace import Mascon, compare
    g = Mascon.read(path)                                   # CSR RL06.x mascons (NetCDF), JPL in the same layout
    g = g.anomalies(reference=('2004-01', '2009-12'))       # relative to the same period as the GNSS data
    s = g.regional_series(lat, lon, inside)                 # area-weighted mean EWH over a region [m]
    stats = compare(months_gnss, s_gnss, g.months, s)       # correlation, RMS, NSE, amplitude ratio

EWH is in metres throughout (the CSR files store centimetres in `lwe_thickness`).
"""
from dataclasses import dataclass

import numpy as np


def _units(var, default=None):
    """the units attribute of a NetCDF variable ('units', or 'Units' as in the CSR mascon files)"""
    for k in ('units', 'Units', 'unit'):
        if k in var.ncattrs():
            return getattr(var, k)
    return default


@dataclass
class Mascon:
    months: list                    # 'YYYY-MM'
    lat: np.ndarray                 # (nlat,) increasing
    lon: np.ndarray                 # (nlon,) -180..180, increasing
    ewh: np.ndarray                 # (n_month, nlat, nlon) [m]

    @classmethod
    def read(cls, path, variable='lwe_thickness', units_to_m=None, box=None):
        """
        Read a mascon NetCDF (time in days since an epoch, lat, lon, lwe_thickness). Several solutions in one
        calendar month are averaged.
        :param units_to_m: scale to metres; default from the variable's units attribute (cm -> 0.01)
        :param box: (lat_min, lat_max, lon_min, lon_max) [deg]; read only this region (a global 0.25-degree file
                    holds about 1 GB per copy in memory). Default: the whole grid.
        """
        import netCDF4 as nc
        with nc.Dataset(path) as f:
            t = f['time']
            dates = nc.num2date(t[:], _units(t), only_use_cftime_datetimes=False, only_use_python_datetimes=True)
            lat = np.asarray(f['lat'][:], dtype=float)
            lon = (np.asarray(f['lon'][:], dtype=float) + 180) % 360 - 180
            ilat, ilon = np.arange(lat.size), np.arange(lon.size)
            if box is not None:
                ilat = np.flatnonzero((lat >= box[0]) & (lat <= box[1]))
                w = (lon >= box[2]) & (lon <= box[3]) if box[2] <= box[3] else (lon >= box[2]) | (lon <= box[3])
                ilon = np.flatnonzero(w)
                if ilat.size == 0 or ilon.size == 0:
                    raise ValueError('no mascon grid node in the box %s' % (box,))
            v = f[variable]
            if units_to_m is None:
                u = _units(v, 'cm').lower()
                units_to_m = 0.01 if u.startswith('c') else (0.001 if u.startswith('mm') else 1.0)
            # contiguous ranges are read as slices; otherwise orthogonal index arrays
            def sel(i):
                return slice(i[0], i[-1] + 1) if np.all(np.diff(i) == 1) else i
            data = np.asarray(v[:, sel(ilat), sel(ilon)], dtype=np.float32) * np.float32(units_to_m)
        lat, lon = lat[ilat], lon[ilon]
        order = np.argsort(lon)
        lon, data = lon[order], data[:, :, order]
        if lat[0] > lat[-1]:
            lat, data = lat[::-1], data[:, ::-1]
        keys = np.array(['%04d-%02d' % (d.year, d.month) for d in dates])
        months = sorted(set(keys))
        ewh = np.stack([data[keys == m].mean(axis=0) for m in months]).astype(float)
        return cls(months, lat, lon, ewh)

    def anomalies(self, reference):
        """relative to the mean over the months of the reference period ('YYYY-MM', 'YYYY-MM'), inclusive"""
        sel = [i for i, m in enumerate(self.months) if reference[0] <= m <= reference[1]]
        if not sel:
            raise ValueError('no GRACE month in the reference period %s' % (reference,))
        return Mascon(self.months, self.lat, self.lon, self.ewh - self.ewh[sel].mean(axis=0))

    def sample(self, lat, lon):
        """EWH at points (nearest grid node), (n_month, N)"""
        i = np.clip(np.searchsorted(self.lat, lat), 1, self.lat.size - 1)
        i -= (np.abs(self.lat[i - 1] - lat) < np.abs(self.lat[i] - lat))
        lon = (np.asarray(lon) + 180) % 360 - 180
        j = np.clip(np.searchsorted(self.lon, lon), 1, self.lon.size - 1)
        j -= (np.abs(self.lon[j - 1] - lon) < np.abs(self.lon[j] - lon))
        return self.ewh[:, i, j]

    def regional_series(self, lat, lon, inside=None):
        """area-weighted mean EWH [m] over a region on the grid (lat, lon), as InversionResult.regional_series"""
        lo2, la2 = np.meshgrid(lon, lat)
        m = np.ones(la2.shape, bool) if inside is None else (inside(la2, lo2) if callable(inside) else inside)
        w = np.cos(np.deg2rad(la2[m]))
        return self.sample(la2[m], lo2[m]) @ w / w.sum()


def far_field(mascon, lat, lon, lln, exclude, component='U', resolution=0.5, months=None):
    """
    Displacement at stations caused by the load OUTSIDE a region, predicted from a mascon solution: to be subtracted
    from GNSS series before a regional inversion (as Shen et al. 2020; Argus et al. 2022). Without it, water outside
    the inverted cells (no stations there to separate it) is absorbed by the cells and inflates the result.
    :param mascon: Mascon (anomalies), read over a box well beyond the region (e.g. 10-15 degrees)
    :param exclude: (lat_min, lat_max, lon_min, lon_max) covered by the inversion cells: load there is NOT predicted
    :param component: 'U', 'N' or 'E'
    :param resolution: cells used for the prediction [deg] (multiple of the mascon grid; 0.5 for 0.25-degree CSR)
    :return: months, (n_month, N) displacement [m]
    """
    from src.basis.green.gf_basis import GreenFunctionBasis
    from src.utils.enums import Displacement
    k = int(round(resolution / abs(mascon.lat[1] - mascon.lat[0])))
    nla, nlo = (mascon.lat.size // k) * k, (mascon.lon.size // k) * k
    la = mascon.lat[:nla].reshape(-1, k).mean(1)
    lo = mascon.lon[:nlo].reshape(-1, k).mean(1)
    E = mascon.ewh[:, :nla, :nlo].reshape(len(mascon.months), la.size, k, lo.size, k).mean((2, 4))
    lo2, la2 = np.meshgrid(lo, la)
    out = ~((la2 > exclude[0]) & (la2 < exclude[1]) & (lo2 > exclude[2]) & (lo2 < exclude[3]))
    cells = GreenFunctionBasis.from_grid(la, lo, resolution, mask=out, kernel='disk')
    obs = {'U': Displacement.Vertical, 'N': Displacement.North, 'E': Displacement.East}[component]
    A = cells.design_matrix(np.asarray(lat), np.asarray(lon), lln, obs)
    pred = E[:, out] @ A.T
    if months is None:
        return list(mascon.months), pred
    idx = {m: i for i, m in enumerate(mascon.months)}
    return list(months), np.array([pred[idx[m]] if m in idx else np.full(pred.shape[1], np.nan) for m in months])


def _harmonics(x):
    return np.c_[np.ones_like(x), x - 2014, np.sin(2 * np.pi * x), np.cos(2 * np.pi * x),
                 np.sin(4 * np.pi * x), np.cos(4 * np.pi * x)]


def _decyear(epochs):
    """'YYYY-MM' or 'YYYY-MM-DD' -> decimal years"""
    from src.preprocessing.series import epoch_to_decyear
    return np.array([epoch_to_decyear(e) for e in epochs])


def to_monthly(epochs, values, min_days=10):
    """daily series ('YYYY-MM-DD') -> calendar-month means (months with fewer than min_days values dropped);
    monthly input is returned unchanged"""
    epochs = list(epochs)
    values = np.asarray(values, dtype=float)
    if not epochs or len(epochs[0]) == 7:
        return epochs, values
    keys = np.array([e[:7] for e in epochs])
    months, out = [], []
    for k in sorted(set(keys)):
        sel = (keys == k) & np.isfinite(values)
        if sel.sum() >= min_days:
            months.append(k)
            out.append(values[sel].mean())
    return months, np.array(out)


def _to_months(t, tg, y, max_gap=3.5 / 12):
    """series y at mascon epochs tg -> GNSS epochs t: interpolated within max_gap of a mascon month, elsewhere (longer
    gaps, outside the mascon period) from a fit of mean, trend, annual and semi-annual terms"""
    p = np.interp(t, tg, y, left=np.nan, right=np.nan)
    p[np.array([np.min(np.abs(tg - x)) for x in t]) > max_gap] = np.nan
    if np.isnan(p).any():
        coef = np.linalg.lstsq(_harmonics(tg), y, rcond=None)[0]
        p[np.isnan(p)] = _harmonics(t[np.isnan(p)]) @ coef
    return p


def geocentre(mascon, lln, resolution=2.0, n_points=400):
    """
    Annual-scale motion of the centre of mass (CM) relative to the centre of the surface figure (CF), from a GLOBAL
    mascon solution (degree 1 of the load, Love numbers of both frames).
    :return: months, T (n_month, 3) [m] in Earth-fixed X, Y, Z, such that a displacement in CM minus the same
             displacement in CF is the projection of T on the direction of that component (unit_vectors)
    """
    import copy
    from src.basis.green.gf_basis import GreenFunctionBasis
    from src.utils.enums import Displacement, Frame
    k = int(round(resolution / abs(mascon.lat[1] - mascon.lat[0])))
    nla, nlo = (mascon.lat.size // k) * k, (mascon.lon.size // k) * k
    if np.ptp(mascon.lat[:nla]) < 170 or np.ptp(mascon.lon[:nlo]) < 350:
        raise ValueError('geocentre needs a global mascon solution (read without box)')
    la = mascon.lat[:nla].reshape(-1, k).mean(1)
    lo = mascon.lon[:nlo].reshape(-1, k).mean(1)
    E = mascon.ewh[:, :nla, :nlo].reshape(len(mascon.months), la.size, k, lo.size, k).mean((2, 4))
    E = E.reshape(len(mascon.months), -1)
    cells = GreenFunctionBasis.from_grid(la, lo, resolution, kernel='disk')
    i = np.arange(n_points) + 0.5                                  # Fibonacci points: even cover of the sphere
    plat = np.degrees(np.arcsin(1 - 2 * i / n_points))
    plon = (np.degrees(np.pi * (1 + 5 ** 0.5) * i) + 180) % 360 - 180
    d = 0
    for frame, sign in ((Frame.CM, 1), (Frame.CF, -1)):
        l_ = copy.deepcopy(lln)
        l_.convert(frame)
        d = d + sign * (E @ cells.design_matrix(plat, plon, l_, Displacement.Vertical).T)
    G = unit_vectors(plat, plon)['U']
    T = np.linalg.lstsq(G, d.T, rcond=None)[0].T
    return list(mascon.months), T


def unit_vectors(lat, lon):
    """{'U', 'N', 'E'}: (N, 3) local unit vectors in Earth-fixed X, Y, Z"""
    f, l = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    return {'U': np.c_[np.cos(f) * np.cos(l), np.cos(f) * np.sin(l), np.sin(f)],
            'N': np.c_[-np.sin(f) * np.cos(l), -np.sin(f) * np.sin(l), np.cos(f)],
            'E': np.c_[-np.sin(l), np.cos(l), np.zeros_like(f)]}


def cm_to_cf(dataset, mascon, lln, detrend=True, resolution=2.0, seasonal=True):
    """
    GNSSDataset converted from CM (as the NGL series) to CF with the geocentre motion of the SAME global mascon
    solution that the results are compared with (geocentre), instead of an external geocentre model. With NGL's own
    translation (Altamimi et al. 2023) the GNSS seasonal cycle in California lagged the GRACE prediction by 3 weeks
    and was 20% larger; with the GRACE-consistent one the lag is 5 days and the ratio about 1.
    :param detrend: remove the linear trend over the GNSS months, as for the station series
    :param seasonal: use only the mean, annual and semi-annual terms of the GRACE geocentre motion. Its month-to-month
                     part (1-3 mm RMS per axis, largest in Z) is mostly GRACE degree-1 noise: kept, it lowered the
                     agreement of regional GNSS water storage with GRACE in the Pacific Northwest from r = 0.81 to 0.73
    """
    import copy
    out = copy.deepcopy(dataset)
    months, T = geocentre(mascon, lln, resolution)
    t, tg = _decyear(dataset.months), _decyear(months)
    if seasonal:
        cols = [0, 2, 3, 4, 5]                                                  # no trend: removed per station anyway
        T = np.array([_harmonics(t)[:, cols] @ np.linalg.lstsq(_harmonics(tg)[:, cols], T[:, j], rcond=None)[0]
                      for j in range(3)]).T
    else:
        T = np.array([_to_months(t, tg, T[:, j]) for j in range(3)]).T      # (n_month, 3)
    u = unit_vectors(dataset.lat, dataset.lon)
    for c in dataset.components:
        p = T @ u[c].T                                                          # (n_month, N)
        for k in range(p.shape[1]):
            if detrend:
                ok = np.isfinite(dataset.values[c][:, k])
                if ok.sum() > 2:
                    p[:, k] -= np.polyval(np.polyfit(t[ok], p[ok, k], 1), t)
            out.values[c][:, k] = dataset.values[c][:, k] - p[:, k]
    return out


def remove_far_field(dataset, mascon, lln, exclude, resolution=0.5, smooth=1, detrend=True):
    """
    GNSSDataset with the far-field displacement (far_field) subtracted from every component. Months without a mascon
    solution are interpolated from the neighbouring months (gaps up to 3 months); in longer gaps (GRACE/GRACE-FO,
    2017-2018) and outside the mascon period the far field is taken from a fit of mean, trend, annual and semi-annual
    terms, so that no GNSS month is lost.
    :param smooth: running mean over this many mascon months (odd; 1 = none): the far field is smooth in time, and
                   smoothing keeps month-to-month GRACE noise out of the GNSS series
    :param detrend: remove a linear trend from the prediction over the GNSS months, as preprocessing removed each
                    station's trend (clean(remove='trend+steps')); otherwise GRACE's trend outside the region would be
                    injected into detrended GNSS series. Use False only if the GNSS trends were kept.
    """
    import copy
    out = copy.deepcopy(dataset)
    t, tg = _decyear(dataset.months), _decyear(mascon.months)
    for c in dataset.components:
        _, pred = far_field(mascon, dataset.lat, dataset.lon, lln, exclude, c, resolution)
        if smooth > 1:
            h = smooth // 2
            pred = np.array([pred[max(0, i - h):i + h + 1].mean(axis=0) for i in range(len(pred))])
        for k in range(pred.shape[1]):
            p = _to_months(t, tg, pred[:, k])
            if detrend:
                ok = np.isfinite(p) & np.isfinite(dataset.values[c][:, k])
                if ok.sum() > 2:
                    p = p - np.polyval(np.polyfit(t[ok], p[ok], 1), t)
            out.values[c][:, k] = dataset.values[c][:, k] - p
    return out


def compare(months_a, a, months_b, b, detrend=False):
    """
    Agreement of two monthly series on their common months (a = estimate, b = reference, e.g. GNSS vs GRACE).
    Daily series ('YYYY-MM-DD') are first averaged to calendar months (to_monthly).
    :param detrend: remove a linear trend from both over the common months. Use it when the station trends were
                    removed in preprocessing (tectonics): the GNSS series then has no trend by construction, and
                    GRACE's (e.g. -10 mm/yr in California, 2006-2022) would otherwise count as disagreement
    :return: dict n, r (correlation), rms (of a - b), nse (Nash-Sutcliffe efficiency of a against b),
             amplitude_ratio (std a / std b), bias (mean a - b), months
    """
    months_a, a = to_monthly(months_a, a)
    months_b, b = to_monthly(months_b, b)
    ib = {m: i for i, m in enumerate(months_b)}
    common = [m for m in months_a if m in ib]
    ia = {m: i for i, m in enumerate(months_a)}
    x = np.array([a[ia[m]] for m in common], dtype=float)
    y = np.array([b[ib[m]] for m in common], dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if detrend and x.size >= 3:
        tc = _decyear([m for m, k in zip(common, ok) if k])
        x = x - np.polyval(np.polyfit(tc, x, 1), tc)
        y = y - np.polyval(np.polyfit(tc, y, 1), tc)
    if x.size < 3:
        return dict(n=int(x.size), r=np.nan, rms=np.nan, nse=np.nan, amplitude_ratio=np.nan, bias=np.nan, months=[])
    return dict(n=int(x.size), r=float(np.corrcoef(x, y)[0, 1]), rms=float(np.sqrt(np.mean((x - y) ** 2))),
                nse=float(1 - np.sum((x - y) ** 2) / np.sum((y - y.mean()) ** 2)),
                amplitude_ratio=float(np.std(x) / np.std(y)), bias=float(np.mean(x - y)),
                months=[m for m, k in zip(common, ok) if k])

"""
Monthly inversion of a GNSS network for the surface load, with any basis:

    from pygwater.inversion.workflow import MonthlyInversion
    inv = MonthlyInversion(basis, lln, dataset, components=('U',), R=regularization.curvature(lat, lon, res))
    alpha = inv.choose_alpha(method='gcv')                    # one alpha for all months (median of monthly choices)
    res = inv.run(alpha, lb=0.0)                              # InversionResult: x per month
    s = res.regional_series(lat, lon, inside)                 # area-weighted mean EWH over a region [m]

The design matrix of every station is computed once; each month uses the rows of the stations with data in that
month. Observation weights are, by default, one sigma per station and component (the median of its monthly sigmas),
so months sharing the same stations share one decomposition of the normal matrix and the whole series is fast.
"""
from contextlib import nullcontext
from dataclasses import dataclass, field

import numpy as np

try:
    from threadpoolctl import threadpool_limits
except ImportError:                                  # optional: without it the BLAS default is used
    threadpool_limits = None

from pygwater.utils.enums import Displacement
from pygwater.utils.constants import EarthConstant
from pygwater.inversion.solver import RegularizedLeastSquares

_OBS = {'U': Displacement.Vertical, 'N': Displacement.North, 'E': Displacement.East}


class MonthlyInversion:

    def __init__(self, basis, lln, dataset, components=None, R=None, weights='station', sigma_floor=1e-4,
                 cache_gb=1.0, threads=1, common_mode=()):
        """
        :param basis: any pygwater Basis (GreenFunctionBasis, SphericalHarmonicBasis, SlepianBasis)
        :param lln: LoveNumber object (choose the frame consistent with the GNSS series, e.g. CF for seasonal loads)
        :param dataset: GNSSDataset (pygwater.preprocessing)
        :param components: subset of dataset.components, e.g. ('U',) or ('U', 'N', 'E')
        :param R: regularisation matrix (pygwater.inversion.regularization); default identity
        :param weights: 'station' (one sigma per station and component, the median over months; fast) or 'month'
                        (the sigma of each monthly value)
        :param sigma_floor: lower bound of the sigmas [m]
        :param cache_gb: memory for cached decompositions (one per distinct set of stations); each takes about
                         3 n^2 x 8 bytes for n parameters, so large grids keep fewer
        :param common_mode: components (e.g. ('N', 'E')) with one unregularised offset per month shared by all
                            stations, estimated with the load. GNSS horizontals carry a network-wide seasonal signal
                            that is not loading (2.3 mm annual in the north component over the Pacific Northwest,
                            against 0.6 mm predicted from GRACE); without this term it is mapped into water
        :param threads: BLAS threads during the inversion (needs threadpoolctl; None keeps the BLAS default). The
                        matrices are small, so extra threads mostly compete: with the default of one thread per
                        logical core a 20-thread desktop took 413 s instead of 10 s, and 4 threads on a 2-core
                        machine stalled. Raise it only for very large grids
        """
        self.basis, self.lln, self.ds = basis, lln, dataset
        self.components = tuple(components or dataset.components)
        self.R = R
        self.weights = weights
        self.sigma_floor = sigma_floor
        self.A = {c: basis.design_matrix(dataset.lat, dataset.lon, lln, _OBS[c]) for c in self.components}
        self.common_mode = tuple(c for c in self.components if c in tuple(common_mode))
        self.n_load = basis.n_param
        if self.common_mode:
            n_cm = len(self.common_mode)
            for c in self.components:
                cols = np.zeros((self.A[c].shape[0], n_cm))
                if c in self.common_mode:
                    cols[:, self.common_mode.index(c)] = 1.0
                self.A[c] = np.hstack([self.A[c], cols])
            n = basis.n_param
            R0 = np.eye(n) if R is None else np.asarray(R, dtype=float)
            self.R = np.zeros((n + n_cm, n + n_cm))
            self.R[:n, :n] = R0
        self.station_sigma = {c: np.maximum(np.nanmedian(dataset.sigma[c], axis=0), sigma_floor)
                              for c in self.components}
        self._solvers = {}
        self._check_resolution()
        n = basis.n_param
        self._max_solvers = max(1, int(cache_gb * 1e9 / (3 * 8.0 * n * n)))
        self.threads = threads

    def _check_resolution(self):
        """warn when a spectral basis (spherical harmonics, Slepian) resolves finer than the network: L > 180/s for a
        median station spacing s [deg] (Yang et al., 2025); the data fit does not reveal this misconfiguration"""
        lmax = getattr(self.basis, 'lmax', None)
        lat, lon = np.radians(self.ds.lat), np.radians(self.ds.lon)
        if lmax is None or lat.size < 2:
            return
        c = np.clip(np.sin(lat)[:, None] * np.sin(lat)[None, :] + np.cos(lat)[:, None] * np.cos(lat)[None, :]
                    * np.cos(lon[:, None] - lon[None, :]), -1, 1)
        d = np.degrees(np.arccos(c))
        np.fill_diagonal(d, np.inf)
        s = float(np.median(d.min(axis=1)))
        self.station_spacing = s
        if s > 0 and lmax > 180.0 / s:
            import warnings
            warnings.warn('basis degree L=%d exceeds 180/s = %d for the median station spacing s = %.2f deg: the '
                          'network cannot resolve it, and the result can show artefacts with an equally good data fit '
                          '(Yang et al., 2025)' % (lmax, int(180.0 / s), s))

    def _limit(self):
        if self.threads is None or threadpool_limits is None:
            return nullcontext()
        return threadpool_limits(limits=self.threads, user_api='blas')

    # --------------------------------------------------------------------------------------------------------------
    def _month(self, month):
        """rows for one month: A (n_comp * N, n), d, sigma, station mask"""
        i = self.ds.months.index(month)
        ok = np.all([np.isfinite(self.ds.values[c][i]) for c in self.components], axis=0)
        A = np.vstack([self.A[c][ok] for c in self.components])
        d = np.concatenate([self.ds.values[c][i, ok] for c in self.components])
        if self.weights == 'station':
            sg = np.concatenate([self.station_sigma[c][ok] for c in self.components])
        else:
            sg = np.concatenate([np.maximum(self.ds.sigma[c][i, ok], self.sigma_floor) for c in self.components])
        return A, d, sg, ok

    def solver(self, month):
        A, d, sg, ok = self._month(month)
        key = ok.tobytes() if self.weights == 'station' else month
        if key not in self._solvers:
            if len(self._solvers) >= self._max_solvers:
                self._solvers.pop(next(iter(self._solvers)))       # drop the oldest
            self._solvers[key] = RegularizedLeastSquares(A, sg, self.R)
        return self._solvers[key], d, ok

    def months(self, min_stations=1):
        return [m for m in self.ds.months if self._month(m)[3].sum() >= min_stations]

    def choose_alpha(self, method='gcv', alphas=None, months=None, min_stations=5, **kwargs):
        """
        One regularisation parameter for the whole series: the median of the choices of the individual months
        (more stable than a single month). See RegularizedLeastSquares.choose_alpha for the methods.
        :return: alpha; the monthly choices are kept in self.alpha_by_month
        """
        alphas = np.logspace(-6, 5, 45) if alphas is None else alphas
        months = months or self.months(min_stations)
        self.alpha_by_month = {}
        with self._limit():
            self._choose(months, alphas, method, kwargs)
        if not self.alpha_by_month:
            raise ValueError('no month with at least %d stations' % min_stations)
        return float(np.exp(np.median(np.log(list(self.alpha_by_month.values())))))

    def _choose(self, months, alphas, method, kwargs):
        for m in months:
            sv, d, ok = self.solver(m)
            self.alpha_by_month[m] = sv.choose_alpha(d, alphas=alphas, method=method, **kwargs)[0]

    def run(self, alpha, months=None, lb=None, ub=None, min_stations=5):
        """
        :param alpha: one value for all months, or 'per_month' to use self.alpha_by_month (run choose_alpha first)
        :param lb, ub: bounds on the parameters (e.g. lb=0 for a non-negative load); bounded solves are slower
        :return: InversionResult
        """
        months = months or self.months(min_stations)
        with self._limit():
            return self._run(alpha, months, lb, ub)

    def _run(self, alpha, months, lb, ub):
        X, alph, nst, fit, dof = [], [], [], [], []
        for m in months:
            sv, d, ok = self.solver(m)
            a = self.alpha_by_month[m] if alpha == 'per_month' else alpha
            if lb is None and ub is None:
                sol = sv.solve(d, a)
            else:
                sol = sv.solve_bounded(d, a, lb=-np.inf if lb is None else lb, ub=np.inf if ub is None else ub)
            X.append(sol.x), alph.append(a), nst.append(int(ok.sum())), dof.append(sol.dof)
            fit.append(float(np.sqrt(np.mean((sv.A @ sol.x - d) ** 2))))
        X = np.array(X)
        meta = {'common_mode': {c: X[:, self.n_load + i] for i, c in enumerate(self.common_mode)}} \
            if self.common_mode else {}
        return InversionResult(basis=self.basis, months=list(months), x=X[:, :self.n_load], alpha=np.array(alph),
                               n_stations=np.array(nst), fit_rms=np.array(fit), dof=np.array(dof),
                               reference=self.ds.reference, components=self.components, meta=meta)


@dataclass
class InversionResult:
    basis: object
    months: list
    x: np.ndarray                   # (n_month, n_param)
    alpha: np.ndarray
    n_stations: np.ndarray
    fit_rms: np.ndarray             # RMS of the data misfit [m]
    dof: np.ndarray                 # effective number of parameters
    reference: tuple
    components: tuple
    meta: dict = field(default_factory=dict)

    def running_mean(self, n=3):
        """
        A new InversionResult with the coefficients averaged over n consecutive epochs (centred; n odd; the ends use
        the epochs available). For daily epochs this lowers the day-to-day noise by about sqrt(n): the inversion is
        linear, so averaging the solutions equals inverting the averaged data when the stations do not change.
        """
        import copy
        h = n // 2
        out = copy.deepcopy(self)
        out.x = np.array([self.x[max(0, i - h):i + h + 1].mean(axis=0) for i in range(len(self.x))])
        out.meta = dict(self.meta, running_mean=n)
        return out

    def ewh(self, lat, lon):
        """EWH anomaly [m] at points, (n_month, N)"""
        return self.basis.load_synthesis(self.x, np.asarray(lat), np.asarray(lon))

    def regional_series(self, lat, lon, inside=None):
        """
        Area-weighted mean EWH [m] over a region, per month, evaluated on a regular grid.
        :param lat, lon: 1d cell-centred grid covering the region (e.g. 0.25 deg)
        :param inside: optional mask (nlat, nlon) or function (lat2d, lon2d) -> bool; default the whole grid
        """
        lo2, la2 = np.meshgrid(lon, lat)
        m = np.ones(la2.shape, bool) if inside is None else (inside(la2, lo2) if callable(inside) else inside)
        w = np.cos(np.deg2rad(la2[m]))
        return self.ewh(la2[m], lo2[m]) @ w / w.sum()

    def regional_mass(self, lat, lon, inside=None):
        """Mass anomaly [Gt] over the region (regional_series x area x water density)"""
        lo2, la2 = np.meshgrid(lon, lat)
        m = np.ones(la2.shape, bool) if inside is None else (inside(la2, lo2) if callable(inside) else inside)
        dlat, dlon = np.deg2rad(abs(lat[1] - lat[0])), np.deg2rad(abs(lon[1] - lon[0]))
        area = np.sum(np.cos(np.deg2rad(la2[m]))) * dlat * dlon * EarthConstant.radiusm ** 2
        return self.regional_series(lat, lon, m) * area * EarthConstant.rhow / 1e12

    def save(self, path):
        np.savez(path, months=np.array(self.months), x=self.x, alpha=self.alpha, n_stations=self.n_stations,
                 fit_rms=self.fit_rms, dof=self.dof, reference=np.array(self.reference),
                 components=np.array(self.components))

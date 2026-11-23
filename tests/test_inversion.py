"""Tests of src.inversion (solver, regularisation) and src.observations (synthetic observations)."""
import numpy as np

from src.utils.geomath import GeoMathKit
from src.inversion.solver import RegularizedLeastSquares
from src.inversion.regularization import difference_operator, smoothness, identity
from src.observations.synthetic import synthetic_network, checkerboard, add_noise


def test_solver_recovers_overdetermined_problem():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(300, 40)); x = rng.normal(size=40)
    d = A @ x + rng.normal(size=300) * 0.01
    s = RegularizedLeastSquares(A, 0.01)
    alpha, curve = s.choose_alpha(d)
    sol = s.solve(d, alpha)
    assert np.abs(sol.x - x).max() < 0.01 and abs(sol.sigma0 - 1) < 0.2 and sol.dof > 39
    # several epochs in one pass
    X = s.solve(np.vstack([d, 2 * d]), alpha).x
    np.testing.assert_allclose(X[1], 2 * X[0], rtol=1e-10)
    # resolution matrix -> identity, covariance -> (A^T W A)^-1 for small alpha
    np.testing.assert_allclose(s.resolution_matrix(1e-10), np.eye(40), atol=1e-6)
    np.testing.assert_allclose(s.covariance(1e-10), np.linalg.inv(s.N), rtol=1e-4)
    assert s.choose_alpha(d, method='lcurve')[0] > 0 and s.choose_alpha(d, method='wgcv', omega=1.2)[0] > 0


def test_difference_operator_wraps_and_annihilates_constants():
    lat = np.array([1., 1., 3., 3.]); lon = np.array([179., -179., 179., -179.])
    D = difference_operator(lat, lon, 2.0)
    assert D.shape == (4, 4)
    np.testing.assert_allclose(D @ np.ones(4), 0, atol=1e-12)
    R = smoothness(lat, lon, 2.0)
    np.testing.assert_allclose(R, R.T) and np.testing.assert_allclose(R @ np.ones(4), 0, atol=1e-12)
    assert identity(3).shape == (3, 3)


def test_synthetic_network_and_checkerboard():
    lat, lon = synthetic_network(lambda a, o: (a > 20) & (a < 40) & (o > 100) & (o < 120), n_target=40,
                                 min_dist=2.0, bbox=(15, 45, 95, 125), seed=3)
    assert len(lat) == 40 and ((lat > 20) & (lat < 40)).all()
    d = np.array([GeoMathKit.angular_distance(lat[i], lon[i], np.delete(lat, i), np.delete(lon, i)).min()
                  for i in range(len(lat))])
    assert d.min() >= 2.0
    c = checkerboard(np.array([1., 1., 5.]), np.array([1., 5., 1.]), 4.0, 0.1)
    assert c[0] == -c[1] and c[0] == -c[2] and abs(c[0]) == 0.1
    assert np.std(add_noise(np.zeros(20000), 2.0) - 0) > 1.9


def test_choose_alpha_rejects_time_series():
    import pytest
    rng = np.random.default_rng(0)
    A = rng.normal(size=(30, 10))
    solver = RegularizedLeastSquares(A, 1.0)
    with pytest.raises(ValueError):
        solver.choose_alpha(rng.normal(size=(5, 30)))      # (ntime, N_obs): ambiguous for GCV
    alpha, _ = solver.choose_alpha(rng.normal(size=30))
    assert alpha > 0


def test_inversion_warns_when_the_basis_is_finer_than_the_network():
    """a misconfigured spectral basis fits the data as well as a correct one, so PyGWater warns instead"""
    import warnings
    from src.inversion.workflow import MonthlyInversion
    from src.preprocessing import GNSSDataset

    class Basis:                               # stands for a spherical-harmonic or Slepian basis of degree L
        n_param = 4

        def __init__(self, lmax):
            self.lmax = lmax

        def design_matrix(self, lat, lon, lln, obs):
            return np.ones((len(lat), self.n_param))

    lat, lon = np.meshgrid(np.arange(40, 46, 1.0), np.arange(-120, -114, 1.0))      # 1-degree spacing
    n = lat.size
    ds = GNSSDataset(list(map(str, range(n))), lat.ravel(), lon.ravel(), ['2020-01'], {'U': np.zeros((1, n))},
                     {'U': np.full((1, n), 1e-3)}, ('U',), ('2020-01', '2020-01'))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        MonthlyInversion(Basis(120), None, ds)                                  # 180/s = 180: fine
        assert not [x for x in w if 'station spacing' in str(x.message)]
        inv = MonthlyInversion(Basis(480), None, ds)                            # beyond what 1-degree spacing resolves
        assert [x for x in w if 'station spacing' in str(x.message)]
    assert 0.75 < inv.station_spacing < 1.15                                     # coverage spacing of a 1-degree grid


def test_inversion_warns_when_green_function_cells_are_finer_than_the_network():
    import warnings
    from src.inversion.workflow import MonthlyInversion
    from src.preprocessing import GNSSDataset

    class Cells:                               # stands for a Green's-function basis of a given cell size
        n_param = 4

        def __init__(self, resolution):
            self.resolution = resolution

        def design_matrix(self, lat, lon, lln, obs):
            return np.ones((len(lat), self.n_param))

    lat, lon = np.meshgrid(np.arange(40, 46, 1.0), np.arange(-120, -114, 1.0))      # 1-degree spacing
    n = lat.size
    ds = GNSSDataset(list(map(str, range(n))), lat.ravel(), lon.ravel(), ['2020-01'], {'U': np.zeros((1, n))},
                     {'U': np.full((1, n), 1e-3)}, ('U',), ('2020-01', '2020-01'))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        MonthlyInversion(Cells(1.0), None, ds)                                  # cells as large as the spacing
        assert not [x for x in w if 'farther than one cell' in str(x.message)]
        inv = MonthlyInversion(Cells(0.1), None, ds)                            # far finer than the network
        assert [x for x in w if 'farther than one cell' in str(x.message)]
    assert inv.uncovered_fraction > 0.9


def test_running_mean_of_an_inversion_result():
    from src.inversion.workflow import InversionResult
    x = np.arange(10.0)[:, None] * np.ones((1, 3))
    r = InversionResult(basis=None, months=[str(i) for i in range(10)], x=x, alpha=np.ones(10),
                        n_stations=np.ones(10), fit_rms=np.ones(10), dof=np.ones(10), reference=('0', '9'),
                        components=('U',))
    m = r.running_mean(3)
    assert np.allclose(m.x[1:-1, 0], x[1:-1, 0]) and m.x[0, 0] == 0.5 and m.meta['running_mean'] == 3

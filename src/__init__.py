"""
PyGWater: elastic loading deformation of the Earth (GNSS displacement, geoid, gravity) from surface water loads,
and the inversion of GNSS observations for those loads, with a choice of basis functions (Yang et al., 2025).

The code lives in the folder src/ of the repository and is imported as `PyGWater` once installed
(`pip install -e .` in the repository folder, see doc/installation.md). The main classes are available directly:

    import PyGWater as pgw
    lln = pgw.LoveNumber().config(lmax=13000, method=pgw.LLN_Data.Wang).get_Love_number()
    basis = pgw.GreenFunctionBasis.from_grid(lat, lon, 0.25, kernel='cluster_disk')
    A = basis.design_matrix(station_lat, station_lon, lln, pgw.Displacement.Vertical)
    sol = pgw.RegularizedLeastSquares(A, sigma, R).solve(d, alpha)

Subpackages: earth (Love numbers, observables), basis (spherical harmonics, Green's functions, Slepian),
forward, inversion, observations (GNSS, synthetic), mesh (grids, masks), utils; config holds the data root.
"""
__version__ = '0.9.0'

from PyGWater.config import FileTool
from PyGWater.utils.enums import Displacement, LLN_Data, LLN_variable, Frame
from PyGWater.utils.constants import EarthConstant
from PyGWater.earth.love_numbers import LoveNumber
from PyGWater.basis.base import Basis
from PyGWater.basis.spherical_harmonic.sh_basis import SphericalHarmonicBasis
from PyGWater.basis.green.gf_basis import GreenFunctionBasis
from PyGWater.basis.slepian import SlepianBasis
from PyGWater.forward.model import ForwardModel
from PyGWater.inversion.solver import RegularizedLeastSquares
from PyGWater.inversion import regularization
from PyGWater.preprocessing import GNSSNetwork, GNSSDataset, read_tenv3, read_ngl_steps
from PyGWater.inversion.workflow import MonthlyInversion, InversionResult
from PyGWater.observations.grace import Mascon

__all__ = ['FileTool', 'Displacement', 'LLN_Data', 'LLN_variable', 'Frame', 'EarthConstant', 'LoveNumber', 'Basis',
           'SphericalHarmonicBasis', 'GreenFunctionBasis', 'SlepianBasis', 'ForwardModel', 'RegularizedLeastSquares',
           'regularization', 'GNSSNetwork', 'GNSSDataset', 'read_tenv3', 'read_ngl_steps', 'MonthlyInversion', 'InversionResult', 'Mascon',
           '__version__']

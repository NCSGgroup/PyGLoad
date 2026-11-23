"""
PyGWater: elastic loading deformation of the Earth (GNSS displacement, geoid, gravity) from surface water loads,
and the inversion of GNSS observations for those loads, with a choice of basis functions (Yang et al., 2025).

The code lives in the folder src/ of the repository and is imported from the repository folder, as in PyGLDA v2:
scripts put that folder on the path (see examples/ and studies/), then `from src.basis ... import ...` or
`import src as pgw`. The main classes are available directly:

    import src as pgw
    lln = pgw.LoveNumber().config(lmax=13000, method=pgw.LLN_Data.Wang).get_Love_number()
    basis = pgw.GreenFunctionBasis.from_grid(lat, lon, 0.25, kernel='cluster_disk')
    A = basis.design_matrix(station_lat, station_lon, lln, pgw.Displacement.Vertical)
    sol = pgw.RegularizedLeastSquares(A, sigma, R).solve(d, alpha)

Subpackages: earth (Love numbers, observables), basis (spherical harmonics, Green's functions, Slepian),
forward, inversion, observations (GNSS, synthetic), mesh (grids, masks), utils; config holds the data root.
"""
__version__ = '0.9.0'

from src.config import FileTool
from src.utils.enums import Displacement, LLN_Data, LLN_variable, Frame
from src.utils.constants import EarthConstant
from src.earth.love_numbers import LoveNumber
from src.basis.base import Basis
from src.basis.spherical_harmonic.sh_basis import SphericalHarmonicBasis
from src.basis.green.gf_basis import GreenFunctionBasis
from src.basis.slepian import SlepianBasis
from src.forward.model import ForwardModel
from src.inversion.solver import RegularizedLeastSquares
from src.inversion import regularization
from src.preprocessing import GNSSNetwork, GNSSDataset, read_tenv3, read_ngl_steps
from src.inversion.workflow import MonthlyInversion, InversionResult
from src.observations.grace import Mascon

__all__ = ['FileTool', 'Displacement', 'LLN_Data', 'LLN_variable', 'Frame', 'EarthConstant', 'LoveNumber', 'Basis',
           'SphericalHarmonicBasis', 'GreenFunctionBasis', 'SlepianBasis', 'ForwardModel', 'RegularizedLeastSquares',
           'regularization', 'GNSSNetwork', 'GNSSDataset', 'read_tenv3', 'read_ngl_steps', 'MonthlyInversion', 'InversionResult', 'Mascon',
           '__version__']

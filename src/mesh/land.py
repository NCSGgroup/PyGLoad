"""
Land fraction of grid cells from the 0.1-degree land mask (<data root>/input/land_mask/land_mask_res0.1.h5; rows from
90N to 90S, columns from 180W to 180E; 1 = land).

    from src.mesh.land import land_fraction
    frac = land_fraction(lat1, lon1, 0.5)          # (nlat, nlon) in [0, 1] for a cell-centred grid
    cells = GreenFunctionBasis.from_grid(lat1, lon1, 0.5, mask=frac >= 0.3)
"""
import numpy as np

from src.config import FileTool

_MASK = None


def _mask():
    global _MASK
    if _MASK is None:
        import h5py
        with h5py.File(FileTool.get_input_dir('land_mask/land_mask_res0.1.h5'), 'r') as f:
            _MASK = np.asarray(f['mask'][:-1, :], dtype=np.float32)     # (1800, 3600), row i: 90 - 0.1 (i + 0.5)
    return _MASK


def land_fraction(lat, lon, resolution):
    """
    Fraction of land in each cell of a cell-centred regular grid.
    :param lat, lon: 1d cell centres [deg]
    :param resolution: cell size [deg] (a multiple of 0.1 is exact; otherwise the nearest 0.1-degree pixels)
    :return: (nlat, nlon)
    """
    m = _mask()
    h = resolution / 2
    out = np.zeros((len(lat), len(lon)))
    for i, la in enumerate(lat):
        r0 = int(np.clip(round((90 - (la + h)) / 0.1), 0, 1800))
        r1 = int(np.clip(round((90 - (la - h)) / 0.1), r0 + 1, 1800))
        for j, lo in enumerate(lon):
            c0 = int(round(((lo - h + 180) % 360) / 0.1))
            n = max(int(round(resolution / 0.1)), 1)
            cols = (c0 + np.arange(n)) % 3600
            out[i, j] = m[r0:r1][:, cols].mean()
    return out

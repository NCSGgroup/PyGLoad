"""
Settings of the GNSS-TWS validation study: monthly vertical GNSS displacements at 4,807 stations (2003-2016) compared
with the elastic loading predicted from monthly terrestrial water storage (TWS) grids. Edit here, or pass the options
of run.py on the command line.

Folders (all under the PyGLoad data root, see doc/installation.md; another root with the option --data-root):

    input/LLN/                       PREM load Love numbers (Wang et al., 2012): PREM-LLNs-complete.dat
    input/gnss/NGL/                  NGL files (download.py --gnss):
        DataHoldings.txt, steps.txt      station catalogue and the offsets file used for the cleaning
        tenv3/IGS20/, tenv3/IGS14/       daily .tenv3 series of the stations, cut after 2016-12-31
        cmcf/                            NGL CM -> CF translation per station (cmcf/<site>)
        ntal.cf/, ntol.cf/               ESMGFZ non-tidal atmosphere / ocean loading per station, CF frame
    input/validation/                land_domain_0.5deg.nc (common 0.5-degree land domain, 0/1), common_months.txt,
                                     stations.csv (candidate stations, selection, region), station_exclusions.csv,
                                     ice_exclusions.csv (glacier basin units and ice fields of the station selection)
    input/mascon/                    CSR and JPL mascon files (download.py)
    input/tws_products/<name>/       other public products: GLWS 2.0, GLWS 4.0, GLDAS-2.2 (download.py)
    processed/GNSS_TWS_validation/   cleaned monthly GNSS series per frame: <frame>/gnss_vertical_monthly.nc (and .csv)
    cache/GNSS_TWS_validation/       product anomaly grids, spherical-harmonic coefficients, model displacements
    output/GNSS_TWS_validation/      scores and summaries per frame: <frame>/
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repository root, so that `src` is importable
from src.config import FileTool

STUDY = 'GNSS_TWS_validation'

# GNSS frames. NGL states that its series are in the CM frame of the reference frame; the loading model is in CF.
FRAMES = {
    'igs20cf': dict(tenv3='IGS20', cmcf=True,
                    label='IGS20 minus the NGL cmcf translation (CM -> CF), minus NTAL + NTOL (CF); main'),
    'igs14': dict(tenv3='IGS14', cmcf=False,
                  label='IGS14 taken as CF, minus NTAL + NTOL (CF); sensitivity'),
}

# public TWS products: reader, file or folder (relative to the data root), description
PRODUCTS = {
    'CSR': dict(kind='csr', path='input/mascon/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc',
                label='CSR RL06.3 mascons'),
    'JPL': dict(kind='jpl', path='input/mascon', pattern='GRCTellus.JPL.*.GLO.RL06.3M.MSCNv04CRI.nc',
                label='JPL RL06.3M mascons, CRI filtered'),
    'GLWS2': dict(kind='glws', path='input/tws_products/GLWS2/GLWS_2_0_2003-2019.nc', label='GLWS 2.0'),
    'GLWS4': dict(kind='glws', path='input/tws_products/GLWS4/GLWS_4_0_2003-2019.nc', label='GLWS 4.0'),
    'GLDAS22': dict(kind='gldas', path='input/tws_products/GLDAS22', pattern='GLDAS_CLSM025_DA1_D.A*.nc4',
                    label='GLDAS-2.2 CLSM with GRACE DA (daily 0.25 deg -> monthly 0.5 deg)'),
}

CONFIG = dict(
    frames=['igs20cf', 'igs14'],
    products=['CSR', 'JPL'],
    gnss='processed',                   # 'processed': read processed/<frame>/; 'raw': clean the NGL daily series
    # GNSS cleaning (the settings of the published series)
    period=(2003.0, 2017.0),            # daily epochs used: 2003.0 <= NGL decimal year < 2017.0
    min_daily=100,                      # at least this many daily solutions in the period
    outlier=5.0,                        # outlier rejection: |residual| > outlier x 1.4826 MAD (all components)
    min_days=10,                        # days needed for a monthly mean
    load_models=('ntal', 'ntol'),       # ESMGFZ loading (CF) subtracted, as the mean over the days of each month
    # forward model
    lmax=360,                           # spherical-harmonic degree (= 180 / 0.5 deg)
    love='PREM',                        # load Love numbers (Wang et al., 2012), converted to the CF frame
    sh_chunk=24,                        # months per spherical-harmonic analysis call
    station_block=190,                  # stations per design-matrix block
    # scores
    n_boot=2000,                        # station-bootstrap draws for the 5-95 % range of the medians
    seed=270,                           # bootstrap seed (each group / product / score has its own stream)
    workers=min(8, os.cpu_count() or 1),
)


def set_data_root(path):
    """use another data root than installation/pygload_paths.cfg (same as setting PYGLOAD_DATA)"""
    os.environ['PYGLOAD_DATA'] = str(Path(path).expanduser().resolve())


def data_path(rel):
    return FileTool.get_data_dir() / rel


def validation_dir(sub=None):
    return FileTool.get_input_dir('validation' if sub is None else 'validation/%s' % sub)


def ngl_dir(sub=None):
    return FileTool.get_input_dir('gnss/NGL' if sub is None else 'gnss/NGL/%s' % sub)


def processed_dir(frame):
    return FileTool.get_data_dir('processed') / STUDY / frame


def cache_dir():
    d = FileTool.get_cache_dir(STUDY)
    d.mkdir(parents=True, exist_ok=True)
    return d


def output_dir(frame=None):
    d = FileTool.get_output_dir(STUDY) if frame is None else FileTool.get_output_dir(STUDY) / frame
    d.mkdir(parents=True, exist_ok=True)
    return d

"""
Settings of the GNSS-TWS validation study: monthly vertical GNSS displacements at about 4800 stations (2003-2016)
compared with the elastic loading predicted from global terrestrial water storage (TWS) products. Edit here, then run
download.py (once, with internet), run.py and verify.py.

Folders (all under the PyGWater data root, see doc/installation.md):

    <data root>/input/GNSS_TWS_validation/   reference files of the study (station list, cleaned monthly GNSS series,
                                             land domain, reference scores); distributed with the study data set
    <data root>/input/tws_products/          the TWS products, one folder each (download.py)
    <data root>/input/gnss/NGL/              NGL series, steps file and load predictions (only to rebuild the GNSS
                                             series from the daily data, see GNSS_SOURCE)
    <data root>/cache/GNSS_TWS_validation/   product grids and spherical-harmonic coefficients (regenerable)
    <data root>/output/GNSS_TWS_validation/  scores and summaries

Environment variables:
    GNSS_TWS_PRODUCTS   comma-separated subset of PRODUCTS, e.g. CSR,GLWS4 (default: all)
    GNSS_TWS_GNSS       'reference' (default): the cleaned monthly series of the reference files;
                        'rebuild': rebuild them from the NGL daily series with pygwater (needs download.py --gnss)
    GNSS_TWS_WORKERS    processes used to read the daily GLDAS files (default 4)
    GNSS_TWS_<NAME>_PATH  use a copy of product <NAME> kept elsewhere instead of input/tws_products/ (a file, or the
                        folder for JPL and GLDAS22), e.g. GNSS_TWS_GLDAS22_PATH=/data/GLDAS/daily
    GNSS_TWS_NGL_DIR    NGL folder used by GNSS_TWS_GNSS=rebuild instead of input/gnss/NGL/
"""
import os

from pygwater.config import FileTool

STUDY = 'GNSS_TWS_validation'

PRODUCTS = {
    # name: reader kind, file or folder under <data root>/input/tws_products/, description
    'CSR': dict(kind='csr', path='CSR/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc',
                label='CSR RL06.3 mascons (GRACE alone)'),
    'JPL': dict(kind='jpl', path='JPL', pattern='GRCTellus.JPL.*.GLO.RL06.3M.MSCNv04CRI.nc',
                label='JPL RL06.3M mascons, CRI filtered (GRACE alone)'),
    'GLWS2': dict(kind='glws', path='GLWS2/GLWS_2_0_2003-2019.nc', label='GLWS 2.0 (WaterGAP + GRACE DA)'),
    'GLWS4': dict(kind='glws', path='GLWS4/GLWS_4_0_2003-2019.nc', label='GLWS 4.0 (WaterGAP + GRACE DA)'),
    'GLDAS22': dict(kind='gldas', path='GLDAS22', pattern='GLDAS_CLSM025_DA1_D.A*.nc4',
                    label='GLDAS-2.2 CLSM GRACE DA (daily 0.25 deg -> monthly 0.5 deg)'),
}

_sel = os.environ.get('GNSS_TWS_PRODUCTS', '')
CONFIG = dict(
    products=[p.strip() for p in _sel.split(',') if p.strip()] or list(PRODUCTS),
    gnss_source=os.environ.get('GNSS_TWS_GNSS', 'reference'),
    period=(2003.0, 2017.0),            # daily GNSS epochs used: 2003-01-01 <= t < 2017-01-01 (decimal years)
    # GNSS cleaning (only for gnss_source='rebuild'; the reference series were made with these settings)
    outlier=5.0,                        # outlier rejection: |residual| > outlier x 1.4826 MAD
    min_span=2.0,                       # drop series shorter than this [years] (first to last day)
    min_daily=100,                      # at least this many daily solutions in the period
    min_days=10,                        # days needed for a monthly mean
    min_coverage=0.6,                   # station selection: share of the common months with a GNSS value
    min_span_months=5.0,                # station selection: first to last common month [years]
    load_models=('ntal', 'ntol'),       # ESMGFZ loading (CF frame, NGL loadpredictions/<model>.cf/) subtracted
    # forward model
    lmax=360,                           # spherical-harmonic degree (= 180 / 0.5 deg grid)
    love='PREM',                        # load Love numbers (Wang et al., 2012), converted to the CF frame
    sh_chunk=24,                        # months per spherical-harmonic analysis call
    station_block=190,                  # stations per design-matrix block (~200 MB)
    # scores
    n_boot=2000,                        # station bootstrap draws for the 5-95 % range of the medians
    seed=270,                           # random seed of the bootstrap
    boundary_km=50.0,                   # stations closer than this to a basin-unit border form the 'boundary' class
    workers=int(os.environ.get('GNSS_TWS_WORKERS', '4')),
)

for p in CONFIG['products']:
    if p not in PRODUCTS:
        raise ValueError('unknown product %s; known: %s' % (p, ', '.join(PRODUCTS)))
if CONFIG['gnss_source'] not in ('reference', 'rebuild'):
    raise ValueError("GNSS_TWS_GNSS must be 'reference' or 'rebuild'")


def reference_dir(sub=None):
    return FileTool.get_input_dir(STUDY if sub is None else '%s/%s' % (STUDY, sub))


def product_path(name):
    """file (or folder for JPL, GLDAS22) of a product: GNSS_TWS_<NAME>_PATH or <data root>/input/tws_products/..."""
    from pathlib import Path
    env = os.environ.get('GNSS_TWS_%s_PATH' % name)
    return Path(env).expanduser() if env else FileTool.get_input_dir('tws_products/%s' % PRODUCTS[name]['path'])


def ngl_dir(sub=None):
    from pathlib import Path
    env = os.environ.get('GNSS_TWS_NGL_DIR')
    root = Path(env).expanduser() if env else FileTool.get_input_dir('gnss/NGL')
    return root if sub is None else root / sub


def cache_dir():
    d = FileTool.get_cache_dir(STUDY)
    d.mkdir(parents=True, exist_ok=True)
    return d


def output_dir():
    d = FileTool.get_output_dir(STUDY)
    d.mkdir(parents=True, exist_ok=True)
    return d

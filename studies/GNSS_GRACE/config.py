"""
Settings of the GNSS-GRACE study; edit here, then run download.py (once, with internet), run.py and verify.py.
The region is chosen with the environment variable GNSS_GRACE_REGION (default PNW), e.g.

    GNSS_GRACE_REGION=CA python studies/GNSS_GRACE/run.py
"""
import os

REGIONS = {
    # name: stations box (lat_min, lat_max, lon_min, lon_max) [deg]; output in <data root>/output/GNSS_GRACE/<name>/
    'PNW': (41.0, 50.0, -125.0, -116.0),        # Washington and Oregon (Fu et al., 2015)
    'CA': (35.0, 42.0, -124.0, -117.0),         # northern and central California (Argus et al., 2014)
    'HARVEY': (27.0, 33.0, -99.0, -91.0),       # south Texas and Louisiana, Hurricane Harvey 2017 (Milliner et al., 2018)
}
# settings that differ from the defaults below, per region
OVERRIDES = {
    'HARVEY': dict(step='day', start='2012-01', end='2021-12', resolution=0.25,
                   reference=('2012-01', '2021-12'),                  # GRACE anomalies (no GRACE in Aug-Oct 2017)
                   gnss_reference=('2017-08-01', '2017-08-24'),       # GNSS: change since the days before landfall
                   window=('2017-08-01', '2017-10-10'),               # days inverted (Milliner et al., 2018)
                   buffer=1.0,                                        # dense network, local event: cells 1 deg beyond it
                   # common mode stacked from the stations NW of the rain, as Milliner et al. (2018) did. The regional
                   # level depends on this choice (peak 28 Gt with it, 62 Gt with a ring of stations 700-1300 km away:
                   # cme_reference=dict(center=(29.8, -95.0), min_km=700, max_km=1300)), because every reference
                   # carries its own slow signal; see the paper's Harvey section
                   cme_reference=(31.0, 33.5, -99.0, -97.0),
                   common_mode=()),                                   # (stacking replaces the inversion's offsets)
}
_region = os.environ.get('GNSS_GRACE_REGION', 'PNW')
_components = os.environ.get('GNSS_GRACE_COMPONENTS', 'U')         # 'U' or 'UNE'
_step = os.environ.get('GNSS_GRACE_STEP', OVERRIDES.get(_region, {}).get('step', 'month'))   # 'month' or 'day'
if _region not in REGIONS:
    raise ValueError('GNSS_GRACE_REGION=%s; known regions: %s' % (_region, ', '.join(REGIONS)))

CONFIG = dict(
    name=_region + ('' if _components == 'U' else '_' + _components) + ('' if _step == 'month' else '_daily'),
                                                # output folder, e.g. PNW, PNW_UNE, PNW_daily
    box=REGIONS[_region],
    start='2006-01', end='2022-12',             # period
    min_years=4.0,                              # at least 4 years of daily data (Schumacher et al., 2018)
    frame='IGS20',                              # NGL archive
    reference=('2006-01', '2022-12'),           # common reference period of GNSS and GRACE anomalies (the whole
                                                # period, so stations installed late are kept)
    components=tuple(_components),              # ('U',) or ('U', 'N', 'E')
    step=_step,                                 # 'month' (monthly means) or 'day' (daily epochs)
    common_mode=('N', 'E'),                     # one network-wide offset per month for these components (used
                                                # only if they are in components): GNSS horizontals carry a
                                                # seasonal common mode that is not loading
    resolution=0.5,                             # cell size of the inversion [deg]
    buffer=2.0,                                 # cells extend this far beyond the station box [deg]
    min_land=0.3,                               # cells with less land are left out (0: keep ocean cells)
    far_field=12.0,                             # subtract the displacement of GRACE load outside the cells, within
                                                # this many degrees of the box (0: none); needs grace_file
    kernel='disk',                              # 'point_mass', 'disk' or 'cluster_disk'
    lmax=13000,                                 # Love numbers (Wang et al., 2012), CF frame
    alpha_method='gcv',                         # 'gcv', 'discrepancy', 'lcurve'
    nonnegative=False,                          # anomalies can be negative
    load_models=('NTAL', 'NTOL'),               # NGL/ESMGFZ load predictions removed from the series (atmosphere,
                                                # non-tidal ocean); () keeps them. Needs download.py's loadpred/ files
    cm_to_cf='grace',                           # series in CM, Love numbers in CF: 'grace' converts with the geocentre
                                                # motion of the global grace_file (consistent with the comparison),
                                                # 'ngl' with NGL's translation (needs download.py's loadpred_cmcf/
                                                # files; Altamimi et al. 2023), False keeps CM
    spatial_filter_km=500,                      # neighbour test radius (Schumacher et al., 2018)
    grace_file='mascon/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc',   # under <data root>/input/
)
CONFIG.update({k: v for k, v in OVERRIDES.get(_region, {}).items() if k != 'step'})

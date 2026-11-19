"""
Download the inputs of the GNSS-TWS validation. Run on a computer with internet access:

    python studies/GNSS_TWS_validation/download.py            # TWS products
    python studies/GNSS_TWS_validation/download.py --gnss     # also the NGL daily series and load predictions
    python studies/GNSS_TWS_validation/download.py --list     # only print the sources, versions and DOIs

Files already present are skipped, so the script can be re-run to complete a download.

TWS products go to <data root>/input/tws_products/<name>/. CSR, GLWS 2.0 and GLWS 4.0 are open downloads. JPL and
GLDAS-2.2 need a free NASA Earthdata login; the script does not handle credentials and prints what to fetch and where
to put it instead.

--gnss is needed only to rebuild the cleaned GNSS series from the daily data (run.py with GNSS_TWS_GNSS=rebuild); the
reference files of the study already hold the cleaned monthly series. It fetches, for the stations used
(<data root>/input/GNSS_TWS_validation/stations/stations.csv), the NGL daily IGS14 .tenv3 series, the steps file and the
ESMGFZ non-tidal atmosphere and ocean loading at the stations in the CF frame (NGL loadpredictions/{ntal,ntol}.cf),
about 7 GB in total. Please keep the download polite (it is sequential, with a pause between files).
"""
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

from pygwater.preprocessing import ngl

from config import CONFIG, PRODUCTS, reference_dir, product_path, ngl_dir

USER_AGENT = 'PyGWater GNSS_TWS_validation study'

SOURCES = {
    'CSR': dict(
        url='https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/'
            'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc',
        version='CSR GRACE/GRACE-FO RL06.3 mascons, all corrections (file updated monthly)',
        cite='Save, H., Bettadpur, S., & Tapley, B. D. (2016). High-resolution CSR GRACE RL05 mascons. J. Geophys. Res. '
             'Solid Earth, 121, 7547-7569, doi:10.1002/2016JB013007; data: doi:10.15781/cgq9-nh24'),
    'JPL': dict(
        url=None,
        version='JPL GRACE/GRACE-FO mascon RL06.3 version 04, CRI filtered (PO.DAAC short name '
                'TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4), file GRCTellus.JPL.<span>.GLO.RL06.3M.MSCNv04CRI.nc',
        login='NASA Earthdata (PO.DAAC): https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4',
        cite='Watkins, M. M., et al. (2015), J. Geophys. Res. Solid Earth, 120, 2648-2671, doi:10.1002/2014JB011547; '
             'Wiese, D. N., et al. (2016), Water Resour. Res., 52, 7490-7502, doi:10.1002/2016WR019344; '
             'data: doi:10.5067/TEMSC-3JC634'),
    'GLWS2': dict(
        url='https://download.pangaea.de/dataset/954742/files/GLWS_2_0_2003-2019.nc', size=1692094275,
        version='GLWS 2.0, ensemble mean, 01/2003-12/2019, 0.5 deg',
        cite='Gerdener, H., Schulze, K., & Kusche, J. (2023). GLWS 2.0 [dataset]. PANGAEA, '
             'doi:10.1594/PANGAEA.954742; Gerdener et al. (2023), J. Geod., 97, 73, doi:10.1007/s00190-023-01763-9'),
    'GLWS4': dict(
        url='https://download.pangaea.de/dataset/993748/files/GLWS_4_0_2003-2019.nc', size=1692093862,
        version='GLWS 4.0, ensemble mean, 01/2003-12/2019, 0.5 deg',
        cite='Gerdener, H., & Kusche, J. (2026). GLWS 4.0 [dataset]. PANGAEA, doi:10.1594/PANGAEA.993748'),
    'GLDAS22': dict(
        url=None,
        version='GLDAS-2.2 CLSM 0.25 deg daily with GRACE DA (GLDAS_CLSM025_DA1_D, version 2.2), '
                'variables TWS_tavg and GWS_tavg, 2003-02-01 .. 2016-12-31 (one file per day)',
        login='NASA Earthdata (GES DISC): https://disc.gsfc.nasa.gov/datasets/GLDAS_CLSM025_DA1_D_2.2/summary; '
              'the subsetter output (TWS_tavg, GWS_tavg) is enough',
        cite='Li, B., et al. (2019). Global GRACE data assimilation for groundwater and drought monitoring. Water '
             'Resour. Res., 55, 7564-7586, doi:10.1029/2018WR024618; data: doi:10.5067/TXBMLX370XX8'),
}
GNSS_SOURCES = (
    'NGL daily IGS14 .tenv3 series and steps.txt, https://geodesy.unr.edu (Blewitt, G., Hammond, W. C., & Kreemer, C. '
    '(2018). Harnessing the GPS data explosion for interdisciplinary science. Eos, 99, doi:10.1029/2018EO104623)',
    'ESMGFZ non-tidal atmospheric and oceanic loading (NTAL, NTOL), CF frame, at the NGL stations '
    '(https://geodesy.unr.edu/gps_timeseries/loadpredictions/; Dill, R., & Dobslaw, H. (2013), J. Geophys. Res. '
    'Solid Earth, 118, 5008-5017, doi:10.1002/jgrb.50353)')
LOADPRED_CF_URL = 'https://geodesy.unr.edu/gps_timeseries/loadpredictions/{model}.cf/{site_lower}'


def fetch(url, path, size=None, pause=0.0):
    """stream url to path (via path.part); skip if present (and of the expected size)"""
    path = Path(path)
    if path.exists() and (size is None or path.stat().st_size == size):
        return 'present'
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + '.part')
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r, open(part, 'wb') as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    if size is not None and part.stat().st_size != size:
        raise IOError('%s: %d bytes, expected %d' % (url, part.stat().st_size, size))
    part.rename(path)
    time.sleep(pause)
    return 'ok'


def list_sources():
    for name, s in SOURCES.items():
        print('%s: %s\n  %s\n  cite: %s' % (name, PRODUCTS[name]['label'], s['version'], s['cite']))
        print('  url: %s' % (s['url'] or 'login needed: ' + s['login']))
    print('GNSS:')
    for s in GNSS_SOURCES:
        print('  ' + s)


def products():
    for name in CONFIG['products']:
        s, p = SOURCES[name], PRODUCTS[name]
        target = product_path(name)
        if s['url'] is None:
            have = sorted(target.glob(p['pattern'])) if target.is_dir() else []
            if have:
                print('%s: present (%d file%s in %s)' % (name, len(have), '' if len(have) == 1 else 's', target))
            else:
                print('%s: NOT present. Download by hand (%s)\n  %s\n  and put the file(s) into %s'
                      % (name, s['login'], s['version'], target))
            continue
        try:
            print('%s: %s' % (name, fetch(s['url'], target, s.get('size'))))
        except Exception as e:
            print('%s: FAILED (%s); url %s' % (name, e, s['url']))


def gnss():
    stations = pd.read_csv(reference_dir('stations/stations.csv'), index_col=0)
    sites = list(stations.index[stations.used.astype(bool)])
    out = ngl_dir()
    out.mkdir(parents=True, exist_ok=True)
    if not (out / 'steps.txt').exists() or not (out / 'DataHoldings.txt').exists():
        ngl.download_catalogue(out)
    print('%d stations' % len(sites))
    status = ngl.download_series(sites, out / 'tenv3_IGS14', frame='IGS14')
    report('IGS14 series', status)
    for model in CONFIG['load_models']:
        status = {}
        for s in sites:
            try:
                status[s] = fetch(LOADPRED_CF_URL.format(model=model, site_lower=s.lower()),
                                  out / ('loadpred_cf/%s/%s' % (model, s.lower())), pause=0.2)
            except Exception as e:
                status[s] = str(e)
        report('%s (CF) load predictions' % model.upper(), status)


def report(what, status):
    failed = {k: v for k, v in status.items() if v not in ('ok', 'present')}
    print('%s: downloaded %d, already present %d, failed %d' % (what, sum(v == 'ok' for v in status.values()),
                                                              sum(v == 'present' for v in status.values()),
                                                              len(failed)))
    for k, v in list(failed.items())[:10]:
        print('  %s: %s' % (k, v))


if __name__ == '__main__':
    if '--list' in sys.argv:
        list_sources()
        sys.exit()
    products()
    if '--gnss' in sys.argv:
        gnss()

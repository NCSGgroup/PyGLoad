"""
Download the inputs of the GNSS-TWS validation (needs internet). Files already present are skipped.

    python studies/GNSS_TWS_validation/download.py                  # CSR mascons; JPL: instructions (login)
    python studies/GNSS_TWS_validation/download.py --products CSR,JPL,GLWS2,GLWS4,GLDAS22
    python studies/GNSS_TWS_validation/download.py --gnss           # NGL series, cmcf, NTAL / NTOL of the stations
    python studies/GNSS_TWS_validation/download.py --check          # md5 of the product files against those we used
    python studies/GNSS_TWS_validation/download.py --list           # sources, versions, citations

Products go to <data root>/input/mascon/ (CSR, JPL) and <data root>/input/tws_products/<name>/ (GLWS 2.0 / 4.0,
GLDAS-2.2). CSR, GLWS 2.0 and GLWS 4.0 are open downloads. JPL and GLDAS-2.2 need a free NASA Earthdata login: the
script does not handle credentials; it prints what to fetch and where to put it.

--gnss fetches, for the stations used (<data root>/input/validation/stations.csv), the NGL daily .tenv3 series in IGS20
and IGS14, the NGL cmcf translation and the ESMGFZ non-tidal atmosphere and ocean loading in the CF frame (NGL
loadpredictions/{ntal,ntol}.cf), each cut after 2016-12-31 with HTTP Range requests (about 14 GB), plus DataHoldings.txt
and steps.txt. NGL reprocesses its series and extends steps.txt over time, so files fetched today can differ from the
ones of the study; to reproduce the published numbers use the files of the study's data package. Please keep the
download polite: it is sequential, with a pause between requests.
"""
import argparse
import datetime as dt
import hashlib
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

import config
from config import PRODUCTS

USER_AGENT = 'PyGLoad GNSS_TWS_validation study'
NGL = 'https://geodesy.unr.edu'

# the files we used (md5); CSR replaces its file every month: the November 2026 file has the same 2003-2016 fields
SOURCES = {
    'CSR': dict(
        url='https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/'
            'CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc',
        md5={'ad180503c3d178281b0d006571474d5c': 'file of October 2024 (102,627,341 bytes), used in the paper',
             '5a837be453345ff5d971fa3cd7064c25': 'file of 9 November 2026 (113,954,516 bytes), same 2003-2016 fields'},
        version='CSR GRACE/GRACE-FO RL06.3 mascons, all corrections, 0.25 deg (file replaced monthly)',
        cite='Save, H., Bettadpur, S., & Tapley, B. D. (2016). High-resolution CSR GRACE RL05 mascons. J. Geophys. Res. '
             'Solid Earth, 121, 7547-7569, doi:10.1002/2016JB013007; data: doi:10.15781/cgq9-nh24'),
    'JPL': dict(
        url=None, file='GRCTellus.JPL.200204_202508.GLO.RL06.3M.MSCNv04CRI.nc',
        md5={'9348bb27dc4392ff24c99e1078002742': 'file used in the paper (44,041,869 bytes, 2002-04 .. 2025-08)'},
        version='JPL GRACE/GRACE-FO mascons RL06.3 version 04, CRI filtered, 0.5 deg (PO.DAAC short name '
                'TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4)',
        login='NASA Earthdata login (free): https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4',
        cite='Watkins, M. M., et al. (2015), J. Geophys. Res. Solid Earth, 120, 2648-2671, doi:10.1002/2014JB011547; '
             'Wiese, D. N., et al. (2016), Water Resour. Res., 52, 7490-7502, doi:10.1002/2016WR019344; '
             'data: doi:10.5067/TEMSC-3JC634'),
    'GLWS2': dict(
        url='https://download.pangaea.de/dataset/954742/files/GLWS_2_0_2003-2019.nc',
        md5={'c30b9605615be04b4c22fab5307da982': 'file used in the paper (1,692,094,275 bytes)'},
        version='GLWS 2.0, ensemble mean, 01/2003-12/2019, 0.5 deg',
        cite='Gerdener, H., Schulze, K., & Kusche, J. (2023). GLWS 2.0 [dataset]. PANGAEA, doi:10.1594/PANGAEA.954742; '
             'Gerdener, H., Kusche, J., Schulze, K., Doell, P., & Klos, A. (2023), J. Geod., 97, 73, '
             'doi:10.1007/s00190-023-01763-9'),
    'GLWS4': dict(
        url='https://download.pangaea.de/dataset/993748/files/GLWS_4_0_2003-2019.nc',
        md5={'5b52ee9b96bb9ef5baf944584839ab67': 'file used in the paper (1,692,093,862 bytes)'},
        version='GLWS 4.0, ensemble mean, 01/2003-12/2019, 0.5 deg',
        cite='GLWS 4.0 [dataset]. PANGAEA, doi:10.1594/PANGAEA.993748'),
    'GLDAS22': dict(
        url=None,
        version='GLDAS-2.2 CLSM 0.25 deg daily with GRACE DA (GLDAS_CLSM025_DA1_D, version 2.2), variable TWS_tavg, '
                '2003-02-01 .. 2016-12-31, one file per day (the subsetter output with TWS_tavg is enough)',
        login='NASA Earthdata login (free): https://disc.gsfc.nasa.gov/datasets/GLDAS_CLSM025_DA1_D_2.2/summary',
        cite='Li, B., et al. (2019). Global GRACE data assimilation for groundwater and drought monitoring. Water '
             'Resour. Res., 55, 7564-7586, doi:10.1029/2018WR024618; data: doi:10.5067/TXBMLX370XX8'),
}
GNSS_SOURCES = (
    'NGL daily .tenv3 series (IGS20 and IGS14), DataHoldings.txt, steps.txt: %s (Blewitt, G., Hammond, W. C., & '
    'Kreemer, C. (2018). Harnessing the GPS data explosion for interdisciplinary science. Eos, 99, '
    'doi:10.1029/2018EO104623)' % NGL,
    'NGL cmcf: CM -> CF translation at the stations from the ITRF2020 geocentre motion model (Altamimi, Z., Rebischung, '
    'P., Collilieux, X., Metivier, L., & Chanard, K. (2023). ITRF2020: an augmented reference frame refining the '
    'modeling of nonlinear station motions. J. Geod., 97, 47, doi:10.1007/s00190-023-01738-w)',
    'ESMGFZ non-tidal atmospheric and oceanic loading (NTAL, NTOL), CF frame, at the stations as distributed by NGL '
    '(%s/gps_timeseries/loadpredictions/; Dill, R., & Dobslaw, H. (2013). Numerical simulations of global-scale '
    'high-resolution hydrological crustal deformations. J. Geophys. Res. Solid Earth, 118, 5008-5017, '
    'doi:10.1002/jgrb.50353)' % NGL)

# NGL station files: url, folder under input/gnss/NGL/, file name, fixed line length, header length
GNSS_FILES = {
    'IGS20': (NGL + '/gps_timeseries/IGS20/tenv3/IGS20/{site}.tenv3', 'tenv3/IGS20', '{site}.tenv3', 204, 204),
    'IGS14': (NGL + '/gps_timeseries/IGS14/tenv3/IGS14/{site}.tenv3', 'tenv3/IGS14', '{site}.tenv3', 204, 204),
    'cmcf': (NGL + '/gps_timeseries/loadpredictions/cmcf/{low}', 'cmcf', '{low}', 58, 83),
    'ntal.cf': (NGL + '/gps_timeseries/loadpredictions/ntal.cf/{low}', 'ntal.cf', '{low}', 58, 83),
    'ntol.cf': (NGL + '/gps_timeseries/loadpredictions/ntol.cf/{low}', 'ntol.cf', '{low}', 58, 83),
}
END = dt.date(2017, 1, 1)
MON = {m: i + 1 for i, m in enumerate('JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC'.split())}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _get(url, a=None, z=None):
    h = {'User-Agent': USER_AGENT}
    if a is not None:
        h['Range'] = 'bytes=%d-%d' % (a, z)
    with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=120) as r:
        return r.read(), r.status


def fetch(url, path):
    """whole file, via path.part"""
    path = Path(path)
    if path.exists():
        return 'present'
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + '.part')
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r, open(part, 'wb') as f:
        for chunk in iter(lambda: r.read(1 << 20), b''):
            f.write(chunk)
    part.rename(path)
    return 'ok'


def _line_date(line):
    s = line.split()[1]
    yy = int(s[:2])
    return dt.date(2000 + yy if yy < 50 else 1900 + yy, MON[s[2:5].upper()], int(s[5:7]))


def fetch_until_2016(url, path, start, line, header, pause=0.3):
    """fixed-width NGL station file from its start to the first line after 2016-12-31 (HTTP Range), cut at a complete
    line: the start of the series up to 2016-12-31 plus a few days, not the full series"""
    path = Path(path)
    if path.exists():
        return 'present'
    data = b''
    z = header + line * ((END - start).days + 3) - 1
    for _ in range(6):
        chunk, status = _get(url, len(data), z)
        data += chunk
        time.sleep(pause)
        if status == 200 or len(chunk) < z - (len(data) - len(chunk)) + 1:
            break                                                   # whole file or end of file reached
        lines = data[:data.rfind(b'\n') + 1].decode().splitlines()
        if len(lines) > 1 and _line_date(lines[-1]) >= END:
            break
        z = len(data) + line * ((END - _line_date(lines[-1])).days + 3) - 1 if len(lines) > 1 else len(data) + line * 400
    data = data[:data.rfind(b'\n') + 1]
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + '.part')
    part.write_bytes(data)
    part.rename(path)
    return 'ok'


def gnss():
    out = config.ngl_dir()
    out.mkdir(parents=True, exist_ok=True)
    for name in ('DataHoldings.txt', 'steps.txt'):
        print('%s: %s' % (name, fetch(NGL + '/NGLStationPages/' + name, out / name)))
    S = pd.read_csv(config.validation_dir('stations.csv'), index_col=0)
    sites = list(S.index[S.used.astype(bool)])
    H = pd.read_csv(out / 'DataHoldings.txt', sep=r'\s+', usecols=[0, 7], skiprows=1, header=None,
                    names=['sta', 'beg']).set_index('sta')
    print('%d stations' % len(sites))
    for kind, (url, folder, fname, line, header) in GNSS_FILES.items():
        status = {}
        for i, s in enumerate(sites):
            beg = pd.Timestamp(H.beg.get(s, '2003-01-01')).date()
            try:
                status[s] = fetch_until_2016(url.format(site=s, low=s.lower()),
                                             out / folder / fname.format(site=s, low=s.lower()), beg, line, header)
            except urllib.error.HTTPError as e:
                status[s] = 'HTTP %d' % e.code
            except Exception as e:
                status[s] = str(e)
            if i % 500 == 0:
                print('  %s %d / %d' % (kind, i, len(sites)), flush=True)
        failed = {k: v for k, v in status.items() if v not in ('ok', 'present')}
        print('%s: downloaded %d, present %d, failed %d %s' % (kind, sum(v == 'ok' for v in status.values()),
                                                             sum(v == 'present' for v in status.values()), len(failed),
                                                             list(failed.items())[:5]))


def target(name):
    p = PRODUCTS[name]
    if 'pattern' in p:
        return config.data_path(p['path']) / SOURCES[name].get('file', '')
    return config.data_path(p['path'])


def products(names):
    for name in names:
        s = SOURCES[name]
        t = target(name)
        if s['url'] is None:
            have = sorted(config.data_path(PRODUCTS[name]['path']).glob(PRODUCTS[name]['pattern']))
            if have:
                print('%s: present (%d file%s)' % (name, len(have), '' if len(have) == 1 else 's'))
            else:
                print('%s: NOT present. Download by hand after logging in (%s):\n  %s%s\n  into %s'
                      % (name, s['login'], s['version'], ('\n  file: %s' % s['file']) if 'file' in s else '',
                         config.data_path(PRODUCTS[name]['path'])))
            continue
        try:
            print('%s: %s (%s)' % (name, fetch(s['url'], t), t))
        except Exception as e:
            print('%s: FAILED (%s); url %s' % (name, e, s['url']))


def check(names):
    ok = True
    for name in names:
        s = SOURCES[name]
        if 'md5' not in s:
            continue
        p = PRODUCTS[name]
        files = sorted(config.data_path(p['path']).glob(p['pattern'])) if 'pattern' in p else [target(name)]
        for f in files:
            if not f.exists():
                print('%s: %s missing' % (name, f)); ok = False
                continue
            h = md5(f)
            print('%s: %s md5 %s -> %s' % (name, f.name, h, s['md5'].get(h, 'NOT a file of the study')))
            ok &= h in s['md5']
    return ok


def list_sources():
    for name, s in SOURCES.items():
        print('%s: %s\n  %s\n  cite: %s' % (name, PRODUCTS[name]['label'], s['version'], s['cite']))
        print('  url: %s' % (s['url'] or 'login needed: ' + s['login']))
        for h, d in s.get('md5', {}).items():
            print('  md5 %s  %s' % (h, d))
    print('GNSS:')
    for s in GNSS_SOURCES:
        print('  ' + s)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='download the inputs of the GNSS-TWS validation')
    ap.add_argument('--data-root', help='PyGLoad data root (default: installation/pygload_paths.cfg or PYGLOAD_DATA)')
    ap.add_argument('--products', default='CSR,JPL', help='comma-separated: %s' % ','.join(SOURCES))
    ap.add_argument('--gnss', action='store_true', help='also the NGL / ESMGFZ station files (about 14 GB)')
    ap.add_argument('--check', action='store_true', help='only check the md5 of the product files')
    ap.add_argument('--list', action='store_true', help='only print the sources')
    a = ap.parse_args()
    if a.data_root:
        config.set_data_root(a.data_root)
    names = [n for n in a.products.split(',') if n]
    if a.list:
        list_sources()
        sys.exit()
    if a.check:
        sys.exit(0 if check(names) else 1)
    products(names)
    if a.gnss:
        gnss()

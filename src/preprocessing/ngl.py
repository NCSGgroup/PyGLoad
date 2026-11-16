"""
Nevada Geodetic Laboratory (NGL, University of Nevada, Reno; Blewitt et al., 2018): station list, station selection
and download of daily .tenv3 time series and the steps (offsets) file.

    from pygwater.preprocessing import ngl
    ngl.download_catalogue(out_dir)                                   # DataHoldings.txt and steps.txt
    st = ngl.read_holdings(out_dir / 'DataHoldings.txt')
    sites = ngl.select_stations(st, box=(42, 49, -125, -116), start='2006-01', end='2020-12', min_years=5)
    ngl.download_series(sites, out_dir / 'tenv3')                     # SITE.tenv3, skips files already present

Downloads need internet access to geodesy.unr.edu; run them on your own computer. Files go to
<data root>/input/gnss/NGL/ by default (see the study studies/GNSS_GRACE/download.py).

NGL series are processed in the IGS20 frame (IGS14 for the older archive). The URL of a station's .tenv3 file is
tried in the order of TENV3_URLS; the first that answers with a valid series is used. The IGS20 layout
(gps_timeseries/IGS20/tenv3/IGS20/SITE.tenv3) was checked in November 2026; the region folders next to it (NA/, EU/,
...) hold plate-fixed series, which are not used here.
"""
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

BASE = 'https://geodesy.unr.edu'
HOLDINGS_URL = BASE + '/NGLStationPages/DataHoldings.txt'
STEPS_URL = BASE + '/NGLStationPages/steps.txt'
LOAD_URL = BASE + '/gps_timeseries/tenv3_loadpredictions/{site}.tenv3'
# CM -> CF translation per station (Altamimi et al. 2023 geocentre motion): NGL series are in CM (IGS20); subtracting
# these values places them in CF, the frame of PyGWater's Love numbers for regional loading (README_loadpred.txt).
# File: lowercase station id, no extension; rows: station YYMMMDD decimal_year east north up [m]
CMCF_URL = BASE + '/gps_timeseries/loadpredictions/cmcf/{site_lower}'
# load predictions appended to the .tenv3 columns: east, north, up [m] of each model, ESMGFZ (atmosphere NTAL, non-tidal
# ocean NTOL, hydrology LSDM HYDL), GRACE mascons MASC and others; the columns are read from the header line
# (__ee_ntal __nn_ntal __uu_ntal ...); these positions are the default when a file has no header
LOAD_MODELS = {'NTAL': 23, 'NTOL': 26, 'HYDL': 29, 'MASC': 32}
TENV3_URLS = {
    'IGS20': (BASE + '/gps_timeseries/IGS20/tenv3/IGS20/{site}.tenv3',
              BASE + '/gps_timeseries/tenv3/IGS20/{site}.tenv3'),
    'IGS14': (BASE + '/gps_timeseries/IGS14/tenv3/IGS14/{site}.tenv3',
              BASE + '/gps_timeseries/tenv3/IGS14/{site}.tenv3'),
}


@dataclass
class StationInfo:
    name: str
    lat: float
    lon: float                      # -180..180
    height: float
    start: str                      # 'YYYY-MM-DD'
    end: str
    n_solutions: int


def read_holdings(path):
    """
    DataHoldings.txt: Sta Lat(deg) Long(deg) Hgt(m) X Y Z Dtbeg Dtend Dtmod NumSol [StaOrigName]
    :return: list of StationInfo
    """
    out = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 11 or p[0] == 'Sta':
                continue
            try:
                lon = (float(p[2]) + 180) % 360 - 180
                out.append(StationInfo(p[0], float(p[1]), lon, float(p[3]), p[7], p[8], int(p[10])))
            except ValueError:
                continue
    return out


def select_stations(stations, box=None, inside=None, start=None, end=None, min_years=3.0, min_solutions=0):
    """
    Stations in a region whose data cover at least `min_years` within [start, end].
    :param box: (lat_min, lat_max, lon_min, lon_max) [degree]
    :param inside: or a function (lat, lon) -> bool array (e.g. a shapefile test)
    :param start, end: 'YYYY-MM' or 'YYYY-MM-DD' (None = open)
    :return: list of station names
    """
    s0 = np.datetime64(start or '1990-01', 'D')
    s1 = np.datetime64(end or '2100-01', 'D')
    out = []
    for s in stations:
        if box is not None:
            lon = (s.lon - box[2]) % 360 + box[2]
            if not (box[0] <= s.lat <= box[1] and box[2] <= lon <= box[3]):
                continue
        if inside is not None and not bool(inside(np.array([s.lat]), np.array([s.lon]))[0]):
            continue
        a, b = max(np.datetime64(s.start, 'D'), s0), min(np.datetime64(s.end, 'D'), s1)
        if (b - a).astype(int) / 365.25 < min_years or s.n_solutions < min_solutions:
            continue
        out.append(s.name)
    return out


def select_reference_ring(stations, center, min_km, max_km, start=None, end=None, min_years=3.0, cell_deg=1.0):
    """
    Reference stations for common-mode removal, far from a load: stations between min_km and max_km from `center`
    with enough data, thinned to the one with most solutions per cell_deg x cell_deg cell (an even ring, few files).
    :param center: (lat, lon) [degree], e.g. the centre of a storm
    :return: list of station names
    """
    la0, lo0 = np.radians(center[0]), np.radians(center[1])
    s0 = np.datetime64(start or '1990-01', 'D')
    s1 = np.datetime64(end or '2100-01', 'D')
    best = {}
    for s in stations:
        la, lo = np.radians(s.lat), np.radians(s.lon)
        d = 6371.0 * np.arccos(np.clip(np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1, 1))
        if not (min_km <= d <= max_km):
            continue
        a, b = max(np.datetime64(s.start, 'D'), s0), min(np.datetime64(s.end, 'D'), s1)
        if (b - a).astype(int) / 365.25 < min_years or s.n_solutions < int(min_years * 365):
            continue
        key = (int(np.floor(s.lat / cell_deg)), int(np.floor(s.lon / cell_deg)))
        if key not in best or s.n_solutions > best[key].n_solutions:
            best[key] = s
    return sorted(v.name for v in best.values())


def _fetch(url, path, timeout=60, check=None):
    req = urllib.request.Request(url, headers={'User-Agent': 'PyGWater'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    if data[:2] == b'\x1f\x8b':                       # gzip-compressed
        import gzip
        data = gzip.decompress(data)
    if check is not None and not check(data):
        raise ValueError('not a valid file: %s' % url)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_bytes(data)
    return len(data)


def _is_tenv3(data):
    """a .tenv3 file starts with the header line 'site YYMMMDD ...' and has data lines"""
    head = data[:200].decode('ascii', 'replace')
    return head.startswith('site') and data.count(b'\n') > 1


def read_load_predictions(path, models=('NTAL', 'NTOL')):
    """
    Load predictions of one station from a NGL tenv3_loadpredictions file, summed over the chosen models, in the form
    taken by GNSSNetwork.clean(corrections=...): (t, {'U': ..., 'N': ..., 'E': ...}) [decimal year, m].
    The default removes the atmospheric and non-tidal ocean loading, leaving continental water in the series.
    """
    rows, cols = [], dict(LOAD_MODELS)
    with open(path) as f:
        for line in f:
            p = line.split()
            if p and p[0].lower() == 'site':            # the header names the columns: __ee_ntal __nn_ntal __uu_ntal ...
                cols.update({h[5:9].strip('_').upper(): i for i, h in enumerate(p) if h.startswith('__ee_')})
                continue
            if p and len(p) >= max(cols[m] for m in models) + 3:
                rows.append(p)
    if not rows:
        raise ValueError('no load predictions in %s' % path)
    t = np.array([float(r[2]) for r in rows])
    out = {'E': 0.0, 'N': 0.0, 'U': 0.0}
    for m in models:
        k = cols[m]
        a = np.array([[float(v) for v in r[k:k + 3]] for r in rows])
        out['E'], out['N'], out['U'] = out['E'] + a[:, 0], out['N'] + a[:, 1], out['U'] + a[:, 2]
    order = np.argsort(t)
    return t[order], {c: v[order] for c, v in out.items()}


def read_cmcf(path):
    """CM -> CF translation of one station (NGL cmcf file): (t, {'U', 'N', 'E'}) [decimal year, m], to be subtracted"""
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 6 and p[3].lower() != 'cmcf':
                try:
                    rows.append([float(v) for v in p[2:6]])
                except ValueError:
                    continue
    if not rows:
        raise ValueError('no CM/CF rows in %s' % path)
    a = np.array(rows)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], {'E': a[:, 1], 'N': a[:, 2], 'U': a[:, 3]}


def load_corrections(folder, sites=None, models=('NTAL', 'NTOL'), cmcf_folder=None):
    """
    {site: (t, {'U', 'N', 'E'})} to subtract from the series (GNSSNetwork.clean(corrections=...)): the load predictions
    of the chosen models (SITE.tenv3 in folder; none if models is empty) plus, with cmcf_folder, the CM -> CF
    translation (lowercase site files). A station is included when every requested part is available for it.
    """
    folder = Path(folder)
    if sites is None:
        sites = sorted(f.stem for f in folder.glob('*.tenv3')) if models else \
            sorted(f.name.upper() for f in Path(cmcf_folder).iterdir() if '.' not in f.name)
    out = {}
    for s in sites:
        parts = []
        if models:
            f = folder / ('%s.tenv3' % s)
            if not f.exists():
                continue
            parts.append(read_load_predictions(f, models))
        if cmcf_folder is not None:
            f = Path(cmcf_folder) / s.lower()
            if not f.exists():
                continue
            parts.append(read_cmcf(f))
        if not parts:
            continue
        t = parts[0][0]
        d = {c: sum(np.interp(t, tp, dp[c]) for tp, dp in parts) for c in ('U', 'N', 'E')}
        out[s] = (t, d)
    return out


def _is_cmcf(data):
    head = data[:300].decode('ascii', 'replace').split()
    return len(head) > 4 and head[3].lower() == 'cmcf'


def download_cmcf(sites, out_dir, overwrite=False, pause=0.2):
    """CM -> CF translation files of the stations, as out_dir/<site lowercase>"""
    out_dir = Path(out_dir)
    status = {}
    for site in sites:
        path = out_dir / site.lower()
        if path.exists() and not overwrite:
            status[site] = 'present'
            continue
        try:
            _fetch(CMCF_URL.format(site_lower=site.lower()), path, check=_is_cmcf)
            status[site] = 'ok'
        except Exception as e:
            status[site] = str(e)
        time.sleep(pause)
    return status


def download_load_predictions(sites, out_dir, overwrite=False, pause=0.2):
    """tenv3 files with load predictions (atmosphere, ocean, hydrology, ...) of the stations, as out_dir/SITE.tenv3"""
    out_dir = Path(out_dir)
    status = {}
    for site in sites:
        path = out_dir / ('%s.tenv3' % site)
        if path.exists() and not overwrite:
            status[site] = 'present'
            continue
        try:
            _fetch(LOAD_URL.format(site=site), path, check=_is_tenv3)
            status[site] = 'ok'
        except Exception as e:
            status[site] = str(e)
        time.sleep(pause)
    return status


def download_catalogue(out_dir):
    """DataHoldings.txt and steps.txt into out_dir"""
    out_dir = Path(out_dir)
    for url in (HOLDINGS_URL, STEPS_URL):
        _fetch(url, out_dir / url.rsplit('/', 1)[1])
    return out_dir


def download_series(sites, out_dir, frame='IGS20', overwrite=False, pause=0.2):
    """
    Daily .tenv3 files of the stations, saved as out_dir/SITE.tenv3.
    :return: dict {site: 'ok' | 'present' | error message}
    """
    out_dir = Path(out_dir)
    status = {}
    for site in sites:
        path = out_dir / ('%s.tenv3' % site)
        if path.exists() and not overwrite:
            status[site] = 'present'
            continue
        err = None
        for pattern in TENV3_URLS[frame]:
            try:
                _fetch(pattern.format(site=site), path, check=_is_tenv3)
                status[site], err = 'ok', None
                break
            except Exception as e:  # try the next URL pattern
                err = str(e)
        if err is not None:
            status[site] = err
        time.sleep(pause)
    return status

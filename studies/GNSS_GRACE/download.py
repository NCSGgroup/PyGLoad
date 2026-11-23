"""
Download the NGL station list, the steps file and the daily .tenv3 series of the stations in a region.
Run on a computer with internet access:

    python studies/GNSS_GRACE/download.py

Files go to <data root>/input/gnss/NGL/ (DataHoldings.txt, steps.txt, tenv3/SITE.tenv3); files already present are
skipped, so the script can be re-run to complete a download. The region and period are those of config.py.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repository root, so that `src` is importable
from src.config import FileTool
from src.preprocessing import ngl

from config import CONFIG

if __name__ == '__main__':
    out = FileTool.get_input_dir('gnss/NGL')
    if not (out / 'DataHoldings.txt').exists() or CONFIG.get('refresh_catalogue', False):
        ngl.download_catalogue(out)
    stations = ngl.read_holdings(out / 'DataHoldings.txt')
    sites = ngl.select_stations(stations, box=CONFIG['box'], start=CONFIG['start'], end=CONFIG['end'],
                                min_years=CONFIG['min_years'], min_solutions=int(CONFIG['min_years'] * 365))
    print('%d stations in %s with at least %g years between %s and %s' % (len(sites), CONFIG['box'],
                                                                         CONFIG['min_years'], CONFIG['start'],
                                                                         CONFIG['end']))
    ref = CONFIG.get('cme_reference')
    if isinstance(ref, dict):                            # far reference stations for the common mode
        ring = ngl.select_reference_ring(stations, ref['center'], ref['min_km'], ref['max_km'], CONFIG['start'],
                                         CONFIG['end'], CONFIG['min_years'])
        print('%d reference stations %g-%g km from %s' % (len(ring), ref['min_km'], ref['max_km'], ref['center']))
        sites = sorted(set(sites) | set(ring))
    status = ngl.download_series(sites, out / 'tenv3', frame=CONFIG['frame'])
    failed = {k: v for k, v in status.items() if v not in ('ok', 'present')}
    print('downloaded %d, already present %d, failed %d' % (sum(v == 'ok' for v in status.values()),
                                                           sum(v == 'present' for v in status.values()), len(failed)))
    for k, v in list(failed.items())[:10]:
        print('  %s: %s' % (k, v))
    if CONFIG.get('load_models'):
        status = ngl.download_load_predictions(sites, out / 'loadpred')
        failed = {k: v for k, v in status.items() if v not in ('ok', 'present')}
        print('load predictions: downloaded %d, already present %d, failed %d'
              % (sum(v == 'ok' for v in status.values()), sum(v == 'present' for v in status.values()), len(failed)))
        for k, v in list(failed.items())[:10]:
            print('  %s: %s' % (k, v))
    if CONFIG.get('cm_to_cf') == 'ngl':
        status = ngl.download_cmcf(sites, out / 'loadpred_cmcf')
        failed = {k: v for k, v in status.items() if v not in ('ok', 'present')}
        print('CM/CF translations: downloaded %d, already present %d, failed %d'
              % (sum(v == 'ok' for v in status.values()), sum(v == 'present' for v in status.values()), len(failed)))
        for k, v in list(failed.items())[:10]:
            print('  %s: %s' % (k, v))

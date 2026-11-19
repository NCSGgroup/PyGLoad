"""
Compare the scores of run.py with the reference scores of the study (<data root>/input/GNSS_TWS_validation/results/).
Writes <data root>/output/GNSS_TWS_validation/verify<suffix>.txt and prints it; exit code 1 if a median differs.

    python studies/GNSS_TWS_validation/verify.py               # run.py with the reference GNSS series
    python studies/GNSS_TWS_validation/verify.py _rebuild      # run.py with GNSS_TWS_GNSS=rebuild

Tolerances: medians and the n of every (group, product, score) must agree to rounding (|diff| <= 0.005, phase
<= 0.5 day, n exact). The station-bootstrap 5-95 % bounds may move (typically 0.001, up to 0.04 in small regions and
2 days for the phase), because the random draws depend on the list of products scored together (the reference run had
more products); they are reported, not tested.
"""
import sys

import numpy as np
import pandas as pd

from config import reference_dir, output_dir

TOL = dict(default=0.005, phase_diff_d=0.5)

if __name__ == '__main__':
    suffix = sys.argv[1] if len(sys.argv) > 1 else ''
    out = output_dir()
    ref = pd.read_csv(reference_dir('results/gnss_summary_public.csv'))
    new = pd.read_csv(out / ('gnss_summary%s.csv' % suffix))
    key = ['group', 'product', 'score']
    M = ref.merge(new, on=key, suffixes=('_ref', '_new'), how='inner')
    M = M[M['product'].isin(new['product'].unique())]
    lines = ['GNSS-TWS validation: run.py%s vs reference scores' % (' (%s)' % suffix if suffix else ''),
             'compared (group, product, score) rows: %d of %d reference rows for the products run'
             % (len(M), int(ref['product'].isin(new['product'].unique()).sum())), '']
    tol = np.where(M.score == 'phase_diff_d', TOL['phase_diff_d'], TOL['default'])
    dmed = (M.median_new - M.median_ref).abs()
    bad = M[(dmed > tol) | (M.n_new != M.n_ref)]
    dlo, dhi = (M.p05_new - M.p05_ref).abs(), (M.p95_new - M.p95_ref).abs()
    lines.append('%-8s %6s %14s %14s %12s' % ('product', 'rows', 'max|d median|', 'n mismatches', 'max|d p05/95|'))
    for p, g in M.groupby('product'):
        i = g.index
        lines.append('%-8s %6d %14.2e %14d %12.3f' % (p, len(g), dmed[i].max(), int((g.n_new != g.n_ref).sum()),
                                                      max(dlo[i].max(), dhi[i].max())))
    lines.append('')
    lines.append('Main table (global, medians): product, reference -> this run')
    for p in sorted(M['product'].unique()):
        x = M[(M.group == 'global') & (M['product'] == p)].set_index('score')
        lines.append('  %-8s n %d -> %d | WRMS red %% %.1f -> %.1f | r_total %.3f -> %.3f | amp err factor %.3f -> %.3f'
                     % (p, x.loc['wrms_red', 'n_ref'], x.loc['wrms_red', 'n_new'],
                        100 * x.loc['wrms_red', 'median_ref'], 100 * x.loc['wrms_red', 'median_new'],
                        x.loc['r_total', 'median_ref'], x.loc['r_total', 'median_new'],
                        np.exp(x.loc['amp_ratio', 'median_abs_log_ref']), np.exp(x.loc['amp_ratio', 'median_abs_log_new'])))
    # per station
    rs = pd.read_csv(reference_dir('results/gnss_scores_per_station_public.csv'), index_col=0)
    ns = pd.read_csv(out / ('gnss_scores_per_station%s.csv' % suffix), index_col=0)
    cols = [c for c in ns.columns if '|' in c and c in rs.columns]
    common = rs.index.intersection(ns.index)
    d = (ns.loc[common, cols] - rs.loc[common, cols]).abs()
    lines.append('')
    lines.append('per station: %d stations in both (reference %d, this run %d); max |diff| over %d score columns %.2e '
                 '(values stored with 4 significant digits)' % (len(common), len(rs), len(ns), len(cols),
                                                                np.nanmax(d.values)))
    lines.append('')
    lines.append('RESULT: %s' % ('PASS (all medians and n agree to rounding)' if bad.empty else
                                 'DIFFERENT in %d rows, e.g.\n%s' % (len(bad), bad[key + ['n_ref', 'n_new', 'median_ref',
                                                                                         'median_new']].head(20).to_string(index=False))))
    txt = '\n'.join(lines) + '\n'
    (out / ('verify%s.txt' % suffix)).write_text(txt)
    print(txt)
    sys.exit(0 if bad.empty else 1)

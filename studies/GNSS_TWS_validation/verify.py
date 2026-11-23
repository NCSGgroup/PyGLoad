"""
Check a run of run.py against the published numbers. Prints PASS / FAIL per check; exit code 1 if any check fails.

    python studies/GNSS_TWS_validation/verify.py                         # after run.py (CSR + JPL, both frames)
    python studies/GNSS_TWS_validation/verify.py --reference DIR         # also per station against reference results

Checks:
  1. global median WRMS reduction of CSR and JPL, rounded to two decimals, equals the published value
     (IGS20-CF: CSR 0.15, JPL 0.13; IGS14: CSR 0.06, JPL 0.06);
  2. with --reference (default: the folder reference_results/ next to the data root, if it exists): every per-station
     score (r_total, r_seas, r_nonseas, r_nonseas_dt, amp_ratio, phase_diff_d, wrms_red, wrms_red_nodt) of every
     product present in both agrees within 1e-6 (phase within 1e-6 day), and the medians and n of the summaries agree.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import config
from config import FRAMES

PUBLISHED = {('igs20cf', 'CSR'): 0.15, ('igs20cf', 'JPL'): 0.13, ('igs14', 'CSR'): 0.06, ('igs14', 'JPL'): 0.06}
SCORES = ['r_total', 'r_seas', 'r_nonseas', 'r_nonseas_dt', 'amp_ratio', 'phase_diff_d', 'wrms_red', 'wrms_red_nodt']
TOL = 1e-6


def main():
    ap = argparse.ArgumentParser(description='check run.py results against the published numbers')
    ap.add_argument('--data-root', help='PyGLoad data root (default: installation/pygload_paths.cfg or PYGLOAD_DATA)')
    ap.add_argument('--reference', help='folder with <frame>/gnss_scores_per_station.csv and gnss_summary.csv')
    ap.add_argument('--frames', default='igs20cf,igs14')
    a = ap.parse_args()
    if a.data_root:
        config.set_data_root(a.data_root)
    ref = Path(a.reference) if a.reference else config.FileTool.get_data_dir().parent / 'reference_results'
    lines, ok = [], True

    def check(cond, text):
        nonlocal ok
        ok &= bool(cond)
        lines.append('%s  %s' % ('PASS' if cond else 'FAIL', text))

    for frame in [f for f in a.frames.split(',') if f]:
        out = config.output_dir(frame)
        f = out / 'gnss_summary.csv'
        if not f.exists():
            check(False, '%s: %s missing (run run.py first)' % (frame, f))
            continue
        SM = pd.read_csv(f)
        for p in ('CSR', 'JPL'):
            x = SM[(SM.group == 'global') & (SM['product'] == p) & (SM.score == 'wrms_red')]
            if x.empty:
                lines.append('----  %s %s not run' % (frame, p))
                continue
            v = float(x['median'].iloc[0])
            check(round(v, 2) == PUBLISHED[(frame, p)], '%s %s global median WRMS reduction %.4f (n %d), published %.2f'
                  % (frame, p, v, int(x.n.iloc[0]), PUBLISHED[(frame, p)]))
        if not ref.is_dir():
            continue
        rs = ref / frame / 'gnss_scores_per_station.csv'
        if rs.exists():
            R = pd.read_csv(rs, index_col=0)
            N = pd.read_csv(out / 'gnss_scores_per_station.csv', index_col=0)
            prods = sorted({c.split('|')[0] for c in N.columns if '|' in c} & {c.split('|')[0] for c in R.columns if '|' in c})
            common = R.index.intersection(N.index)
            check(len(common) == len(R) == len(N), '%s stations: %d in both (reference %d, this run %d)'
                  % (frame, len(common), len(R), len(N)))
            for p in prods:
                cols = ['%s|%s' % (p, s) for s in SCORES]
                a_, b_ = N.loc[common, cols].values, R.loc[common, cols].values
                same_nan = np.array_equal(np.isnan(a_), np.isnan(b_))
                d = np.nanmax(np.abs(a_ - b_))
                check(same_nan and d <= TOL, '%s %s per station: max |diff| over %d scores x %d stations %.2e (tol %.0e)'
                      % (frame, p, len(SCORES), len(common), d, TOL))
        rm = ref / frame / 'gnss_summary.csv'
        if rm.exists():
            R = pd.read_csv(rm)
            M = R.merge(SM, on=['group', 'product', 'score'], suffixes=('_ref', '_new'))
            if len(M):
                dm = (M.median_new - M.median_ref).abs().max()
                check(dm <= 1e-4 and (M.n_new == M.n_ref).all(),
                      '%s summary: %d rows (products %s), max |median diff| %.1e, n equal %s'
                      % (frame, len(M), ','.join(sorted(M['product'].unique())), dm, (M.n_new == M.n_ref).all()))
    lines.append('')
    lines.append('RESULT: %s' % ('PASS' if ok else 'FAIL'))
    txt = '\n'.join(lines) + '\n'
    (config.output_dir() / 'verify.txt').write_text(txt)
    print(txt)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()

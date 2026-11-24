import pandas as pd, numpy as np

df = pd.read_csv('outputs/control/run_20261119_205739/control_results.csv')
print('=== RESULTS SUMMARY ===')
print(df[['scenario','source','status','J_true','damage_ratio','effort_norm','peak_ratio','blocked']].to_string(index=False))
print()
rb = pd.read_csv('outputs/control/run_20261119_205739/timing_robustness.csv')
print('=== TIMING ROBUSTNESS ===')
print(rb.round(4).to_string(index=False))
print()
ok = df[df.status == 'ok']
print('Scenarios OK: %d, blocked: %d, not blocked: %d' % (len(ok), ok.blocked.sum(), (~ok.blocked).sum()))
oracle_J = float(df.loc[df.scenario == 'oracle', 'J_true'].iloc[0])
print('Oracle J_true: %.4f' % oracle_J)
print('Mean J_true (ok, non-oracle): %.4f' % float(ok[ok.source != 'oracle'].J_true.mean()))

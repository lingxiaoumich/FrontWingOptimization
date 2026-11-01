"""Compare test_solver_settings.py results against the reference setting.

  python scripts/compare_solver_settings.py optimizations/solver_tolerance_test_20261031
Reads runs/<case>_<setting>/step_test.json and writes comparison.json and
comparison.txt beside runs/. Gradient error is |g - g_ref| / |g_ref| over the
full design vector (incidence + 20 shape variables), per function.
"""
import json
import sys
from pathlib import Path

import numpy as np

T = Path(sys.argv[1])
res = {}
for f in sorted(T.glob('runs/*/step_test.json')):
    case, cfg = f.parent.name.split('_', 1)
    res.setdefault(case, {})[cfg] = json.loads(f.read_text())
def vec(g): return np.concatenate([np.ravel(g[k]) for k in sorted(g)])
rows, lines = [], []
for case in sorted(res):
    ref = res[case].get('ref')
    if ref is None: continue
    lines.append(f'case {case}')
    lines.append(f"  {'setting':15s} {'flow s':>7s} {'adj s':>6s} {'total':>6s} {'speedup':>7s} {'dCL rel':>9s} {'dCD rel':>9s} "
                 f"{'grad CL err':>11s} {'grad CD err':>11s}  status")
    for cfg, r in sorted(res[case].items(), key=lambda kv: (kv[0] != 'ref', kv[0])):
        row = {'case': case, 'setting': cfg, 'step_failed': r['step_failed'], 'flow_time_s': r['flow_time_s'],
               'adjoint_time_s': r.get('adjoint_time_s'), 'settings': r['settings']}
        status = 'FLOW FAILED' if r['step_failed'] else ('ADJOINT FAILED' if r.get('adjoint_failed') else 'ok')
        if not r['step_failed'] and 'gradient' in r:
            tot = r['flow_time_s']+r['adjoint_time_s']; reft = ref['flow_time_s']+ref['adjoint_time_s']
            row.update({'total_s': tot, 'speedup': reft/tot, 'cl_rel_diff': (r['cl']-ref['cl'])/abs(ref['cl']),
                        'cd_rel_diff': (r['cd']-ref['cd'])/abs(ref['cd'])})
            for fn in ['objective', 'drag']:
                g, g0 = vec(r['gradient'][fn]), vec(ref['gradient'][fn])
                row[f'grad_{fn}_rel_err'] = float(np.linalg.norm(g-g0)/np.linalg.norm(g0))
            lines.append(f"  {cfg:15s} {r['flow_time_s']:7.0f} {r['adjoint_time_s']:6.0f} {tot:6.0f} {row['speedup']:6.2f}x "
                         f"{row['cl_rel_diff']:9.1e} {row['cd_rel_diff']:9.1e} {row['grad_objective_rel_err']:11.1e} "
                         f"{row['grad_drag_rel_err']:11.1e}  {status}")
        else:
            lines.append(f"  {cfg:15s} {r['flow_time_s']:7.0f} {'':>6s} {'':>6s} {'':>7s} {'':>9s} {'':>9s} {'':>11s} {'':>11s}  {status}")
        rows.append(row)
(T/'comparison.json').write_text(json.dumps(rows, indent=1)+'\n')
(T/'comparison.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))

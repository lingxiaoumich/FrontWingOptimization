"""One line per run of kkt_queue.py: its current round and where PSQP is (python3 status.py)."""
import json
import time
from pathlib import Path

O = Path('/home/lingxiao77/AER588/aero588_project/optimizations')
RUNS = ['run_036_area050MinLossClc10CpMin4CapStart', 'run_037_area050MinLossClc08CpMin4CapStart',
        'run_038_area050MinLossClc14CpMin4CapStart', 'run_039_area050MinLossClc12CpMin4']
LIMITS = {'run_036': -1.0, 'run_037': -0.8, 'run_038': -1.4, 'run_039': -1.2}
for name in RUNS:
    link = O/name
    rounds = sorted((r for r in link.glob('r[0-9]*') if (r/'command.txt').exists() and 'KKT continuation' in (r/'command.txt').read_text()),
                    key=lambda r: r.stat().st_mtime)
    if not rounds: print(name[:7], 'no KKT round yet'); continue
    r = rounds[-1]
    rows = [json.loads(l) for l in (r/'history.jsonl').read_text().splitlines() if l.strip()] if (r/'history.jsonl').exists() else []
    its = sorted(int(d.name[5:]) for d in (r/'iterates').glob('eval_*') if (d/'design.json').exists())
    ok = [x for x in rows if not x['cfd_failed']]
    done = (r/'run_complete.json').exists()
    line = f"{name[:7]} {r.name}: {len(rows)} evals ({sum(x['cfd_failed'] for x in rows)} failed), {len(its)} iterates"
    if its:
        q = json.loads((r/'iterates'/f'eval_{its[-1]:03d}'/'result.json').read_text())
        line += (f"; last iterate eval {its[-1]}: loss {q['cpt_loss_outlet']:.6f}, C_L c {.5*q['cl']:.4f} (<= {LIMITS[name[:7]]}), "
                 f"Cp_min {q['cp_min_ks']:.3f} (>= -4), feasible {q['feasible']}")
    if ok: line += f"; {(time.time()-(r/'history.jsonl').stat().st_mtime)/60:.0f} min since the last evaluation"
    if done: line += f"; PSQP exit {json.loads((r/'run_complete.json').read_text()).get('inform')}"
    print(line)

"""Refresh the post-processing plots of a running optimization as its PSQP iterates are saved.

Give the chain of run folders in order; the last one is the run being watched (it needs
--save-iterate-flows). Every few new iterates this re-runs
  plot_iterate_flows.py --cache   on the watched run: <run>/plots/cpt_iterates_*, cp_iterates_*
  plot_run_chain.py               on the whole chain:  <out>/plots/chain_history, chain_sections
and once more when the run has stopped, then exits. A failed refresh is logged and left to
the next one. Each refresh takes some 10 s of a core from the solver; the defaults keep
that near 5 % of the run. Run with the MACH-Aero env, in the background:
  python scripts/live_plots.py --out optimizations/<chain>_summary_<date> RUN1 RUN2 ... > live_plots.log 2>&1 &
"""
import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

S = Path(__file__).resolve().parent
p = argparse.ArgumentParser()
p.add_argument('runs', type=Path, nargs='+')
p.add_argument('--out', type=Path, required=True, help='summary folder for the chain plots')
p.add_argument('--chain-every', type=int, default=2, help='new iterates between refreshes of the chain plots')
p.add_argument('--flows-every', type=int, default=5, help='new iterates between refreshes of the flow renders')
p.add_argument('--poll', type=float, default=20., help='s between looks at the run folder')
a = p.parse_args()
runs = [r.resolve() for r in a.runs]; run = runs[-1]
name = re.sub(r'_gap\w+_\d{8}$', '', run.name)  # as the earlier runs' plot names

# numpy's BLAS would otherwise start a thread per core, on cores the solver is using
ONE_THREAD = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
def saved():
    """Iterates written so far (the optimizer writes design.json last)."""
    return sum((d/'design.json').exists() for d in (run/'iterates').glob('eval_*'))
def running():
    return subprocess.run(['pgrep', '-f', f'[p]ython -u {run}/source'], capture_output=True).returncode == 0
def refresh(script, *args):
    done = subprocess.run([sys.executable, str(S/script), *map(str, args)], capture_output=True, text=True, env=ONE_THREAD)
    why = '' if done.returncode == 0 else ': FAILED, ' + (done.stderr.strip().splitlines() or ['no message'])[-1]
    print(f"{time.strftime('%H:%M:%S')} {script} at {saved()} iterates{why}", flush=True)

flows = chain = 0  # iterates at the last refresh
while True:
    alive, n = running(), saved()
    if n and (n-flows >= a.flows_every or not alive):
        refresh('plot_iterate_flows.py', '--cache', '--out', run, '--name', name, run); flows = n
    if n and (n-chain >= a.chain_every or not alive):  # after the flows: it reads their recirculation areas
        refresh('plot_run_chain.py', '--out', a.out, *runs); chain = n
    if not alive: break
    time.sleep(a.poll)
print(f"{time.strftime('%H:%M:%S')} the run has stopped; plots are final at {n} iterates", flush=True)

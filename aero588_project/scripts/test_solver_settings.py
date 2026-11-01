"""Time and compare one optimizer step under different flow/adjoint solver settings.

Runs exactly what the optimizer runs (the setup, objfun and sensfun of
run_slip_optimization.py): solve the design --restart-design from its converged
flow --restart-flow, then solve the step design --step-to and evaluate both
adjoints there. Writes CL, CD, residuals, timings and the full gradients to
<root>/step_test.json. All other arguments go to run_slip_optimization.py, so
--solver/--nk-switch-tol/--l2-convergence/--adjoint-l2-convergence select the
settings under test.
  mpirun -np 32 python scripts/test_solver_settings.py --step-to B.json --root R \
      --ffd-stations 5 --shape-directions xy --drag-cap 0.1 --min-local-thickness 0.002 \
      --solver ank_deep --ncycles 20000 --reuse-checks --restart-design A.json --restart-flow A_vol.cgns
"""
import argparse
import json
import sys
import time
from pathlib import Path

h = argparse.ArgumentParser(add_help=False)
h.add_argument('--step-to', type=Path, required=True)
own, rest = h.parse_known_args()
src = Path(__file__).with_name('run_slip_optimization.py')
sys.argv = [str(src)]+rest
sys.path.insert(0, str(src.parent))
code = src.read_text().split('# Check one independent CFD derivative')[0]
exec(compile(code, str(src), 'exec'))

x0 = load_design(a.restart_design); x1 = load_design(own.step_to)
f0, fail0 = objfun(x0)
tic = time.perf_counter(); f1, fail1 = objfun(x1); t_flow = time.perf_counter()-tic
rec = cache['record']
out = {'settings': {k: options[k] for k in ['NKSwitchTol', 'L2Convergence', 'L2ConvergenceRel', 'adjointL2Convergence']},
       'solver': a.solver, 'ranks': comm.size, 'start_failed': fail0, 'step_failed': fail1,
       'flow_time_s': t_flow, 'cl': rec['cl'], 'cd': rec['cd'], 'residual_ratio_initial': rec['residual_ratio_initial'],
       'residual_start': rec['residual_start'], 'residual_final': rec['residual_final']}
if not fail1:
    tic = time.perf_counter(); g, gfail = sensfun(x1, f1); t_adj = time.perf_counter()-tic
    out.update({'adjoint_time_s': t_adj, 'adjoint_failed': gfail,
                'gradient': {f: realdict(g[f]) for f in ['objective', 'drag']}})
dump(root/'step_test.json', out)
if comm.rank == 0: print('STEP_TEST', json.dumps({k: v for k, v in out.items() if k != 'gradient'}), flush=True)

"""Run the optimizer's aerodynamic derivative check at any design, with a chosen incidence step.

Runs aerodynamic_derivative_check of run_slip_optimization.py at --at, starting from
--restart-flow if given, with incidence step --step-deg (the optimizer uses 0.01 deg), and
writes <root>/derivative_check_at.json. Other arguments go to run_slip_optimization.py.
  mpirun -np 16 python scripts/check_derivatives_at.py --at D.json --step-deg 0.05 --root R \
      --ffd-stations 5 --shape-directions xy --drag-cap 0.4 --nk-switch-tol 1e-6 ...
"""
import argparse
import json
import sys
from pathlib import Path

h_ = argparse.ArgumentParser(add_help=False)
h_.add_argument('--at', type=Path, required=True)
h_.add_argument('--step-deg', type=float, default=1e-2)
own, rest = h_.parse_known_args()
src = Path(__file__).with_name('run_slip_optimization.py')
sys.argv = [str(src)]+rest
sys.path.insert(0, str(src.parent))
code = src.read_text().split('# Baseline always starts from the same NACA shape')[0]  # the setup and functions
exec(compile(code.replace('h=1e-2', f'h={own.step_deg!r}'), str(src), 'exec'))

x = load_design(own.at)
try:
    aerodynamic_derivative_check(x)
    passed = True
except RuntimeError:
    passed = False
report = json.loads((root/'source/aerodynamic_derivative_check.json').read_text())
report.update({'design': str(own.at.resolve()), 'step_deg': own.step_deg, 'passed_5pct': passed and report['verdict'] == 'passed',
               'cl_base': cache['record'].get('cl')})
dump(root/'derivative_check_at.json', report)
if comm.rank == 0: print('DERIV_CHECK_AT', json.dumps(report), flush=True)

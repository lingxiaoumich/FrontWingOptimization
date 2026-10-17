"""Smoke-check the MACH-Aero environment; does not run CFD."""
import importlib
import json
import platform
from pathlib import Path

from mpi4py import MPI

names = ['baseclasses', 'pyspline', 'pygeo', 'idwarp', 'adflow', 'pyhyp', 'pyoptsparse', 'multipoint', 'prefoil']
modules = {}
for name in names:
    module = importlib.import_module(name)
    modules[name] = {'version': str(getattr(module, '__version__', 'unknown')), 'path': module.__file__}
from pyoptsparse import OPT
OPT('SLSQP')
comm = MPI.COMM_WORLD
assert comm.allreduce(1) == comm.size
if comm.rank == 0:
    report = {'python': platform.python_version(), 'modules': modules, 'optimizer': 'SLSQP instantiated', 'mpi_ranks': comm.size, 'status': 'imports and MPI passed', 'cfd_run': False}
    path = Path('/project/reports/environment_check.json')
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))

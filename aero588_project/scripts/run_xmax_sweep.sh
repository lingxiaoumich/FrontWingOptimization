#!/bin/bash
# PSQP step-limit (XMAX) sweep: short restarts from the October 29 final design.
# Run inside the MDO Lab container with aero588_project mounted at /project:
#   CHECKS=<passed aerodynamic_derivative_check.json> [SOLVER=main|ank_uncoupled|ank_deep] \
#   [SWEEP_NAME=xmax_sweep_20261031] [NCYCLES=20000] [TIMEOUT=150m] \
#   bash /project/scripts/run_xmax_sweep.sh [maxit] [xmax ...]
# Each XMAX gets a fresh run folder built from the prepared five-pair x/y setup and
# starts from the converged October 29 final flow (--restart-flow), so no run has to
# solve the NACA baseline and then jump to the restart design: from the baseline flow
# that jump stalled at a 9e-7 residual ratio in 4000 cycles and was still running
# after 30 min with 20000. NCYCLES counts NK linear iterations too. Fortran output is
# unbuffered so the residual trace survives a kill. A run that exceeds TIMEOUT is
# killed; its history.jsonl keeps the evaluations completed up to that point.
set -uo pipefail
source /home/mdolabuser/.bashrc_mdolab
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

MAXIT=${1:-10}; shift || true
XMAXES=("${@:-0.05 0.2 0.5 1.0}"); read -r -a XMAXES <<< "${XMAXES[*]}"
P=/project
SWEEP=$P/optimizations/${SWEEP_NAME:-xmax_sweep_20261031}
RUN=$P/optimizations/slip_ground_xy5_cd005_20261029
CHECKS=${CHECKS:?set CHECKS to a passed aerodynamic_derivative_check.json}
SOLVER=${SOLVER:-main}
NCYCLES=${NCYCLES:-20000}
TIMEOUT=${TIMEOUT:-150m}
mkdir -p "$SWEEP"; cp "$RUN/final/design.json" "$SWEEP/start_design.json"
echo "solver=$SOLVER ncycles=$NCYCLES timeout=$TIMEOUT maxit=$MAXIT xmax=${XMAXES[*]}" > "$SWEEP/sweep_settings.txt"

for x in "${XMAXES[@]}"; do
  root=$SWEEP/xmax_$x
  if [[ -f $root/run_complete.json ]]; then echo "skip $root (complete)"; continue; fi
  mkdir -p "$root/source"
  cp "$P"/setups/slip_ground_xy5_cd005/source/* "$root/source/"
  cp "$P/scripts/run_slip_optimization.py" "$CHECKS" "$root/source/"
  echo "$(date '+%F %T') start XMAX=$x solver=$SOLVER"
  timeout "$TIMEOUT" mpirun -np 4 -x GFORTRAN_UNBUFFERED_PRECONNECTED=y python -u \
      "$root/source/run_slip_optimization.py" --root "$root" \
      --ffd-stations 5 --shape-directions xy --drag-cap 0.05 --optimizer PSQP \
      --maxit "$MAXIT" --xmax "$x" --ncycles "$NCYCLES" --solver "$SOLVER" --reuse-checks \
      --restart-design "$SWEEP/start_design.json" --restart-flow "$RUN/final/flow_000_vol.cgns" \
      > "$root/run.log" 2>&1
  rc=$?
  echo "$(date '+%F %T') end XMAX=$x exit=$rc (124 = timed out)"
done

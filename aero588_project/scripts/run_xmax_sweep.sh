#!/bin/bash
# PSQP step-limit (XMAX) sweep: short restarts from the October 29 final design.
# Run inside the MDO Lab container with aero588_project mounted at /project:
#   CHECKS=<passed aerodynamic_derivative_check.json> [NCYCLES=20000] \
#   bash /project/scripts/run_xmax_sweep.sh [maxit] [xmax ...]
# Each XMAX gets a fresh run folder built from the prepared five-pair x/y setup and
# starts from the converged October 29 final flow (--restart-flow), so no run has to
# solve the NACA baseline and then jump to the restart design: from the baseline flow
# that jump stalled at a 9e-7 residual ratio in 4000 cycles and was still running
# after 30 min with 20000. NCYCLES counts NK linear iterations too.
set -uo pipefail
source /home/mdolabuser/.bashrc_mdolab
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

MAXIT=${1:-10}; shift || true
XMAXES=("${@:-0.05 0.2 0.5 1.0}"); read -r -a XMAXES <<< "${XMAXES[*]}"
P=/project
SWEEP=$P/optimizations/xmax_sweep_20261031
RUN=$P/optimizations/slip_ground_xy5_cd005_20261029
CHECKS=${CHECKS:?set CHECKS to a passed aerodynamic_derivative_check.json}
NCYCLES=${NCYCLES:-20000}
mkdir -p "$SWEEP"; cp "$RUN/final/design.json" "$SWEEP/start_design.json"

for x in "${XMAXES[@]}"; do
  root=$SWEEP/xmax_$x
  if [[ -f $root/run_complete.json ]]; then echo "skip $root (complete)"; continue; fi
  mkdir -p "$root/source"
  cp "$P"/setups/slip_ground_xy5_cd005/source/* "$root/source/"
  cp "$P/scripts/run_slip_optimization.py" "$CHECKS" "$root/source/"
  echo "$(date '+%F %T') start XMAX=$x"
  mpirun -np 4 python "$root/source/run_slip_optimization.py" --root "$root" \
      --ffd-stations 5 --shape-directions xy --drag-cap 0.05 --optimizer PSQP \
      --maxit "$MAXIT" --xmax "$x" --ncycles "$NCYCLES" --reuse-checks \
      --restart-design "$SWEEP/start_design.json" --restart-flow "$RUN/final/flow_000_vol.cgns" \
      > "$root/run.log" 2>&1
  rc=$?
  echo "$(date '+%F %T') end XMAX=$x exit=$rc"
done

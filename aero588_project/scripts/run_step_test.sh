#!/bin/bash
# Re-solve one saved optimizer trial design from the converged October 29 final flow,
# with a chosen solver variant and live (unbuffered) residual output.
# Run inside the container, aero588_project mounted at /project:
#   bash /project/scripts/run_step_test.sh <trial design.json> <solver> <run name> [minutes]
set -uo pipefail
source /home/mdolabuser/.bashrc_mdolab
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
P=/project
SWEEP=$P/optimizations/xmax_sweep_20261031
RUN=$P/optimizations/slip_ground_xy5_cd005_20261029
CHECKS=$SWEEP/failed_ncycles4000_xmax_0.05/source/aerodynamic_derivative_check.json
design=$1; solver=$2; root=$SWEEP/$3; minutes=${4:-45}
mkdir -p "$root/source"
cp "$P"/setups/slip_ground_xy5_cd005/source/* "$root/source/"
cp "$P/scripts/run_slip_optimization.py" "$CHECKS" "$root/source/"
cp "$design" "$root/trial_design.json"
echo "$(date '+%F %T') start $3 solver=$solver"
timeout "${minutes}m" mpirun -np 4 -x GFORTRAN_UNBUFFERED_PRECONNECTED=y python -u \
    "$root/source/run_slip_optimization.py" --root "$root" \
    --ffd-stations 5 --shape-directions xy --drag-cap 0.05 --optimizer PSQP \
    --ncycles 20000 --solver "$solver" --reuse-checks \
    --restart-design "$SWEEP/start_design.json" --restart-flow "$RUN/final/flow_000_vol.cgns" \
    --evaluate-design "$root/trial_design.json" > "$root/run.log" 2>&1
rc=$?
echo "$(date '+%F %T') end $3 exit=$rc (124 = timed out)"

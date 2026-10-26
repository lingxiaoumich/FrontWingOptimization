#!/bin/bash
# Solver (ANK) test + mesh study. Run inside the container from /project (env/fwo.sh):
#   bash study/run_campaign.sh "0 -3" "L2 L15 L1" "repo ank_uncoupled ank_off"
# Resumable: meshes and solves that already have output are skipped.
# Logs in study/logs, results in study/results/solve.
cd /project
ALPHAS=${1:-"0 -3"}
LEVELS=${2:-"L2 L15 L1"}
SOLVERS=${3:-"repo ank_uncoupled ank_off"}
NP=${NP:-8}
for a in $ALPHAS; do
  for L in $LEVELS; do
    C=$(python -c "import sys; sys.path.insert(0,'study'); from common import case_name; print(case_name('$L', $a))")
    if [ ! -f study/grids/$C/case.json ]; then
      timeout 1800 mpirun -np 4 python study/build_mesh.py --level $L --alpha $a > study/logs/build_$C.log 2>&1 \
        || { echo "$(date +%T) BUILD FAILED $C"; continue; }
      echo "$(date +%T) built $C: $(tail -1 study/logs/build_$C.log | cut -c1-80)"
    fi
    for s in $SOLVERS; do
      [ -f study/results/solve/${C}__$s.json ] && continue
      NC=20000; [ $s = ank_off ] && NC=30000   # nCycles counts total (linear) iterations
      timeout ${TMAX:-10800} mpirun -np $NP python study/solve.py --case study/grids/$C --solver $s --ncycles $NC \
        > study/logs/${C}__$s.log 2>&1 || echo "$(date +%T) SOLVE FAILED/TIMEOUT $C $s"
      echo "$(date +%T) $C $s: $(grep -E '"(solve_s|residual_reduction|cl|cd|solve_failed)"' study/logs/${C}__$s.log | tr -d ' \n')"
    done
  done
done

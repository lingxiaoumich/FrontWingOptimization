"""Queued after runs 031-034 (comparisons/minloss_20261121/minloss_queue.py): run_016_area050 out of ground effect
(user, 2026-11-21: "based on 16, same box, but not in ground effect (no ground)"):
  run_035_area050FreeAir   run_016's setup (maximize downforce; level NACA 0012 start, 20 FFD pairs moving only up and
                           down, tilt locked, the same box, section area >= 50 % of the start, no drag or loss cap),
                           with no ground: the slip wall at y = 0 is replaced by more of the domain, down to a far
                           field 10 chords below (as the top is 10 chords above).
Grid (remesh_design.py --free-air, kept by every remesh of the chain through baseline.json 'mesh'): the original near
grid and the original background lines above y = 0, then 2.5 mm cells down to y = -0.15 m, growing 8 % a cell to
y = -5 m, a far field there; 128,616 cells before blanking (88,548 with the ground). The driver runs with --no-ground:
no y = 0 symmetry plane in the mesh warp, and the outlet loss plane spans y = -5..5 m, so its mass-averaged loss is
over twice the height of the ground runs' (compare the flux form cpt_loss_outlet_flux_c_m). No initial flow fits this
grid: the first solve starts from freestream (--first-ncycles 60000), at the start tilted 1 deg (--first-solve-tilt 1),
and the level start is reached from that flow by the stalled-solve walk. The 4-rank test (scratchpad testC_freeair,
2026-11-21) found 0 orphans, but solved directly from freestream the level start's symmetric zero-lift flow stalled
between 2e-9 and 2e-7 of its initial residual for 35,000 steps (target 1e-10; lift flipping sign at 1e-5).
Since 2026-11-21 21:39 this job is run by comparisons/priority_20261121/priority_queue.py (build() only).
Lanes are shared with the other queues (lane_A/B.busy in comparisons/floor_20261121); the lane is taken only once
run_034 has started, and after READY (written once the free-air test passed).
  setsid nohup taskset -c 12-63 python -u freeair_queue.py > queue.log 2>&1 &
"""
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # chain(), finish(), claim() and the lane locks; imports lanes.py as FQ.L and run_chain as FQ.RC
L = FQ.L
READY = HERE/'READY'
AFTER = (HERE.parent/'minloss_20261121/queue.log', 'run_034_area050MinLossClc12: on lane')
NAME = 'run_035_area050FreeAir'

def build(name, lane):
    link, real = L.O/name, L.BIG/name
    if link.exists(): raise RuntimeError(f'{link} exists already')
    real.mkdir(parents=True)
    if real != link: link.symlink_to(real)
    grid = link/'grids/r00_level_freeair'; grid.mkdir(parents=True); (link/'summary').mkdir()
    shutil.copy(L.START/'grids/r00_level_start/source_design.json', grid/'source_design.json')
    (grid/'source.txt').write_text('run_006 level start (grids/r00_level_start/source_design.json), on a free-air grid\n')
    made = subprocess.run([FQ.PY, str(L.S/'remesh_design.py'), '--design', str(grid/'source_design.json'), '--out', str(grid),
                           '--hyp-options', FQ.RC.HYP[0], '--free-air'], capture_output=True, text=True, env=FQ.RC.ONE_THREAD, timeout=600)
    (grid/'remesh.log').write_text(made.stdout+made.stderr)
    if made.returncode: raise RuntimeError(f'the start grid failed; see {grid}/remesh.log')
    (grid/'offset_design.json').write_text(json.dumps({'incidence': [-3.02], 'shape': [0.]*40}))
    root = link/'r00'; src = root/'source'; src.mkdir(parents=True)
    for f in L.INPUTS:
        if f not in ('restart_vol.cgns', 'ground.cgns', 'airfoil_coordinates.csv'): shutil.copy(L.START/'r00/source'/f, src/f)
    for f in ('ground.cgns', 'airfoil_coordinates.csv'): shutil.copy(grid/f, src/f)
    for f in (L.START/'r00/source').glob('*.py'): shutil.copy(L.S/f.name, src/f.name)
    shutil.copy(FQ.RANKFILE[lane], root/'rankfile')
    cmd = (f'mpirun -np 24 --map-by rankfile:file=rankfile --bind-to core --mca btl self,sm -x GFORTRAN_UNBUFFERED_PRECONNECTED=y '
           f'python -u {src}/run_slip_optimization.py --root {root} {L.COMMON.format(cap="--drag-cap 10 --no-ground --first-solve-tilt 1", amin=L.AMIN)} '
           f'--baseline {grid}/baseline.json --check-design {grid}/offset_design.json')
    cores = '12-35' if lane == 'A' else '36-59'
    (root/'command.txt').write_text(
        cmd + f"\n# {name}: as run_016_area050 r00 ({L.SETUP}, no drag or loss cap), out of ground effect: no ground, the domain "
        f"goes on to a far field 10 chords below y = 0 (remesh_design.py --free-air; --no-ground: no y = 0 symmetry plane in the "
        f"mesh warp, outlet loss plane over y = -5..5 m). The box is unchanged. First solve from freestream (no initial flow fits "
        f"this grid) at the start tilted 1 deg (--first-solve-tilt 1), then walked to the level start: in the 4-rank test the "
        f"level start's symmetric zero-lift flow, solved directly from freestream, stalled near 2e-9 of its initial residual. 24 ranks on cores {cores}, lane {lane} of comparisons/freeair_20261121/freeair_queue.py, {time.strftime('%c')}. "
        f"Folder {real}\n")
    return link, root

if __name__ == '__main__':
    FQ.log('queue:', NAME)
    while AFTER[1] not in AFTER[0].read_text() or not READY.exists(): time.sleep(60)
    lane = None
    while lane is None:  # the first lane to come free
        for k in 'AB':
            if not (FQ.HERE/f'lane_{k}.busy').exists() and subprocess.run(['pgrep', '-f', f'lanes.py {k}$'], capture_output=True).returncode:
                lane = k; break
        else: time.sleep(60)
    FQ.claim(lane, NAME)
    FQ.log(f'{NAME}: on lane {lane}')
    try:
        link, root = build(NAME, lane)
        FQ.log(f'{NAME}: started r00 (free air, --no-ground)')
        FQ.log(f'{NAME}: {FQ.chain(link, [root], lane)}')
        FQ.finish(link, [root], NAME)
    except Exception as e: FQ.log(f'{NAME}: queue error {type(e).__name__}: {e}')
    finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log('done')

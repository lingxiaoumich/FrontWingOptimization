"""run_035_area050FreeAir: retry of its second grid, on the next free lane ahead of priority_queue.py (run_035 is the
user's first priority, 2026-11-21).
r00 went well (C_L c -0.001 -> -1.120 in 45 evaluations, ~1 min each, until the warp gave negative volumes). r01, the
remesh of r00 evaluation 41, died at its first solve: from r00's flow mapped onto the new grid, NK stalled at ~1.4e-9 of
the initial residual for 50,000 cycles (target 1e-10) with C_L c settled at -1.122, like the level start from freestream
in test C. In free air only such rough starts stall: solves from a nearby converged flow reach ~1e-11, and the 1-deg
tilted solve from freestream converged in r00. So r01b remeshes the same iterate (same pyHyp setting, so the same
grid as r01) and starts with no initial flow: --first-solve-tilt 1 solves it tilted from freestream, then walks to it.
Later grids of the chain keep the mapped initial flow, and the driver (since 2026-11-22 04:40) falls back to the
tilted start from freestream if that first solve stalls.
Lane: two 'lanes.py A'/'lanes.py B' hold processes stop priority_queue.py from taking a lane (its lane_free() pgreps
for them) until this script has the lane lock of the first lane to come free; the holds are then dropped.
  setsid nohup taskset -c 12-63 python -u retry_035.py > retry_035.log 2>&1 &
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # chain(), finish(), log(), RANKFILE and the lane locks; FQ.RC is run_chain, FQ.L lanes.py
NAME = 'run_035_area050FreeAir'

def take_lane():
    holds = [subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(1e9)', 'lanes.py', k]) for k in 'AB']
    try:
        while True:
            for k in 'AB':
                try:
                    fd = os.open(FQ.HERE/f'lane_{k}.busy', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                    os.write(fd, f'{NAME} r01b {os.getpid()}\n'.encode()); os.close(fd); return k
                except FileExistsError: pass
            time.sleep(20)
    finally:
        for h in holds: h.kill()

if __name__ == '__main__':
    FQ.log(f'{NAME}: r01b waits for the next free lane')
    lane = take_lane()
    FQ.log(f'{NAME}: r01b on lane {lane}')
    try:
        link = FQ.L.O/NAME; r00 = link/'r00'
        new = FQ.RC.prepare(r00, link, 'r01b', [], 0)
        if new is None: raise RuntimeError('no usable grid for r01b')
        root, _, hyp = new
        (root/'source/restart_vol.cgns').unlink()
        shutil_rank = FQ.RANKFILE[lane]
        subprocess.run(['cp', str(shutil_rank), str(root/'rankfile')], check=True)
        with (root/'command.txt').open('a') as f:
            f.write(f"# retry of r01 (its first solve, from r00's mapped flow, stalled near 1.4e-9): the same remesh, no initial "
                    f"flow, so --first-solve-tilt 1 starts it tilted from freestream and walks to it; pyHyp {hyp}, lane {lane} "
                    f"of comparisons/priority_20261121/retry_035.py, {time.strftime('%c')}\n")
        FQ.log(f'{NAME}: started r01b (pyHyp {hyp}, no initial flow)')
        FQ.log(f'{NAME}: {FQ.chain(link, [r00, root], lane, append=True)}')
        FQ.finish(link, [r00, root], NAME)
    except Exception as e: FQ.log(f'{NAME}: retry error {type(e).__name__}: {e}')
    finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log('done')

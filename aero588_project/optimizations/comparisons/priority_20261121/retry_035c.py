"""run_035_area050FreeAir: second retry of its second grid (r01c), on the next free lane ahead of priority_queue.py
(run_035 is the user's first priority, 2026-11-21).
r01 and r01b (both the remesh of r00 evaluation 41, the same grid) died at their first solve, which never met
ADflow's target of 1e-10 of the grid's freestream residual:
  r01  from r00's flow mapped onto the grid: below 3e-9 after 29,700 steps, at C_L c -1.122000, then between 7e-10
       and 1.6e-6 (spikes) to the 60,000-step limit (-1.121995 at the end);
  r01b from freestream, tilted 1 deg (retry_035.py): below 1.5e-9 after 34,800 steps (C_L c -1.153594), then between
       1.2e-9 and 5e-4 until run_chain.py stopped it after 60 min with no evaluation.
So in free air at this downforce the residual stalls near 1e-9 whatever the start, while the forces settle to ~5e-6.
r01c remeshes the same iterate (same pyHyp setting, same grid), starts from the mapped flow like r01, and runs with
--l2-convergence 3e-9: a solve ends the first time the residual is below 3e-9 of the freestream residual. Later
grids of the chain keep it (run_chain.py copies the command). --first-solve-tilt 1 stays as the fallback.
The derivative check's 5 % agreement is far above the force noise this leaves (~5e-6 in C_L c).
Lane: as retry_035.py, two 'lanes.py A'/'lanes.py B' hold processes stop priority_queue.py from taking a lane until
this script has the lane lock of the first lane to come free.
  setsid nohup taskset -c 12-63 python -u retry_035c.py > retry_035c.log 2>&1 &
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
L2 = '3e-9'

def take_lane():
    holds = [subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(1e9)', 'lanes.py', k]) for k in 'AB']
    try:
        while True:
            for k in 'AB':
                try:
                    fd = os.open(FQ.HERE/f'lane_{k}.busy', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                    os.write(fd, f'{NAME} r01c {os.getpid()}\n'.encode()); os.close(fd); return k
                except FileExistsError: pass
            time.sleep(20)
    finally:
        for h in holds: h.kill()

if __name__ == '__main__':
    FQ.log(f'{NAME}: r01c waits for the next free lane')
    lane = take_lane()
    FQ.log(f'{NAME}: r01c on lane {lane}')
    try:
        link = FQ.L.O/NAME; r00 = link/'r00'
        new = FQ.RC.prepare(r00, link, 'r01c', [], 0)
        if new is None: raise RuntimeError('no usable grid for r01c')
        root, _, hyp = new
        subprocess.run(['cp', str(FQ.RANKFILE[lane]), str(root/'rankfile')], check=True)
        lines = (root/'command.txt').read_text().split('\n')
        lines[0] += f' --l2-convergence {L2}'
        lines.insert(-1, f"# retry of r01 and r01b (their first solves stalled near 1e-9 of the freestream residual, the target being "
                         f"1e-10, with the forces settled): the same remesh, from r00's mapped flow, with --l2-convergence {L2}; "
                         f"pyHyp {hyp}, lane {lane} of comparisons/priority_20261121/retry_035c.py, {time.strftime('%c')}")
        (root/'command.txt').write_text('\n'.join(lines))
        FQ.log(f'{NAME}: started r01c (pyHyp {hyp}, mapped flow, --l2-convergence {L2})')
        FQ.log(f'{NAME}: {FQ.chain(link, [r00, root], lane, append=True)}')
        FQ.finish(link, [r00, root], NAME)
    except Exception as e: FQ.log(f'{NAME}: retry error {type(e).__name__}: {e}')
    finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log('done')

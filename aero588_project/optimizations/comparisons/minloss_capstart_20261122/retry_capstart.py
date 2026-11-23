"""Retries of the two CapStart rounds whose first solve stalled (2026-11-22), with --l2-convergence 3e-9 as run_035 r01c:
  run_037_area050MinLossClc08CpMin4CapStart r00 (start: run_030 r00 eval 30) -> r00b, the same remesh
  run_032_area050MinLossClc16CapStart       r01 (remesh of its r00 eval 23)   -> r01b, the same remesh
Both first solves, from the mapped flow, ended at 60,000 steps between 8e-10 and 1.3e-9 of the grid's freestream
residual (target 1e-10), with C_L c settled to 6 digits and spikes up to 1e-3 of it: the stall seen in free air
(run_035 r01/r01b). With 3e-9, run_035 r01c and r02 ran through (derivative checks within 6-7 %, a warning, vs
0.009 % at 1e-10). Later grids of these chains keep the option (run_chain.py copies the command).
  setsid nohup taskset -c 12-63 python -u retry_capstart.py <name> <lane> >> retry_capstart.log 2>&1 &
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import minloss_capstart_queue as Q
FQ, L, RC = Q.FQ, Q.L, Q.RC
L2 = '3e-9'
JOBS = {'run_037_area050MinLossClc08CpMin4CapStart': (-0.8, True, 'r00b'),
        'run_032_area050MinLossClc16CapStart': (-1.6, False, 'r01b')}

if __name__ == '__main__':
    name, lane = sys.argv[1], sys.argv[2]
    t, cp, rnd = JOBS[name]
    fd = os.open(FQ.HERE/f'lane_{lane}.busy', os.O_CREAT | os.O_EXCL | os.O_WRONLY)  # fails if the lane is taken
    os.write(fd, f'{name} {rnd} {os.getpid()}\n'.encode()); os.close(fd)
    FQ.log(f'{name}: {rnd} on lane {lane}')
    try:
        link = L.O/name
        why = (f"# retry: the first solve of the round before it, from the mapped flow, stalled near 1e-9 of the freestream residual "
               f"(target 1e-10) with the forces settled; the same remesh with --l2-convergence {L2}; lane {lane} of "
               f"comparisons/minloss_capstart_20261122/retry_capstart.py, {time.strftime('%c')}")
        if rnd == 'r00b':
            _, root, what = Q.build(name, t, cp, lane, rnd=rnd, l2=L2)
            rounds = [root]
        else:
            r00 = link/'r00'
            new = RC.prepare(r00, link, rnd, [], 0)
            if new is None: raise RuntimeError(f'no usable grid for {rnd}')
            root, _, hyp = new
            if (link/'grids'/f'{rnd}_from_r00_eval023'/'source.txt').read_text().split()[-1] != '23': raise RuntimeError('not evaluation 23')
            subprocess.run(['cp', str(FQ.RANKFILE[lane]), str(root/'rankfile')], check=True)
            what = f'r00 evaluation 23, pyHyp {hyp}'
            rounds = [r00, root]
        lines = (root/'command.txt').read_text().split('\n')
        if '--l2-convergence' not in lines[0]: lines[0] += f' --l2-convergence {L2}'
        lines.insert(-1, why)
        (root/'command.txt').write_text('\n'.join(lines))
        FQ.log(f'{name}: started {rnd} ({what}, --l2-convergence {L2})')
        FQ.log(f'{name}: {FQ.chain(link, rounds, lane, append=True)}')
        Q.finish(link, name, t, cp, rounds)
    except Exception as e: FQ.log(f'{name}: retry error {type(e).__name__}: {e}')
    finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log(f'{name}: done')

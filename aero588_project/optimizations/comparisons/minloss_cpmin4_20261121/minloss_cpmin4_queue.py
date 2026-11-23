"""Queued after runs 031-034 and run_035 (user, 2026-11-21: "add on another 4 runs, same as 31 to 34, but limit to
Cp>=-4"): runs 031-034 with run_027_area050CpMin4's pressure limit added:
  run_036_area050MinLossClc18CpMin4   minimize the outlet Cpt loss, C_L c <= -1.8, wing surface Cp >= -4
  run_037_area050MinLossClc16CpMin4   C_L c <= -1.6, Cp >= -4
  run_038_area050MinLossClc14CpMin4   C_L c <= -1.4, Cp >= -4
  run_039_area050MinLossClc12CpMin4   C_L c <= -1.2, Cp >= -4
Otherwise as runs 031-034 (run_016/021's setup, level start, no drag cap). The chain merit is unchanged (loss + 0.1 x
the C_L c shortfall; the Cp limit is held by PSQP, as in runs 025-027). Each job may find no strictly feasible design
if Cp >= -4 cannot carry that much downforce: run_027 (max downforce with Cp >= -4) reached C_L c -1.393 at r00.
Not started on its own: comparisons/priority_20261121/priority_queue.py runs these jobs with run() (user's priority
order, 2026-11-21: 035, 034, 039, then the rest), after READY here, written by the scratchpad test D
(--max-clc -1.4 --cp-min -4 on run_021 r01's grid). The worker/claim code below is the unused stand-alone queue.
"""
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # chain() and the lane locks; imports lanes.py as FQ.L
L = FQ.L
READY = HERE/'READY'
AFTER = (HERE.parent/'minloss_20261121/queue.log', 'run_034_area050MinLossClc12: on lane')
FREEAIR = (HERE.parent/'freeair_20261121', 'run_035_area050FreeAir: on lane')
CP = -4
JOBS = [('run_036_area050MinLossClc18CpMin4', -1.8), ('run_037_area050MinLossClc16CpMin4', -1.6),
        ('run_038_area050MinLossClc14CpMin4', -1.4), ('run_039_area050MinLossClc12CpMin4', -1.2)]
CP_NOTE = ('--cp-min: ADflow cpmin, a KS aggregate over the wing surface with cpMinRho 100, at most ~0.06 below the true '
           'minimum; in runs 006-021 the lowest Cp of the whole field was always on the wing')
lock = threading.Lock()

def freeair_waiting():
    """run_035 is ready and waiting for a lane: it goes first."""
    log = FREEAIR[0]/'queue.log'
    return ((FREEAIR[0]/'READY').exists() and FREEAIR[1] not in log.read_text()
            and not subprocess.run(['pgrep', '-f', 'freeair_queue.py'], capture_output=True).returncode)

def claim(lane, name):
    """FQ.claim, but yielding to run_035 while it waits."""
    while True:
        if not freeair_waiting() and subprocess.run(['pgrep', '-f', f'lanes.py {lane}$'], capture_output=True).returncode:
            try:
                fd = os.open(FQ.HERE/f'lane_{lane}.busy', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, f'{name} {os.getpid()}\n'.encode()); os.close(fd); return
            except FileExistsError: pass
        time.sleep(60)

def finish(link, name, t):
    """Final chain plots and the best strictly feasible design: the lowest loss with C_L c <= t and Cp >= CP."""
    text = (link/'summary/run_chain.log').read_text()
    rounds = [link/'r00'] + [link/l.split()[-1] for l in text.splitlines() if ' started r' in l]
    rounds = [r for r in dict.fromkeys(rounds) if (r/'history.jsonl').exists()]
    time.sleep(120)  # live_plots.py makes its last refresh when the run stops
    with (link/'summary/final_plots.log').open('w') as f:
        subprocess.run([FQ.PY, str(FQ.S/'plot_run_chain.py'), '--out', str(link/'summary'), *map(str, rounds)], stdout=f, stderr=subprocess.STDOUT)
    rows = [json.loads((r/'best_feasible/result.json').read_text()) | {'round': r.name} for r in rounds if (r/'best_feasible/result.json').exists()]
    b = min(rows, key=lambda q: q['cpt_loss_outlet']) if rows else None
    if b: FQ.log(f"{name}: best loss {b['cpt_loss_outlet']:.6f} at C_L c {.5*b['cl']:.3f} (limit {t:g}), Cp min (KS) "
                 f"{b.get('cp_min_ks', float('nan')):.2f} (limit {CP}), C_D c {.5*b['cd']:.4f}, area {b['section_area_m2']/L.A0:.1%} "
                 f"at {b['round']} evaluation {b['evaluation']}")
    else: FQ.log(f'{name}: no strictly feasible design')

def run(name, t, lane):
    note = (f'as run_031-034 (run_021_area050CptLoss025 r00, {L.SETUP}, no drag cap, the outlet Cpt loss minimized with '
            f'C_L c <= {t:g}: --max-clc, downforce at least {-t:g}) with the wing surface Cp >= {CP} of run_027_area050CpMin4 '
            f'({CP_NOTE})')
    root, _ = L.make(name, lane, f'--drag-cap 10 --max-clc {t:g} --cp-min {CP}', note)
    text = (root/'command.txt').read_text().replace('comparisons/area050_20261120/lanes.py', 'comparisons/minloss_cpmin4_20261121/minloss_cpmin4_queue.py')
    (root/'command.txt').write_text(text)
    FQ.log(f'{name}: started r00 (--max-clc {t:g} --cp-min {CP})')
    FQ.log(f'{name}: {FQ.chain(root.parent, [root], lane)}')
    finish(root.parent, name, t)

def worker(lane):
    while AFTER[1] not in AFTER[0].read_text() or not READY.exists(): time.sleep(60)
    while True:
        with lock:
            if not JOBS: break
            name, t = JOBS[0]
        claim(lane, name)
        with lock:
            if not JOBS or JOBS[0][0] != name:  # the other lane took it meanwhile
                (FQ.HERE/f'lane_{lane}.busy').unlink(); continue
            JOBS.pop(0)
        FQ.log(f'{name}: on lane {lane}')
        try: run(name, t, lane)
        except Exception as e: FQ.log(f'{name}: queue error {type(e).__name__}: {e}')
        finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log(f'lane {lane} done')

if __name__ == '__main__':
    FQ.log('queue:', ', '.join(j[0] for j in JOBS))
    threads = [threading.Thread(target=worker, args=(lane,)) for lane in 'AB']
    for t in threads: t.start()
    for t in threads: t.join()

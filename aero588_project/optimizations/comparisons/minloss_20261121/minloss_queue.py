"""Queued after runs 025-030 (comparisons/cpmin_20261121/next_queue.py): four runs on run_021_area050CptLoss025's
setup with the roles of downforce and loss swapped (user, 2026-11-21): minimize the outlet Cpt loss with C_L c held
at or below a limit (--max-clc), instead of maximizing downforce under a loss cap:
  run_031_area050MinLossClc18   C_L c <= -1.8 (downforce at least 1.8)
  run_032_area050MinLossClc16   C_L c <= -1.6
  run_033_area050MinLossClc14   C_L c <= -1.4
  run_034_area050MinLossClc12   C_L c <= -1.2
Otherwise as run_016/021 r00: level NACA 0012 start on the 50 mm floor, 20 FFD pairs moving only up and down, tilt
locked, the reduced box, section area >= 50 % of the start, run_006 r00's initial flow, no drag cap. Each is a full
chain; run_chain.py judges a loss run on the merit loss + 0.1 x its C_L c shortfall (see run_chain.objective()).
Lanes are shared with floor_queue.py and next_queue.py (lane_A/B.busy in comparisons/floor_20261121); a lane is taken
only once run_030, the last job of next_queue.py, has started, and after READY (written once the driver tests passed).
  setsid nohup taskset -c 12-63 python -u minloss_queue.py > queue.log 2>&1 &
"""
import json
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # chain(), claim() and the lane locks; imports lanes.py as FQ.L
L = FQ.L
READY = HERE/'READY'
AFTER = (HERE.parent/'cpmin_20261121/queue.log', 'run_030_area050CptLoss015: on lane')
JOBS = [('run_031_area050MinLossClc18', -1.8), ('run_032_area050MinLossClc16', -1.6),
        ('run_033_area050MinLossClc14', -1.4), ('run_034_area050MinLossClc12', -1.2)]
lock = threading.Lock()

def finish(link, name, t):
    """Final chain plots and the best strictly feasible design: the lowest loss with C_L c <= t."""
    text = (link/'summary/run_chain.log').read_text()
    rounds = [link/'r00'] + [link/l.split()[-1] for l in text.splitlines() if ' started r' in l]
    rounds = [r for r in dict.fromkeys(rounds) if (r/'history.jsonl').exists()]
    time.sleep(120)  # live_plots.py makes its last refresh when the run stops
    with (link/'summary/final_plots.log').open('w') as f:
        subprocess.run([FQ.PY, str(FQ.S/'plot_run_chain.py'), '--out', str(link/'summary'), *map(str, rounds)], stdout=f, stderr=subprocess.STDOUT)
    rows = [json.loads((r/'best_feasible/result.json').read_text()) | {'round': r.name} for r in rounds if (r/'best_feasible/result.json').exists()]
    b = min(rows, key=lambda q: q['cpt_loss_outlet']) if rows else None
    if b: FQ.log(f"{name}: best loss {b['cpt_loss_outlet']:.6f} at C_L c {.5*b['cl']:.3f} (limit {t:g}), C_D c {.5*b['cd']:.4f}, "
                 f"area {b['section_area_m2']/L.A0:.1%} at {b['round']} evaluation {b['evaluation']}")
    else: FQ.log(f'{name}: no strictly feasible design')

def run(name, t, lane):
    note = (f'as run_021_area050CptLoss025 r00 ({L.SETUP}, no drag cap) with the roles of downforce and loss swapped: the outlet Cpt '
            f'loss (mass-averaged (p0_inf - p0)/q_inf on the plane x = 7.3 m over the full domain height) is minimized with '
            f'C_L c <= {t:g} (--max-clc: downforce at least {-t:g}) instead of maximizing downforce under a loss cap')
    root, _ = L.make(name, lane, f'--drag-cap 10 --max-clc {t:g}', note)
    text = (root/'command.txt').read_text().replace('comparisons/area050_20261120/lanes.py', 'comparisons/minloss_20261121/minloss_queue.py')
    (root/'command.txt').write_text(text)
    FQ.log(f'{name}: started r00 (--max-clc {t:g})')
    FQ.log(f'{name}: {FQ.chain(root.parent, [root], lane)}')
    finish(root.parent, name, t)

def worker(lane):
    while AFTER[1] not in AFTER[0].read_text() or not READY.exists(): time.sleep(60)
    while True:
        with lock:
            if not JOBS: break
            name, t = JOBS[0]
        FQ.claim(lane, name)
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

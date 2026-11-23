"""Queued after the lower-floor runs (comparisons/floor_20261121): six more runs on run_016_area050's setup
(level NACA 0012 start on the 50 mm floor, 20 FFD pairs moving only up and down, tilt locked, the reduced box,
section area >= 50 % of the start, run_006 r00's initial flow), each a full chain:
  run_025_area050CpMin8   wing surface Cp >= -8   (--cp-min, ADflow's KS cpmin; the field minimum sits on the wing)
  run_026_area050CpMin6   Cp >= -6
  run_027_area050CpMin4   Cp >= -4
  run_028_area050DragCap005   C_D c <= 0.05 (--drag-cap 0.10)
  run_029_area050DragCap002   C_D c <= 0.02 (--drag-cap 0.04)
  run_030_area050CptLoss015   outlet Cpt loss <= 0.15 x run_016's best (0.020445 -> 0.003067)
Lanes are shared with floor_queue.py (its lane_A/B.busy locks); a lane is taken only once run_024, the last floor
run, has started, and the Cp runs wait for READY (written once the cpmin derivative test passed).
  setsid nohup taskset -c 12-63 python -u next_queue.py > queue.log 2>&1 &
"""
import json
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # chain(), finish(), claim() and its lane locks; imports lanes.py as FQ.L
L = FQ.L
READY = HERE/'READY'
base = json.loads((HERE.parent/'area050_20261120/base_best.json').read_text())['cpt_loss_outlet']
JOBS = [('run_025_area050CpMin8', '--drag-cap 10 --cp-min -8', 'wing surface Cp >= -8'),
        ('run_026_area050CpMin6', '--drag-cap 10 --cp-min -6', 'wing surface Cp >= -6'),
        ('run_027_area050CpMin4', '--drag-cap 10 --cp-min -4', 'wing surface Cp >= -4'),
        ('run_028_area050DragCap005', '--drag-cap 0.10', 'the drag limit C_D c <= 0.05 (--drag-cap 0.10 caps cd; C_D c = 0.5 cd)'),
        ('run_029_area050DragCap002', '--drag-cap 0.04', 'the drag limit C_D c <= 0.02 (--drag-cap 0.04 caps cd; C_D c = 0.5 cd)'),
        ('run_030_area050CptLoss015', f'--drag-cap 10 --loss-cap {.15*base:.6f}',
         f'an outlet Cpt loss cap instead of the drag cap: mass-averaged (p0_inf - p0)/q_inf on the plane x = 7.3 m over the '
         f'full domain height <= {.15*base:.6f} = 0.15 x run_016\'s best ({base:.9f})')]
CP_NOTE = (' (--cp-min: ADflow cpmin, a KS aggregate over the wing surface with cpMinRho 100, at most ~0.06 below the true '
           'minimum; in runs 006-021 the lowest Cp of the whole field was always on the wing, 1.7-3.3 below the ground\'s)')
lock = threading.Lock()

def run(name, cap, what, lane):
    note = f'as run_016_area050 r00 ({L.SETUP}, no drag or loss cap) with {what}' + (CP_NOTE if 'cp-min' in cap else '')
    root, _ = L.make(name, lane, cap, note)
    text = (root/'command.txt').read_text().replace('comparisons/area050_20261120/lanes.py', 'comparisons/cpmin_20261121/next_queue.py')
    (root/'command.txt').write_text(text)
    FQ.log(f'{name}: started r00 ({cap})')
    FQ.log(f'{name}: {FQ.chain(root.parent, [root], lane)}')
    FQ.finish(root.parent, [root], name)

def worker(lane):
    while 'run_024_area050Floor12p5mm: on lane' not in (FQ.HERE/'queue_floor.log').read_text(): time.sleep(60)
    while True:
        with lock:
            if not JOBS: break
            name, cap, what = JOBS[0]
        while 'cp-min' in cap and not READY.exists(): time.sleep(60)
        FQ.claim(lane, name)
        with lock:
            if not JOBS or JOBS[0][0] != name:
                (FQ.HERE/f'lane_{lane}.busy').unlink(); continue
            JOBS.pop(0)
        FQ.log(f'{name}: on lane {lane}')
        try: run(name, cap, what, lane)
        except Exception as e: FQ.log(f'{name}: queue error {type(e).__name__}: {e}')
        finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log(f'lane {lane} done')

if __name__ == '__main__':
    FQ.log('queue:', ', '.join(j[0] for j in JOBS))
    threads = [threading.Thread(target=worker, args=(lane,)) for lane in 'AB']
    for t in threads: t.start()
    for t in threads: t.join()

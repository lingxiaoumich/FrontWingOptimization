"""One queue in the user's priority order (2026-11-21: "priority orders: 35, 34, 34 + Cp>=-4, the rest"), replacing
the waiting queues minloss_20261121/minloss_queue.py (031-034), freeair_20261121/freeair_queue.py (035) and
minloss_cpmin4_20261121/minloss_cpmin4_queue.py (036-039), and the unstarted jobs 028-030 of
cpmin_20261121/next_queue.py. That process still runs 026/027; it skips 028-030 because each has a placeholder folder
in optimizations/ (its make() refuses an existing folder), which is removed here before the job starts.
jobs.json holds the order and is re-read every time a lane comes free, so it can be edited while this runs. A free lane
takes the first job not yet started whose test has passed:
  freeair    freeair_20261121/READY          (free-air driver test)       built by freeair_queue.build()
  minloss    minloss_20261121/READY          (--max-clc test)             run by minloss_queue.run()
  minlosscp  minloss_cpmin4_20261121/READY   (--max-clc + --cp-min test)  run by minloss_cpmin4_queue.run()
  next       none                                                        run by next_queue.run()
  minlosscap minloss_capstart_20261122/READY (start-from-design build)   run by minloss_capstart_queue.run()
             (2026-11-22: the min-loss runs from a loss-cap design; job fields t = the C_L c limit, cp = with Cp >= -4)
A job whose test has not passed yet is passed over, not waited for. Those modules are re-imported for each job, so
fixes to them reach jobs that have not started. Lanes are the shared lane_A/B.busy locks in comparisons/floor_20261121.
  setsid nohup taskset -c 12-63 python -u priority_queue.py > queue.log 2>&1 &      (--dry: show the next picks)
"""
import importlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).parent
C = HERE.parent
for d in ('floor_20261121', 'cpmin_20261121', 'minloss_20261121', 'freeair_20261121', 'minloss_cpmin4_20261121', 'minloss_capstart_20261122'):
    sys.path.insert(0, str(C/d))
import floor_queue as FQ   # chain(), finish(), log() and the lane locks; imports lanes.py as FQ.L
L = FQ.L
JOBS, STARTED = HERE/'jobs.json', HERE/'started.txt'
GATE = {'freeair': C/'freeair_20261121/READY', 'minloss': C/'minloss_20261121/READY',
        'minlosscp': C/'minloss_cpmin4_20261121/READY', 'next': None,
        'minlosscap': C/'minloss_capstart_20261122/READY'}
NQ_LOG = C/'cpmin_20261121/queue.log'
lock = threading.Lock()

def module(name):
    return importlib.reload(importlib.import_module(name))

def nq_done():
    """next_queue.py has no unstarted job left (it went through the placeholders), or has ended."""
    return ('run_030_area050CptLoss015: queue error' in NQ_LOG.read_text()
            or subprocess.run(['pgrep', '-f', 'python -u next_queue.py'], capture_output=True).returncode != 0)

def placeholder(name):
    d = L.O/name
    return d.is_dir() and not d.is_symlink() and [p.name for p in d.iterdir()] == ['PLACEHOLDER.txt']

def started():
    return set(STARTED.read_text().split()) if STARTED.exists() else set()

def ready(job):
    gate = GATE[job['kind']]
    if gate is not None and not gate.exists(): return False
    if job['kind'] == 'next' and placeholder(job['name']) and not nq_done(): return False
    return not (L.O/job['name']).exists() or placeholder(job['name'])

def pick():
    """(the first startable job, whether any job is left)."""
    done = started(); left = [j for j in json.loads(JOBS.read_text()) if j['name'] not in done]
    return next((j for j in left if ready(j)), None), bool(left)

def lane_free(lane):
    return (nq_done() and not (FQ.HERE/f'lane_{lane}.busy').exists()
            and subprocess.run(['pgrep', '-f', f'lanes.py {lane}$'], capture_output=True).returncode != 0)

def run(job, lane):
    name, kind = job['name'], job['kind']
    if kind == 'freeair':
        FA = module('freeair_queue')
        link, root = FA.build(name, lane)
        FQ.log(f'{name}: started r00 (free air, --no-ground)')
        FQ.log(f'{name}: {FQ.chain(link, [root], lane)}')
        FQ.finish(link, [root], name)
    elif kind == 'minloss': module('minloss_queue').run(name, job['t'], lane)
    elif kind == 'minlosscp': module('minloss_cpmin4_queue').run(name, job['t'], lane)
    elif kind == 'minlosscap': module('minloss_capstart_queue').run(name, job['t'], job.get('cp', False), lane)
    else:
        if placeholder(name): shutil.rmtree(L.O/name)
        NQ = module('next_queue')
        _, cap, what = next(j for j in NQ.JOBS if j[0] == name)
        NQ.run(name, cap, what, lane)

def worker(lane):
    while True:
        while not lane_free(lane): time.sleep(60)
        with lock:
            job, left = pick()
            if not left: break
            if job is None: job = 'wait'
            else:
                try:
                    fd = os.open(FQ.HERE/f'lane_{lane}.busy', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                    os.write(fd, f"{job['name']} {os.getpid()}\n".encode()); os.close(fd)
                except FileExistsError: job = 'wait'
                else:
                    with STARTED.open('a') as f: f.write(job['name'] + '\n')
        if job == 'wait': time.sleep(60); continue
        FQ.log(f"{job['name']}: on lane {lane}")
        try: run(job, lane)
        except Exception as e: FQ.log(f"{job['name']}: queue error {type(e).__name__}: {e}")
        finally: (FQ.HERE/f'lane_{lane}.busy').unlink()
    FQ.log(f'lane {lane} done')

if __name__ == '__main__':
    if '--dry' in sys.argv:
        done = started()
        for j in json.loads(JOBS.read_text()):
            g = GATE[j['kind']]
            print(f"{j['name']:38s} {j['kind']:9s} started={j['name'] in done} test_passed={g is None or g.exists()} startable_now={ready(j)}")
        print('next_queue done with its jobs:', nq_done(), '| lanes free:', {k: lane_free(k) for k in 'AB'})
        sys.exit()
    FQ.log('queue:', ', '.join(j['name'] for j in json.loads(JOBS.read_text())))
    threads = [threading.Thread(target=worker, args=(lane,)) for lane in 'AB']
    for t in threads: t.start()
    for t in threads: t.join()

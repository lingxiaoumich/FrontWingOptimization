"""Lower box floor campaign (2026-11-21), queued behind the area050 lanes (comparisons/area050_20261120).

Jobs, in order, each taken by the first free 24-rank lane (A: cores 12-35, B: 36-59):
  1. run_017_area050DragCap003: retry of its third grid. r02's first solve stalled (NK at ~1e-4 for 55,000
     cycles, C_L c -1.622 settled), so r02b remeshes r01's last iterate with the next pyHyp setting; if that
     also gives no iterate, r02c uses the iterate before with the first setting. The chain then goes on.
  2-4. run_016_area050's limits (area >= 50 % of the start, no drag or loss cap), level NACA 0012 start on a
     lower box floor: 37.5 mm (run_022), 25 mm (run_023), 12.5 mm (run_024), i.e. 25/50/75 % less clearance.
Grid family for floor h (remesh_design.py --ground-gap/--ground-spacing/--fine-height, kept by every remesh
of the chain through baseline.json 'mesh'): the original 25 mm near grid, its lines that would reach below h/2
shortened about their wall point (smoothly from line to line), and background cells 0.05 h high from the ground
up to h, growing 10 % a cell to the original 2.5 mm. Under the wing that keeps the original's 10 background
cells across the near grid's depth and 10 in the gap below it; elsewhere the grid is the original. At h = 50 mm
it is exactly the original grid, so run_016 is the 50 mm member. (A first try, the whole near grid h/2 deep,
gave 1438 orphans at 12.5 mm: too little overlap above the wing.) No initial flow fits these grids: the first solve starts
from freestream (--first-ncycles 60000).
  setsid nohup taskset -c 12-63 python -u floor_queue.py [retry|floor] > queue_<kind>.log 2>&1 &
A lane is held by lane_<A|B>.busy while a job of any copy runs on it.
"""
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'area050_20261120'))
sys.path.insert(0, '/home/lingxiao77/AER588/aero588_project/scripts')
import lanes as L          # paths, the run_016 command, chain helpers
import run_chain as RC     # prepare(), HYP, orphans(), stop()

O, S, BIG, START = L.O, L.S, L.BIG, L.START
PY = sys.executable
RANKFILE = {k: O/v/'r00/rankfile' for k, v in L.RANKS.items()}
JOBS = [('retry', None, 'run_017_area050DragCap003'),
        ('floor', .0375, 'run_022_area050Floor37p5mm'),
        ('floor', .025, 'run_023_area050Floor25mm'),
        ('floor', .0125, 'run_024_area050Floor12p5mm')]
lock = threading.Lock()
def log(*w):
    with lock: L.log(*w)

def chain(link, rounds, lane, append=False):
    """Start the last round's command (already written), then live plots and run_chain.py on the rounds;
    returns run_chain.py's last line once the chain has ended."""
    root = rounds[-1]
    cmd = (root/'command.txt').read_text().split('\n')[0]
    subprocess.Popen(cmd, shell=True, cwd=root, stdout=(root/'run.log').open('w'), stderr=subprocess.STDOUT, start_new_session=True)
    t = time.time()
    while not L.ranks(root) and time.time()-t < 600: time.sleep(10)
    if not L.ranks(root): return 'the ranks never came up; see run.log'
    n = RC.orphans(root)
    if n:
        RC.stop(root, f'{n} orphan cells in the overset grid')
        return f'{n} orphan cells in the overset grid'
    time.sleep(30)
    subprocess.Popen([PY, '-u', str(S/'live_plots.py'), '--out', str(link/'summary'), *map(str, rounds)],
                     stdout=(link/'summary/live_plots.log').open('a'), stderr=subprocess.STDOUT, start_new_session=True)
    with (link/'summary/run_chain.log').open('a' if append else 'w') as f:
        subprocess.run([PY, '-u', str(S/'run_chain.py'), '--live-plots', str(link/'summary'), '--stall-minutes', '15', *map(str, rounds)],
                       stdout=f, stderr=subprocess.STDOUT)
    return (link/'summary/run_chain.log').read_text().strip().splitlines()[-1]

def finish(link, rounds, name):
    """Final chain plots on the rounds that form the chain, and the best strictly feasible design."""
    text = (link/'summary/run_chain.log').read_text()
    rounds = rounds + [link/l.split()[-1] for l in text.splitlines() if ' started r' in l and link/l.split()[-1] not in rounds]
    rounds = [r for r in rounds if (r/'history.jsonl').exists()]
    time.sleep(120)  # live_plots.py makes its last refresh when the run stops
    with (link/'summary/final_plots.log').open('w') as f:
        subprocess.run([PY, str(S/'plot_run_chain.py'), '--out', str(link/'summary'), *map(str, rounds)], stdout=f, stderr=subprocess.STDOUT)
    rows = [json.loads((r/'best_feasible/result.json').read_text()) | {'round': r.name} for r in rounds if (r/'best_feasible/result.json').exists()]
    b = min(rows, key=lambda q: q['cl']) if rows else None
    if b: log(f"{name}: best C_L c {.5*b['cl']:.3f}, C_D c {.5*b['cd']:.4f}, loss {b['cpt_loss_outlet']:.6f}, area "
              f"{b['section_area_m2']/L.A0:.1%}, clearance {b['minimum_clearance_m']*1e3:.2f} mm at {b['round']} evaluation {b['evaluation']}")
    else: log(f'{name}: no strictly feasible design')

def retry(name, lane):
    link = O/name; r01 = link/'r01'
    for rname, skip, back in [('r02b', [RC.HYP[0]], 0), ('r02c', [], 1)]:
        new = RC.prepare(r01, link, rname, skip, back)
        if new is None: log(f'{name}: no usable grid for {rname}'); continue
        root, _, hyp = new
        shutil.copy(RANKFILE[lane], root/'rankfile')
        with (root/'command.txt').open('a') as f:
            f.write(f"# retry of r02 (its first solve stalled) with pyHyp {hyp}, lane {lane} of comparisons/floor_20261121/floor_queue.py, {time.strftime('%c')}\n")
        log(f'{name}: started {rname} (pyHyp {hyp}, {back} iterate(s) before the newest)')
        end = chain(link, [link/'r00', r01, root], lane, append=True)
        log(f'{name}: {end}')
        if list((root/'iterates').glob('eval_*/design.json')): break
    finish(link, [link/'r00', r01, root], name)  # the retry round is started here, so run_chain.log does not list it

def floor(name, h, lane):
    link, real = O/name, BIG/name
    if link.exists(): raise RuntimeError(f'{link} exists already')
    real.mkdir(parents=True); link.symlink_to(real)
    grid = link/f'grids/r00_level_floor{h*1e4:.0f}'; grid.mkdir(parents=True); (link/'summary').mkdir()
    design = json.loads((START/'grids/r00_level_start/source_design.json').read_text())
    design['shape_y'] = [v-(.05-h) for v in design['shape_y']]  # every control point down: a rigid shift
    (grid/'source_design.json').write_text(json.dumps(design)+'\n')
    (grid/'source.txt').write_text(f'run_006 level start (grids/r00_level_start/source_design.json) lowered by {(.05-h)*1e3:g} mm to the {h*1e3:g} mm floor\n')
    gap = h/2
    made = subprocess.run([PY, str(S/'remesh_design.py'), '--design', str(grid/'source_design.json'), '--out', str(grid),
                           '--hyp-options', RC.HYP[0], '--ground-gap', f'{gap:g}', '--ground-spacing', f'{.05*h:g}', '--fine-height', f'{h:g}'],
                          capture_output=True, text=True, env=RC.ONE_THREAD, timeout=600)
    (grid/'remesh.log').write_text(made.stdout+made.stderr)
    if made.returncode: log(f'{name}: the start grid failed; see {grid}/remesh.log'); return
    (grid/'offset_design.json').write_text(json.dumps({'incidence': [-3.02], 'shape': [0.]*40}))
    root = link/'r00'; src = root/'source'; src.mkdir(parents=True)
    for f in L.INPUTS:
        if f not in ('restart_vol.cgns', 'ground.cgns', 'airfoil_coordinates.csv'): shutil.copy(START/'r00/source'/f, src/f)
    for f in ('ground.cgns', 'airfoil_coordinates.csv'): shutil.copy(grid/f, src/f)
    for f in (START/'r00/source').glob('*.py'): shutil.copy(S/f.name, src/f.name)
    shutil.copy(RANKFILE[lane], root/'rankfile')
    common = L.COMMON.format(cap='--drag-cap 10', amin=L.AMIN)
    assert '--box -0.0123 0.512 0.05 0.1881' in common
    common = common.replace('--box -0.0123 0.512 0.05 0.1881', f'--box -0.0123 0.512 {h:g} 0.1881')
    cmd = (f'mpirun -np 24 --map-by rankfile:file=rankfile --bind-to core --mca btl self,sm -x GFORTRAN_UNBUFFERED_PRECONNECTED=y '
           f'python -u {src}/run_slip_optimization.py --root {root} {common} '
           f'--baseline {grid}/baseline.json --check-design {grid}/offset_design.json')
    cores = '12-35' if lane == 'A' else '36-59'
    (root/'command.txt').write_text(
        cmd + f"\n# {name}: as run_016_area050 r00 (level NACA 0012 stretched to the box x-ends, 20 FFD pairs moving only up and down, "
        f"tilt locked, section area >= {L.AMIN:.9f} m^2 = 0.50 x the start's, no drag or loss cap), with the box floor at {h*1e3:g} mm "
        f"instead of 50 mm ({1-h/.05:.0%} less clearance) and the start on that floor. Grid family for the floor: the 25 mm near grid "
        f"with its lines shortened to stay {gap*1e3:g} mm above the ground, background cells {.05*h*1e3:g} mm high up to {h*1e3:g} mm, "
        f"growing 10 % a cell to 2.5 mm. First solve from freestream (no initial flow fits this grid). 24 ranks on cores {cores}, "
        f"lane {lane} of comparisons/floor_20261121/floor_queue.py, {time.strftime('%c')}. Folder on the 4T HD ({real}), linked from optimizations/\n")
    log(f'{name}: started r00 (floor {h*1e3:g} mm, near grid kept {gap*1e3:g} mm above the ground)')
    log(f'{name}: {chain(link, [root], lane)}')
    finish(link, [root], name)

def claim(lane, name):
    """Take the lane once the area050 lane driver has left it and no other copy of this queue holds it."""
    while True:
        if subprocess.run(['pgrep', '-f', f'lanes.py {lane}$'], capture_output=True).returncode:
            try:
                fd = os.open(HERE/f'lane_{lane}.busy', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, f'{name} {os.getpid()}\n'.encode()); os.close(fd); return
            except FileExistsError: pass
        time.sleep(60)

def worker(lane):
    while True:
        with lock:
            if not JOBS: break
            kind, h, name = JOBS[0]
        claim(lane, name)
        with lock:
            if not JOBS or JOBS[0][2] != name:  # the other lane took it meanwhile
                (HERE/f'lane_{lane}.busy').unlink(); continue
            JOBS.pop(0)
        log(f'{name}: on lane {lane}')
        try: retry(name, lane) if kind == 'retry' else floor(name, h, lane)
        except Exception as e: log(f'{name}: queue error {type(e).__name__}: {e}')
        finally: (HERE/f'lane_{lane}.busy').unlink()
    log(f'lane {lane} done')

if __name__ == '__main__':
    # optional arguments: the kinds of job to run (retry, floor); default all
    if sys.argv[1:]: JOBS[:] = [j for j in JOBS if j[0] in sys.argv[1:]]
    log('queue:', ', '.join(j[2] for j in JOBS))
    threads = [threading.Thread(target=worker, args=(lane,)) for lane in 'AB']
    for t in threads: t.start()
    for t in threads: t.join()

"""Area-floor campaign (2026-11-20): run_016_area050 = run_006_levelNoTilt with the section area held at
>= 50 % of the start's, then the cap sweeps of runs 008, 009, 012, 013, 014 again on that setup.

Two 24-rank lanes, each run a full chain (run_chain.py remeshes as before):
  lane A, cores 12-35: run_016_area050 -> run_019_area050CptLoss075 -> run_021_area050CptLoss025
  lane B, cores 36-59: run_017_area050DragCap003 -> run_018_area050DragCap008 -> run_020_area050CptLoss050
The loss caps are fractions of run_016's best design's outlet Cpt loss (as runs 012-014 were of
run_006's), so those runs wait for run_016's chain to end; lane A writes base_best.json then.
Every run starts from run_006's start (level NACA 0012 stretched to the box, run_006 r00's flow).

Runs 016-039 were made on the 4T HD and linked from optimizations/. Since 2026-11-23 the whole repo is on the 4T HD
(/home/lingxiao77/AER588 links to it), so BIG = O and new runs are plain folders. The real path has a space, so
everything here and in the chain uses the /home/lingxiao77 paths.
Run with the MACH-Aero env, detached, helpers on cores 12-63:
  setsid nohup taskset -c 12-63 python -u lanes.py A > lane_A.log 2>&1 &
"""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

P = Path('/home/lingxiao77/AER588/aero588_project')
O, S, HERE = P/'optimizations', P/'scripts', P/'optimizations/comparisons/area050_20261120'
BIG = O   # was '/media/lingxiao77/4T HD/WORK/PHD/AER588/optimizations' (linked runs) before the repo moved there
START = O/'run_012_cptLoss075'   # run_006's start: the same grid, coordinates and initial flow (checked)
A0 = 0.021549464275835817        # m^2, area of that start section
AMIN = 0.5*A0
BASE = 'run_016_area050'
RANKS = {'A': 'run_012_cptLoss075', 'B': 'run_013_cptLoss050'}  # their rankfiles: cores 12-35, 36-59
COMMON = ('--ffd-stations 20 --shape-directions y --fix-incidence {cap} --min-area {amin:.9f} --min-local-thickness 0.002 '
          '--te-min-angle 8 --te-min-thickness 0.0012 --box -0.0123 0.512 0.05 0.1881 --optimizer PSQP --maxit 500 '
          '--max-evaluations 2000 --xmax 1.0 --ncycles 10000 --nk-switch-tol 1e-6 --warp-err-tol 1e-5 '
          '--split-stalled-solves 4 --save-iterate-flows --iterate-volume-only --first-ncycles 60000')
SETUP = ('level NACA 0012 stretched to the box x-ends, 20 FFD pairs moving only up and down, tilt locked, the reduced box, '
         f'section area >= {AMIN:.9f} m^2 = 0.50 x the start\'s {A0:.9f} m^2 (--min-area)')
# name: (lane, cap flags or a loss fraction, note)
RUNS = {
    BASE: ('A', '--drag-cap 10', f'as run_006_levelNoTilt r00 ({SETUP.split(", section")[0]}, no drag cap), with the section area '
           f'held at >= {AMIN:.9f} m^2 = 0.50 x the start\'s {A0:.9f} m^2 (--min-area: area inside the 289 surface vertices)'),
    'run_017_area050DragCap003': ('B', '--drag-cap 0.06', f'as {BASE} r00 ({SETUP}) with the drag limit C_D c <= 0.03 '
                                  '(--drag-cap 0.06 caps cd; C_D c = 0.5 cd), the cap of run_008_dragCap003'),
    'run_018_area050DragCap008': ('B', '--drag-cap 0.16', f'as {BASE} r00 ({SETUP}) with the drag limit C_D c <= 0.08 '
                                  '(--drag-cap 0.16 caps cd; C_D c = 0.5 cd), the cap of run_009_dragCap008'),
    'run_019_area050CptLoss075': ('A', .75, ''), 'run_020_area050CptLoss050': ('B', .50, ''),
    'run_021_area050CptLoss025': ('A', .25, ''),
}
LANES = {'A': [BASE, 'run_019_area050CptLoss075', 'run_021_area050CptLoss025'],
         'B': ['run_017_area050DragCap003', 'run_018_area050DragCap008', 'run_020_area050CptLoss050']}
INPUTS = ['airfoil_coordinates.csv', 'ground.cgns', 'restart_vol.cgns', 'mesh_check.json', 'sha256_manifest.json',
          'software_versions.json']

def log(*w): print(time.strftime('%m-%d %H:%M:%S'), *w, flush=True)
def ranks(r00): return subprocess.run(['pgrep', '-f', f'[p]ython -u {r00}/source/'], capture_output=True).returncode == 0

def base_best():
    """Wait for lane A's verdict on run_016; returns its best design's record, or None if it has none."""
    while True:
        for f, ok in [(HERE/'base_best.json', True), (HERE/'base_failed.json', False)]:
            if f.exists(): return json.loads(f.read_text()) if ok else None
        time.sleep(60)

def chain_rounds(run):
    """r00 and every round run_chain.py started (failed orphan retries are not in the chain)."""
    lines = (run/'summary/run_chain.log').read_text().splitlines() if (run/'summary/run_chain.log').exists() else []
    return [run/'r00'] + [run/l.split()[-1] for l in lines if ' started r' in l]

def best_feasible(run):
    rows = [json.loads((r/'best_feasible/result.json').read_text()) | {'round': r.name}
            for r in chain_rounds(run) if (r/'best_feasible/result.json').exists()]
    return min(rows, key=lambda q: q['cl']) if rows else None

def make(name, lane, cap, note):
    link, real = O/name, BIG/name
    if link.exists(): raise RuntimeError(f'{link} exists already')
    real.mkdir(parents=True)
    if real != link: link.symlink_to(real)
    shutil.copytree(START/'grids/r00_level_start', link/'grids/r00_level_start')
    src = link/'r00/source'; src.mkdir(parents=True); (link/'summary').mkdir()
    for f in INPUTS: shutil.copy(START/'r00/source'/f, src/f)
    for f in (START/'r00/source').glob('*.py'): shutil.copy(S/f.name, src/f.name)  # the current scripts
    shutil.copy(O/RANKS[lane]/'r00/rankfile', link/'r00/rankfile')
    root, grid = link/'r00', link/'grids/r00_level_start'
    cmd = (f'mpirun -np 24 --map-by rankfile:file=rankfile --bind-to core --mca btl self,sm -x GFORTRAN_UNBUFFERED_PRECONNECTED=y '
           f'python -u {root}/source/run_slip_optimization.py --root {root} ' + COMMON.format(cap=cap, amin=AMIN) +
           f' --baseline {grid}/baseline.json --check-design {grid}/offset_design.json')
    cores = '12-35' if lane == 'A' else '36-59'
    (root/'command.txt').write_text(cmd + f"\n# {name}: {note}. Initial flow source/restart_vol.cgns = run_006 r00's (converged on this "
                                    f"grid). 24 ranks on cores {cores}, lane {lane} of comparisons/area050_20261120/lanes.py, "
                                    f"{time.strftime('%c')}. Folder {real}\n")
    return root, cmd

def run(name):
    lane, cap, note = RUNS[name]
    if isinstance(cap, float):
        b = base_best()
        if b is None: log(f'{name}: skipped, {BASE} has no best design'); return
        value = cap*b['cpt_loss_outlet']
        note = (f"as {BASE} r00 ({SETUP}), but the drag cap is replaced by an outlet Cpt loss cap: mass-averaged (p0_inf - p0)/q_inf "
                f"on the plane x = 7.3 m over the full domain height <= {value:.6f} = {cap:.2f} x {BASE}'s best "
                f"({b['cpt_loss_outlet']:.9f}, {b['round']} evaluation {b['evaluation']}), as {'run_012 run_013 run_014'.split()[[.75, .5, .25].index(cap)]} "
                f"did on run_006's")
        cap = f'--drag-cap 10 --loss-cap {value:.6f}'
    root, cmd = make(name, lane, cap, note)
    rlink = root.parent
    subprocess.Popen(cmd, shell=True, cwd=root, stdout=(root/'run.log').open('w'), stderr=subprocess.STDOUT, start_new_session=True)
    log(f'{name}: started r00 ({cap})')
    t = time.time()
    while not ranks(root) and time.time()-t < 600: time.sleep(10)
    if not ranks(root): log(f'{name}: the ranks never came up; see r00/run.log'); return
    time.sleep(30)
    py = sys.executable
    subprocess.Popen([py, '-u', str(S/'live_plots.py'), '--out', str(rlink/'summary'), str(root)],
                     stdout=(rlink/'summary/live_plots.log').open('a'), stderr=subprocess.STDOUT, start_new_session=True)
    with (rlink/'summary/run_chain.log').open('w') as f:
        subprocess.run([py, '-u', str(S/'run_chain.py'), '--live-plots', str(rlink/'summary'), '--stall-minutes', '15', str(root)],
                       stdout=f, stderr=subprocess.STDOUT)
    end = (rlink/'summary/run_chain.log').read_text().strip().splitlines()[-1]
    log(f'{name}: {end}')
    rounds = chain_rounds(rlink)
    time.sleep(120)  # live_plots.py makes its last refresh when the run stops
    with (rlink/'summary/final_plots.log').open('w') as f:
        subprocess.run([py, str(S/'plot_run_chain.py'), '--out', str(rlink/'summary'), *map(str, rounds)], stdout=f, stderr=subprocess.STDOUT)
    b = best_feasible(rlink)
    if b: log(f"{name}: best C_L c {.5*b['cl']:.3f}, C_D c {.5*b['cd']:.4f}, loss {b['cpt_loss_outlet']:.6f}, area "
              f"{b['section_area_m2']:.6f} m^2 ({b['section_area_m2']/A0:.1%}) at {b['round']} evaluation {b['evaluation']}")
    if name == BASE:
        out = {k: b[k] for k in ['round', 'evaluation', 'cl', 'cd', 'cpt_loss_outlet', 'cpt_loss_outlet_flux_c_m', 'section_area_m2']} if b else None
        (HERE/('base_best.json' if b else 'base_failed.json')).write_text(json.dumps(out or {'chain_end': end}, indent=1)+'\n')

if __name__ == '__main__':
    lane = sys.argv[1]
    log(f'lane {lane}:', ' -> '.join(LANES[lane]))
    for name in LANES[lane]:
        try: run(name)
        except Exception as e: log(f'{name}: lane driver error {type(e).__name__}: {e}')
        if name == BASE and not ((HERE/'base_best.json').exists() or (HERE/'base_failed.json').exists()):
            (HERE/'base_failed.json').write_text(json.dumps({'error': 'the lane driver stopped before a verdict'})+'\n')
    log(f'lane {lane} done')

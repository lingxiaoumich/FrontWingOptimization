"""Rotated-start campaign (2026-11-25): can a start with the trailing edge already high in the box avoid remeshing?

Start: the level start of runs 016-029 (NACA 0012 stretched to the box x-ends, run_006's grids/r00_level_start)
mapped by an exact affine map, a rotation and a uniform scale about its leading edge and then a shift. The map is
chosen so the section again spans the box x-ends (-0.0123 to 0.512 m), its trailing-edge midpoint sits at 2/3 of
the box height (y = 0.05 + 2/3 x 0.1381 = 0.1421 m), and its lowest point is on the 50 mm floor (rotation 8.28 deg
nose down, scale 1.0102). FFD is affine invariant, so the same map on every control point of the level start's
design gives exactly this section, and remesh_design.py builds its grid and lattice as for any start.
Two runs, each the limits of a campaign run on this start, at the same time on the two 24-rank lanes:
  lane A, cores 12-35: run_040_area050RotStart          = run_016_area050's limits (no drag cap)
  lane B, cores 36-59: run_041_area050DragCap002RotStart = run_029_area050DragCap002's limits (C_D c <= 0.02)
No initial flow fits the new grid: the first solve starts from freestream (--first-ncycles 60000, as before).
If the solver reports orphan cells on the start grid, the start is rebuilt with the next pyHyp setting as r00b, ...
run_chain.py then remeshes as in every chain (--stall-minutes 15, as runs 016-029).
Run with the MACH-Aero env, detached, helpers on cores 12-63:
  setsid nohup taskset -c 12-63 python -u rotstart_queue.py > queue.log 2>&1 &
"""
import json
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'area050_20261120'))
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
sys.path.insert(0, '/home/lingxiao77/AER588/aero588_project/scripts')
import lanes as L          # paths, the run_016 command, the start's inputs
import run_chain as RC     # HYP, orphans(), stop()
import floor_queue as F    # chain(), finish()

O, S = L.O, L.S
PY = sys.executable
LEVEL = O/'run_016_area050/grids/r00_level_start'
BOX = dict(xmin=-.0123, xmax=.512, ymin=.05, ymax=.1881)
TE_FRACTION = 2/3
RANKFILE = {k: O/v/'r00/rankfile' for k, v in L.RANKS.items()}
JOBS = {'A': ('run_040_area050RotStart', '--drag-cap 10', 'run_016_area050', 'no drag or loss cap'),
        'B': ('run_041_area050DragCap002RotStart', '--drag-cap 0.04', 'run_029_area050DragCap002',
              'the drag limit C_D c <= 0.02 (--drag-cap 0.04 caps cd; C_D c = 0.5 cd)')}
lock = threading.Lock()
def log(*w):
    with lock: L.log(*w)

def rotated_design():
    """The level start's design with every control point moved by the affine map; returns (design, map record)."""
    xy = np.loadtxt(LEVEL/'airfoil_coordinates.csv', delimiter=',', skiprows=1)
    P = xy[136].copy()  # leading-edge vertex
    te_y = BOX['ymin']+TE_FRACTION*(BOX['ymax']-BOX['ymin'])
    def amap(p, th, s, tx, ty):
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        return (p-P)@R.T*s+P+[tx, ty]
    def eqs(v):
        q = amap(xy, *v)
        return [q[:, 0].min()-BOX['xmin'], q[:, 0].max()-BOX['xmax'], (q[0, 1]+q[272, 1])/2-te_y, q[:, 1].min()-BOX['ymin']]
    v = fsolve(eqs, [np.deg2rad(9), 1., 0., 0.], xtol=1e-14)
    assert np.abs(eqs(v)).max() < 1e-12, eqs(v)
    lines = (LEVEL/'ffd.xyz').read_text().split('\n'); dims = [int(t) for t in lines[1].split()]
    vals = np.array(' '.join(lines[2:]).split(), float).reshape(3, -1)
    lat = np.stack([a.reshape(dims, order='F') for a in vals], -1)[:, :, 0, :2]  # (station, lower/upper, xy)
    new = amap(lat.reshape(-1, 2), *v).reshape(lat.shape)
    d = json.loads((LEVEL/'source_design.json').read_text())
    assert d['incidence'] == [-3.0] and len(d['shape_x']) == lat.shape[0]*2
    d['shape_x'] = (np.array(d['shape_x'])+(new[..., 0]-lat[..., 0]).ravel()).tolist()
    d['shape_y'] = (np.array(d['shape_y'])+(new[..., 1]-lat[..., 1]).ravel()).tolist()
    q = amap(xy, *v)
    record = {'rotation_deg': float(np.degrees(v[0])), 'scale': float(v[1]), 'shift_m': [float(v[2]), float(v[3])],
              'pivot_leading_edge_m': P.tolist(), 'te_mid_y_m': te_y, 'te_fraction_of_box_height': TE_FRACTION,
              'section_y_range_m': [float(q[:, 1].min()), float(q[:, 1].max())],
              'section_area_m2': float(.5*abs(q[:, 0]@np.roll(q[:, 1], -1)-q[:, 1]@np.roll(q[:, 0], -1)))}
    return d, record

def start_grid(grid, design, record, hyp):
    grid.mkdir(parents=True)
    (grid/'source_design.json').write_text(json.dumps(design)+'\n')
    (grid/'rotation_map.json').write_text(json.dumps(record, indent=1)+'\n')
    (grid/'source.txt').write_text(
        f"run_006 level start (run_016_area050/grids/r00_level_start/source_design.json) rotated {record['rotation_deg']:.3f} deg nose down "
        f"and scaled {record['scale']:.5f} about its leading edge, then shifted, so it spans the box x-ends with its trailing-edge "
        f"midpoint at y = {record['te_mid_y_m']:.4f} m (2/3 of the box height) and its lowest point on the 50 mm floor "
        f"(comparisons/rotstart_20261125/rotstart_queue.py)\n")
    made = subprocess.run([PY, str(S/'remesh_design.py'), '--design', str(grid/'source_design.json'), '--out', str(grid), '--hyp-options', hyp],
                          capture_output=True, text=True, env=RC.ONE_THREAD, timeout=600)
    (grid/'remesh.log').write_text(made.stdout+made.stderr)
    if made.returncode: return False
    (grid/'offset_design.json').write_text(json.dumps({'incidence': [-3.02], 'shape': [0.]*40}))
    return True

def setup_round(link, root, grid, lane, name, cap, like, what, hyp):
    src = root/'source'; src.mkdir(parents=True)
    for f in L.INPUTS:
        if f not in ('restart_vol.cgns', 'ground.cgns', 'airfoil_coordinates.csv'): shutil.copy(L.START/'r00/source'/f, src/f)
    for f in ('ground.cgns', 'airfoil_coordinates.csv'): shutil.copy(grid/f, src/f)
    for f in (L.START/'r00/source').glob('*.py'): shutil.copy(S/f.name, src/f.name)
    shutil.copy(RANKFILE[lane], root/'rankfile')
    common = L.COMMON.format(cap=cap, amin=L.AMIN)
    cmd = (f'mpirun -np 24 --map-by rankfile:file=rankfile --bind-to core --mca btl self,sm -x GFORTRAN_UNBUFFERED_PRECONNECTED=y '
           f'python -u {src}/run_slip_optimization.py --root {root} {common} '
           f'--baseline {grid}/baseline.json --check-design {grid}/offset_design.json')
    cores = '12-35' if lane == 'A' else '36-59'
    (root/'command.txt').write_text(
        cmd + f"\n# {name}: as {like} r00 (20 FFD pairs moving only up and down, tilt locked, the reduced box, section area >= "
        f"{L.AMIN:.9f} m^2, {what}), but from a rotated start: the level NACA 0012 start rotated nose down so its trailing edge "
        f"is at 2/3 of the box height and its lowest point on the floor (see {grid.name}/source.txt). Start grid pyHyp {hyp}. "
        f"First solve from freestream (no initial flow fits this grid). 24 ranks on cores {cores}, lane {lane} of "
        f"comparisons/rotstart_20261125/rotstart_queue.py, {time.strftime('%c')}. Folder {link}\n")
    return cmd

def job(lane, key=None):
    """Run JOBS[key] (default: the lane's own job) on `lane`."""
    name, cap, like, what = JOBS[key or lane]
    link = O/name
    if link.exists() and any(p.name != 'failed_attempts' for p in link.iterdir()): raise RuntimeError(f'{link} exists already')
    (link/'summary').mkdir(parents=True)
    design, record = rotated_design()
    log(f"{name}: rotated start {record['rotation_deg']:.3f} deg, scale {record['scale']:.5f}, TE mid y {record['te_mid_y_m']:.4f} m")
    for k, hyp in enumerate(RC.HYP):
        suffix = '' if k == 0 else chr(ord('a')+k)
        grid, root = link/f'grids/r00{suffix}_rot_start', link/f'r00{suffix}'
        if not start_grid(grid, design, record, hyp):
            log(f'{name}: start grid with pyHyp {hyp} failed; see {grid}/remesh.log'); continue
        cmd = setup_round(link, root, grid, lane, name, cap, like, what, hyp)
        subprocess.Popen(cmd, shell=True, cwd=root, stdout=(root/'run.log').open('w'), stderr=subprocess.STDOUT, start_new_session=True)
        log(f'{name}: started {root.name} (pyHyp {hyp})')
        t = time.time()
        while not L.ranks(root) and time.time()-t < 600: time.sleep(10)
        if not L.ranks(root): log(f'{name}: the ranks never came up; see {root.name}/run.log'); return
        n = RC.orphans(root)
        if n:
            RC.stop(root, f'{n} orphan cells in the overset grid; the adjoint cannot converge with them')
            log(f'{name}: {n} orphan cells on the start grid; rebuilding it with the next pyHyp setting'); continue
        # F.chain starts the round's command itself; it is running already, so give it a no-op start
        end = chain_running(link, root)
        log(f'{name}: {end}')
        F.finish(link, [root], name)
        return
    log(f'{name}: no pyHyp setting gave a usable start grid')

def chain_running(link, root):
    """live_plots.py and run_chain.py on the round that is already running; run_chain.py's last line."""
    time.sleep(30)
    subprocess.Popen([PY, '-u', str(S/'live_plots.py'), '--out', str(link/'summary'), str(root)],
                     stdout=(link/'summary/live_plots.log').open('a'), stderr=subprocess.STDOUT, start_new_session=True)
    with (link/'summary/run_chain.log').open('w') as f:
        subprocess.run([PY, '-u', str(S/'run_chain.py'), '--live-plots', str(link/'summary'), '--stall-minutes', '15', str(root)],
                       stdout=f, stderr=subprocess.STDOUT)
    return (link/'summary/run_chain.log').read_text().strip().splitlines()[-1]

def worker(lane):
    try: job(lane)
    except Exception as e: log(f'{JOBS[lane][0]}: queue error {type(e).__name__}: {e}')
    log(f'lane {lane} done')

def rerun_041():
    """run_041's first attempt (05:45) was ended after 20 evaluations by run_chain.py's stall rule: its start is 76 % over
    the drag cap, PSQP traded downforce for drag, and the rule saw no C_L c gain in 15 min. run_chain.py now counts only
    iterates within 1 % of the caps (user, 2026-11-25: "fix rule, rerun"). The attempt is kept in failed_attempts/; the
    rerun takes lane A once this queue's first run (run_040) has left it."""
    log_text = lambda: (HERE/'queue.log').read_text()
    while 'lane A done' not in log_text(): time.sleep(60)
    name = JOBS['B'][0]; link = O/name
    if (link/'r00').exists():
        keep = link/'failed_attempts'/'attempt_0545_stall_rule'; keep.mkdir(parents=True)
        for p in list(link.iterdir()):
            if p.name != 'failed_attempts': shutil.move(str(p), keep/p.name)
        (keep/'README.txt').write_text(rerun_041.__doc__.strip()+'\n')
    log(f'{name}: rerun on lane A')
    try: job('A', 'B')
    except Exception as e: log(f'{name}: queue error {type(e).__name__}: {e}')
    log('rerun done')

if __name__ == '__main__':
    if sys.argv[1:] == ['--rerun-041']:
        rerun_041(); sys.exit()
    lanes = sys.argv[1:] or ['A', 'B']
    log('queue:', ', '.join(JOBS[k][0] for k in lanes))
    threads = [threading.Thread(target=worker, args=(k,)) for k in lanes]
    for t in threads: t.start()
    for t in threads: t.join()

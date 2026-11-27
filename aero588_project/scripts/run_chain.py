"""Keep an optimization going across worn-out grids: watch the run, and when its grid needs
replacing, stop it, remesh at its last PSQP iterate and start the next run there.

An optimization lives in optimizations/run_NNN_<name>/: r00 is the run on the original
grid, r01, r02, ... the runs on each new grid (a letter, r03b, marks a retry), grids/ the
grids and summary/ the chain plots. Give the run folders in order; the last one is the
current run. The hand-over:
  1. stop the run and note why in its command.txt;
  2. remesh_design.py at its last saved iterate (iterates/eval_NNN) into
     grids/rNN_from_<run>_evalNNN, with the first pyHyp setting in HYP that gives
     no cell with a negative corner (remesh_attempts.json lists the tries);
  3. the next run folder rNN with the previous run's source/, the current
     run_slip_optimization.py, the new grid and the iterate's flow as the initial flow;
     the same command with the new --baseline;
  4. start it (it runs its own derivative check) and, with --live-plots, live_plots.py. If the
     solver reports orphan cells in the overset grid, stop it and go back to step 2 with the
     next pyHyp setting, as rNNb, rNNc, ...; when no setting gives a usable grid, with the
     iterate before (up to 3 back).

The grid needs replacing when
  - the last --max-failures evaluations all failed, or 6 of the last 8 did,
  - the best C_L c of the PSQP iterates improved by under --min-gain in --stall-minutes (for a run that
    minimizes the outlet loss under --max-clc: its merit, by under --min-loss-gain; for one that minimizes the
    suction peak, --max-clc --min-suction: its merit, by under --min-cp-gain; see objective()), counting only
    iterates within 1 % of every drag, loss or Cp cap (a run restoring feasibility is not stalled; see iterates()),
  - PSQP exited by itself (a fresh grid also restarts PSQP), or
  - no evaluation finished in --hang-minutes.
The chain ends when a run gained under --min-gain in all, or its best is not --min-gain better
than the best of the earlier runs (a fresh grid moves C_L c by some 2 %, which a run can regain
without real progress), a run died by itself, no pyHyp setting gives a usable grid, the disk is
nearly full, or --max-runs runs have started.
Killing this driver leaves the optimization running. Run with the MACH-Aero env:
  cd optimizations/run_NNN_<name>
  python ../../scripts/run_chain.py --live-plots summary r00 [r01 ...] > summary/run_chain.log 2>&1 &
--dry-run DIR only does steps 2-3 for the last run, under DIR.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

S = Path(__file__).resolve().parent
# pyHyp settings, tried in this order. Volume smoothing can leave cells with a negative
# corner at a leaning trailing-edge base; the solver accepts them, but a clean grid is preferred.
HYP = ['{"cMax":1.0,"volSmoothIter":30}', '{"cMax":1.0,"volSmoothIter":0}', '{"cMax":1.0,"volSmoothIter":30,"epsI":4.0}',
       '{"cMax":1.0,"volSmoothIter":30,"marchDist":0.03}']  # last: a deeper near grid, which removes orphan cells for some shapes
C0 = .5  # reference chord: C_L c = C0 * cl
LOSS_PER_CLC = .1  # merit of a run that minimizes the outlet loss under --max-clc: loss + this x its C_L c shortfall
CP_PER_CLC = 10.   # merit of a run that minimizes the suction peak (--min-suction): -Cp_min + this x its C_L c shortfall
# numpy's BLAS would otherwise start a thread per core, on cores the solver is using
ONE_THREAD = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')

def unresolved(p):
    """Absolute path of p with links kept: a run folder may link to a disk whose path has a space,
    which the commands (split on spaces, run by a shell) and the process matching must not see.
    Relative paths go against the shell's $PWD, since os.getcwd() resolves links."""
    pwd = os.environ.get('PWD')
    base = pwd if pwd and os.path.isdir(pwd) and os.path.samefile(pwd, '.') else os.getcwd()
    return Path(os.path.normpath(os.path.join(base, p)))

def log(*words):
    print(time.strftime('%m-%d %H:%M:%S'), *words, flush=True)

def history(run):
    f = run/'history.jsonl'
    return [json.loads(l) for l in f.read_text().splitlines() if l.strip()] if f.exists() else []

def pids(pattern):
    return [int(v) for v in subprocess.run(['pgrep', '-f', pattern], capture_output=True, text=True).stdout.split()]
def ranks(run): return pids(f'[p]ython -u {run}/source/')
def launcher(run): return pids(f'^(mpirun|prterun) .*--root {run} ')

def objective(run):
    """(score of a result record, its name, the option of its least gain): what progress is judged on, lower
    is better. C_L c; or, for a run that minimizes the outlet Cpt loss under --max-clc, the merit loss +
    LOSS_PER_CLC x how far C_L c is short of (above) its limit. Near these designs a unit of C_L c costs some
    0.01-0.05 of loss, so the merit is lowest at the constrained optimum, and an iterate still on its way to
    the limit counts as progress."""
    cmd = (run/'command.txt').read_text().split('\n')[0].split()
    t = option(cmd, '--max-clc')
    if t is None: return (lambda r: C0*r['cl']), 'C_L c', 'min_gain'
    t = float(t)
    if '--min-suction' in cmd:  # -Cp_min scales with the loading (some 2-4 per unit of C_L c), so 10 x the shortfall dominates
        return (lambda r: -r['cp_min_ks']+CP_PER_CLC*max(0., C0*r['cl']-t)), f'the suction merit (C_L c <= {t:g})', 'min_cp_gain'
    return (lambda r: r['cpt_loss_outlet']+LOSS_PER_CLC*max(0., C0*r['cl']-t)), f'the loss merit (C_L c <= {t:g})', 'min_loss_gain'

NEAR_CAP = 1.01  # an iterate counts for progress only within 1 % of every cap the run imposes

def caps(run):
    """The caps the run imposes, each as a function of a result record scaled as run_slip_optimization.py
    scales it (<= 1 within it): drag (--drag-cap caps cd; 1 or more is off), --loss-cap, --cp-min."""
    cmd = (run/'command.txt').read_text().split('\n')[0].split()
    out = []
    d = float(option(cmd, '--drag-cap') or .05)
    if d < 1 and option(cmd, '--loss-cap') is None and option(cmd, '--max-clc') is None: out.append(lambda r: r['cd']/d)
    if option(cmd, '--loss-cap') is not None: out.append(lambda r, c=float(option(cmd, '--loss-cap')): r['cpt_loss_outlet']/c)
    if option(cmd, '--cp-min') is not None: out.append(lambda r, c=float(option(cmd, '--cp-min')): r['cp_min_ks']/c)
    return out

def iterates(run, near=True):
    """(elapsed s, score) of the saved PSQP iterates. Progress is judged on these, not on strictly
    feasible evaluations: an iterate can sit microns outside a limit until PSQP converges. With `near`, only
    those within NEAR_CAP of every cap: from a start well past a cap (run_041's rotated start had 76 % too much
    drag) PSQP first trades the objective for feasibility, which is no stall (the lift limit of a --max-clc run
    is in its merit instead)."""
    score, _, least = objective(run); out = []
    within = caps(run) if near else []
    for d in (run/'iterates').glob('eval_*'):
        if (d/'design.json').exists():
            r = json.loads((d/'result.json').read_text())
            # a loss run skips the round's start on its unwarped grid: any warp adds 1-4 % to the loss
            if least == 'min_loss_gain' and r.get('zero_deformation', unwarped(r['design'])): continue
            if any(f(r) > NEAR_CAP for f in within): continue
            out.append((r['elapsed_wall_time_s'], score(r)))
    return sorted(out)

def unwarped(design):
    return design['incidence'][0] == -3. and not any(any(v) for k, v in design.items() if k != 'incidence')

def verdict(run, rows, a):
    """Why the grid needs replacing, from the evaluations so far; None while the run is healthy."""
    failed = [r['cfd_failed'] for r in rows]
    if len(rows) >= a.max_failures and all(failed[-a.max_failures:]):
        return f'the last {a.max_failures} evaluations failed'
    if len(rows) >= 8 and sum(failed[-8:]) >= 6:
        return '6 or more of the last 8 evaluations failed'
    its = iterates(run)
    every = iterates(run, near=False)
    if not its and every and every[-1][0]-every[0][0] > 4*60*a.stall_minutes:
        its = every  # still not near the caps after four stall windows: judge on every iterate, as before
    if its:
        now = rows[-1]['elapsed_wall_time_s']
        _, name, least = objective(run); least = getattr(a, least)
        before = [s for t, s in its if t < now-60*a.stall_minutes]
        # Judged once PSQP has stepped: the first two iterates are the derivative check and PSQP's
        # start, which a slow machine can push 15 min apart from the first step.
        if len(its) > 2 and before and min(before)-min(s for _, s in its) < least:
            return f'the PSQP iterates improved {name} by under {least:g} in {a.stall_minutes:g} min'
    return None

def gain(run):
    its = [s for _, s in iterates(run)]
    return its[0]-min(its) if its else 0.

def last_iterate(run, back=0):
    """The newest saved PSQP iterate, or the one `back` before it (the optimizer writes design.json last)."""
    saved = sorted(d for d in (run/'iterates').glob('eval_*') if (d/'design.json').exists())
    d = saved[max(0, len(saved)-1-back)]
    return d, json.loads((d/'result.json').read_text())

def next_name(name):
    """r00 -> r01; r03b, a retry on grid 3, -> r04."""
    n = int(re.match(r'r(\d+)', name).group(1))+1
    return f'r{n:02d}'

def option(cmd, name):
    return cmd[cmd.index(name)+1] if name in cmd else None

def retry_name(name):
    """r01 -> r01b -> r01c: another grid for the same run."""
    m = re.fullmatch(r'(r\d+)([a-z]?)', name)
    return m.group(1)+(chr(ord(m.group(2))+1) if m.group(2) else 'b')

def prepare(run, where, name=None, skip=(), back=0):
    """Steps 2-3 under the folder `where`, as run `name` (default: the next one), at the iterate
    `back` before the newest, leaving out the pyHyp settings in `skip`. Returns (root, command,
    pyHyp setting), or None without a usable grid."""
    if all(h in skip for h in HYP): return None
    it, r = last_iterate(run, back)
    cmd = (run/'command.txt').read_text().split('\n')[0].split()
    name = name or next_name(run.name)
    mesh = where/'grids'/f"{name}_from_{run.name}_{it.name.replace('_', '')}"
    root = where/name
    mesh.mkdir(parents=True)
    shutil.copy(it/'design.json', mesh/'source_design.json')
    (mesh/'source.txt').write_text(f"{run.name} evaluation {r['evaluation']}\n")
    design = json.loads((it/'design.json').read_text())
    # the derivative check runs just off the new grid's zero point, a kink of the warp
    (mesh/'offset_design.json').write_text(json.dumps({k: [-3.02] if k == 'incidence' else [0.]*len(v) for k, v in design.items()}))
    old = option(cmd, '--baseline')
    tries, hyps = [], [h for h in HYP if h not in skip]
    for k, hyp in enumerate(hyps):
        try:  # pyHyp can hang on a shape: a grid taking over 10 min counts as failed
            made = subprocess.run([sys.executable, str(S/'remesh_design.py'), '--design', str(mesh/'source_design.json'), '--out', str(mesh/f'try{k}'),
                                   '--hyp-options', hyp] + (['--baseline', old] if old else []), capture_output=True, text=True, env=ONE_THREAD, timeout=600)
        except subprocess.TimeoutExpired as e:
            made = subprocess.CompletedProcess(e.cmd, 1, str(e.stdout or ''), 'remesh_design.py timed out after 10 min\n')
        (mesh/f'remesh_try{k}.log').write_text(made.stdout+made.stderr)
        failed = bool(made.returncode)  # remesh_design.py exits with an error on inverted cells or a lost lattice
        corners = None if failed else sum(b['nonpositive_corners'] for b in json.loads((mesh/f'try{k}/baseline.json').read_text())['blocks'])
        tries.append({'pyhyp': json.loads(hyp), 'failed': failed, 'negative_corners': corners})
        if corners == 0: break
    valid = [k for k, t in enumerate(tries) if not t['failed']]
    use = min(valid, key=lambda k: tries[k]['negative_corners'], default=None)  # a clean grid, else the fewest negative corners
    for k, t in enumerate(tries): t['used'] = k == use
    (mesh/'remesh_attempts.json').write_text(json.dumps(tries, indent=1)+'\n')
    if use is None:
        log('remesh failed with every pyHyp setting (see remesh_try*.log) -', mesh); return None
    for f in (mesh/f'try{use}').iterdir(): shutil.move(str(f), mesh/f.name)
    shutil.move(str(mesh/f'remesh_try{use}.log'), mesh/'remesh.log')
    for k in range(len(tries)): shutil.rmtree(mesh/f'try{k}', ignore_errors=True)
    if tries[use]['negative_corners']:
        log(f"warning: the new grid has {tries[use]['negative_corners']} cells with a negative corner; no setting gave a clean one")
    shutil.copytree(run/'source', root/'source', ignore=shutil.ignore_patterns('__pycache__', 'aerodynamic_derivative_check.json'))
    shutil.copy(S/'run_slip_optimization.py', root/'source/run_slip_optimization.py')  # so a fix made during a chain reaches its next run
    for f in ('airfoil_coordinates.csv', 'ground.cgns'): shutil.copy(mesh/f, root/'source'/f)
    shutil.copy(it/'flow_000_vol.cgns', root/'source/restart_vol.cgns')
    shutil.copy(run/'rankfile', root/'rankfile')
    # The same command on the new grid, without the previous run's start. --solver and
    # --thickness-measure are retired: their values (ank_deep, gap) are now built in.
    for flag, values in [('--baseline', 1), ('--check-design', 1), ('--restart-design', 1), ('--restart-flow', 1), ('--reuse-checks', 0),
                         ('--hotstart', 1), ('--solver', 1), ('--thickness-measure', 1)]:
        if flag in cmd: del cmd[cmd.index(flag):cmd.index(flag)+1+values]
    cmd[cmd.index('--root')+1] = str(root)
    cmd[[i for i, w in enumerate(cmd) if w.endswith('/source/run_slip_optimization.py')][0]] = str(root/'source/run_slip_optimization.py')
    if '--first-ncycles' not in cmd: cmd += ['--first-ncycles', '60000']  # the first solve is on a new grid
    cmd += ['--baseline', str(mesh/'baseline.json'), '--check-design', str(mesh/'offset_design.json')]
    (root/'command.txt').write_text(' '.join(cmd)+f"\n# grid: remesh of {run.name} evaluation {r['evaluation']} (C_L c {C0*r['cl']:.3f}); "
                                    "design variables restart at zero; initial flow source/restart_vol.cgns = that evaluation's "
                                    'converged flow; set up and started by run_chain.py\n')
    log(f"remeshed {run.name} at evaluation {r['evaluation']} (C_L c {C0*r['cl']:.3f}) with pyHyp {json.dumps(tries[use]['pyhyp'])} -> {root.name}")
    return root, ' '.join(cmd), hyps[use]

def orphans(run):
    """Orphan cells (no donor) left in the overset grid; None if the solver reports nothing in 5 min.
    The solver assembles the overlap in passes and its first count is often large: the grid is bad
    only when the passes end non-zero ('Orphans present in grid'), clean once a pass counts 0. The
    flow solve tolerates orphans, but the adjoint then fails."""
    t = time.time()
    while time.time()-t < 300:
        text = (run/'run.log').read_text(errors='ignore')
        found = re.findall(r'Total number of orphans:\s*(\d+)', text)
        if 'Orphans present in grid' in text: return int(found[-1]) if found else 1
        if found and found[-1] == '0': return 0
        time.sleep(5)
    return None

def stop(run, why):
    for pid in launcher(run): subprocess.run(['kill', str(pid)])
    t = time.time()
    while ranks(run):
        if time.time()-t > 60:
            for pid in ranks(run)+launcher(run): subprocess.run(['kill', '-9', str(pid)])
        time.sleep(2)
    with (run/'command.txt').open('a') as f:
        f.write(f"stopped by run_chain.py {time.strftime('%c')}: {why}\n")

def check_verdict(run):
    """The verdict of the run's derivative check; None until it is written."""
    try: return json.loads((run/'source/aerodynamic_derivative_check.json').read_text())['verdict']
    except (OSError, ValueError): return None

def watch(run, a):
    """Follow the run until its grid needs replacing. Returns (why, whether the chain can go on)."""
    warned = False
    while True:
        rows = history(run)
        v = check_verdict(run)
        if not warned and v not in (None, 'passed', 'failed'):
            warned = True; log(f'warning: {run.name}: derivative check {v}; the run goes on')
        if not ranks(run):
            if (run/'run_complete.json').exists():
                return 'PSQP exited by itself: '+str(json.loads((run/'run_complete.json').read_text())['inform']), True
            if 'stopped by' in (run/'command.txt').read_text():
                return 'it was already stopped', True
            return 'the run died by itself; see the end of its run.log', False
        why = verdict(run, rows, a)
        newest = max(f.stat().st_mtime for f in (run/'history.jsonl', run/'command.txt') if f.exists())
        if not why and time.time()-newest > 60*a.hang_minutes:
            why = f'no evaluation finished in {a.hang_minutes:g} min'
        if why:
            stop(run, why)
            return why, True
        time.sleep(a.poll)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('runs', type=Path, nargs='+')
    p.add_argument('--live-plots', type=Path, metavar='OUT', help='start live_plots.py on every new run, with the chain plots in OUT')
    p.add_argument('--max-failures', type=int, default=3, help='failed evaluations in a row that mean a worn-out grid')
    p.add_argument('--stall-minutes', type=float, default=30.)
    p.add_argument('--min-gain', type=float, default=.005, help='in C_L c')
    p.add_argument('--min-loss-gain', type=float, default=3e-5, help='in outlet Cpt loss (the merit of a run with --max-clc)')
    p.add_argument('--min-cp-gain', type=float, default=.01, help='in -Cp_min (the merit of a run with --max-clc --min-suction)')
    p.add_argument('--hang-minutes', type=float, default=60.)
    p.add_argument('--max-runs', type=int, default=12, help='runs this driver may start')
    p.add_argument('--min-free-gb', type=float, default=3., help='free disk below which no new run starts (a run writes up to 1 GB)')
    p.add_argument('--poll', type=float, default=30., help='s between looks at the run')
    p.add_argument('--dry-run', type=Path, metavar='DIR', help='only set up the next grid and run folder of the last run, under DIR')
    a = p.parse_args()
    runs = [unresolved(r) for r in a.runs]
    if a.dry_run:
        print(prepare(runs[-1], unresolved(a.dry_run))); return
    log('chain driver started on', runs[-1].name, '(running)' if ranks(runs[-1]) else '(stopped)')
    for started in range(a.max_runs+1):
        run = runs[-1]
        why, go_on = watch(run, a)
        its = [s for _, s in iterates(run)]
        _, name, least = objective(run); n = 6 if least == 'min_loss_gain' else 3; least = getattr(a, least)
        log(f"{run.name}: {why}. Its PSQP iterates went from {name} {its[0]:.{n}f} to {min(its):.{n}f} "
            f"in {len(history(run))} evaluations." if its else f'{run.name}: {why}. No iterate was saved.')
        if not go_on:
            log('chain ended: the run needs a look.'); return
        if gain(run) < least:
            log(f'chain ended: this run gained under {least:g}, so another grid would not help.'); return
        before = [s for r in runs[:-1] for _, s in iterates(r)]
        if before and min(before)-min(its) < least:
            log(f'chain ended: this run did not beat the best of the earlier runs by {least:g}.'); return
        if started == a.max_runs:
            log(f'chain ended: --max-runs {a.max_runs} runs started.'); return
        free = shutil.disk_usage(run).free/1e9
        if free < a.min_free_gb:
            log(f'chain ended: {free:.1f} GB of disk left, under --min-free-gb {a.min_free_gb:g}.'); return
        skip, name, back = [], None, 0
        while True:
            new = prepare(run, run.parent, name, skip, back)
            if new is None and back < 3:  # no usable grid at this iterate: try the one before
                back, skip = back+1, []
                log(f'no usable grid at that iterate; trying the iterate {back} before the newest'); continue
            if new is None:
                log('chain ended: no new grid.'); return
            root, cmd, hyp = new
            subprocess.Popen(cmd, shell=True, cwd=root, stdout=(root/'run.log').open('w'), stderr=subprocess.STDOUT, start_new_session=True)
            n = orphans(root)
            if not n: break
            stop(root, f'{n} orphan cells in the overset grid; the adjoint cannot converge with them')
            log(f'{root.name}: {n} orphan cells in the overset grid; trying the next pyHyp setting')
            skip.append(hyp); name = retry_name(root.name)
        runs.append(root); log('started', root.name)
        if a.live_plots:
            out = unresolved(a.live_plots); out.mkdir(parents=True, exist_ok=True)
            subprocess.Popen([sys.executable, '-u', str(S/'live_plots.py'), '--out', str(out), *map(str, runs)],
                             stdout=(out/'live_plots.log').open('a'), stderr=subprocess.STDOUT, start_new_session=True)
        time.sleep(60)  # the ranks are up well before the first look

if __name__ == '__main__':
    main()

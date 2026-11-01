"""Chained PSQP continuation until the surface reaches the box x limits.

Each stage is a fresh run_slip_optimization.py run restarted from the previous
stage's best design and a converged flow near it. Between stages:
  - a drag constraint at its cap doubles the cap. The cap counts as hit once the
    last two geometrically feasible evaluations both have C_D >= DRAG_ACTIVE*cap
    (PSQP rides the cap from either side, often slightly above it); the stage is
    then ended early;
  - the next stage starts from the evaluation with the most downforce that is
    feasible under the next stage's cap, with the stage's best_feasible flow as
    the initial flow (the only saved flow; it is a nearby design);
  - the chain stops when both x limits are active, or a stage gains less than
    STALL_GAIN without reaching the drag cap.
Run outside the container with the MACH-Aero environment active:
  python scripts/run_continuation.py --chain optimizations/<name> \
      --start-design D.json --start-flow F.cgns --xmax 2.0 --ranks 24 --first-core 8
Ranks are pinned to cores first-core.. via an OpenMPI rankfile, so concurrent
chains do not share cores. Rerunning the same command resumes the chain; a stage
left without stage_summary.json (driver killed) is summarised from its history.
"""
import argparse
import numpy as np
import json
import shutil
import subprocess
import time
from pathlib import Path

P = Path(__file__).resolve().parents[1]
BOX = {'xmin': -.125, 'xmax': .625, 'ymin': .05, 'ymax': .35, 'tmin': .03}
p = argparse.ArgumentParser()
p.add_argument('--chain', type=Path, required=True)
p.add_argument('--start-design', type=Path, required=True)
p.add_argument('--start-flow', type=Path, required=True,
               help="converged flow at --start-design, or 'none' for a remeshed baseline: stage 00 then "
                    "solves it from freestream and runs its own aerodynamic derivative check")
p.add_argument('--source-dir', type=Path, default=P/'setups/slip_ground_xy5_cd005/source',
               help='optimizer inputs copied into every stage (grid, coordinates, ...)')
p.add_argument('--baseline', type=Path, help='remesh_design.py baseline.json, passed to the optimizer')
p.add_argument('--checks', type=Path, help='passed aerodynamic_derivative_check.json to reuse; with '
               "--start-flow none it is the remeshed grid's own check (taken off the zero-deformation point, "
               'where the IDWarp response has a kink), and stage 00 then skips the baseline check')
p.add_argument('--xmax', type=float, default=2.)
p.add_argument('--maxit', type=int, default=15)
p.add_argument('--drag-cap', type=float, default=.05, help='C_D cap on the 0.5 m reference chord')
p.add_argument('--min-local-thickness', type=float, default=0., help='m; passed to the optimizer (0 = off)')
p.add_argument('--save-iterate-flows', action='store_true', help='passed to the optimizer')
p.add_argument('--nk-switch-tol', type=float, help='passed to the optimizer')
p.add_argument('--l2-convergence', type=float, help='passed to the optimizer')
p.add_argument('--adjoint-l2-convergence', type=float, help='passed to the optimizer')
p.add_argument('--nk-fallback-tol', type=float, help='passed to the optimizer')
p.add_argument('--warp-err-tol', type=float, help='passed to the optimizer')
p.add_argument('--solver', default='ank_deep')
p.add_argument('--ncycles', type=int, default=20000)
p.add_argument('--first-stage-ncycles', type=int,
               help='nCycles for stage 00 of a remeshed chain, which starts from freestream; later stages '
                    'restart from converged flows, where a smaller cap catches a stalled solve sooner')
p.add_argument('--ranks', type=int, default=24)
p.add_argument('--first-core', type=int, default=8)
p.add_argument('--stage-timeout', default='6h')
p.add_argument('--max-stages', type=int, default=12)
p.add_argument('--x-tol', type=float, default=1e-3, help='m; x limit counts as active within this')
p.add_argument('--drag-active', type=float, default=.995)
p.add_argument('--stall-gain', type=float, default=.002)
a = p.parse_args()
reuse_without_flow = a.checks is not None and str(a.start_flow) == 'none'
if a.checks is None:
    a.checks = P/'optimizations/xmax_sweep_20261031/failed_ncycles4000_xmax_0.05/source/aerodynamic_derivative_check.json'
chain = a.chain.resolve(); chain.mkdir(parents=True, exist_ok=True)
log = chain/'chain_log.jsonl'
TOL = 5e-7  # m; the optimizer's feasibility tolerance (1e-6 of the 0.5 m reference chord)
# PSQP's linearised local-thickness floor is overshot by ~1-2 um on large steps; a restart
# design may undercut the floor by this much (0.25 % of 2 mm) instead of losing the step.
LOCAL_TOL = 5e-6

def records(stage):
    f = stage/'history.jsonl'
    return [json.loads(l) for l in f.read_text().splitlines()] if f.exists() else []

def geometric_ok(r):
    return (not r['cfd_failed'] and r['minimum_x_m'] >= BOX['xmin']-TOL and r['maximum_x_m'] <= BOX['xmax']+TOL
            and r['minimum_clearance_m'] >= BOX['ymin']-TOL and r['maximum_height_m'] <= BOX['ymax']+TOL
            and r['maximum_section_thickness_m'] >= BOX['tmin']-TOL
            and r.get('minimum_local_thickness_m', np.inf) >= a.min_local_thickness-LOCAL_TOL)

def drag_hit(h, cap):
    g = [r for r in h[1:] if geometric_ok(r)]
    return len(g) >= 2 and all(r['cd']/cap >= a.drag_active for r in g[-2:])

def summarise(stage, k, cap, rc, wall, drag_stop):
    h = records(stage)
    if not h or not (stage/'best_feasible/flow_000_vol.cgns').exists():
        return {'stage': k, 'exit': rc, 'stop': 'no feasible design saved'}
    drag_active = drag_stop or drag_hit(h, cap)
    nextcap = cap*2 if drag_active else cap
    cands = [r for r in h if geometric_ok(r) and r['cd'] <= nextcap*(1+1e-5)]
    best = min(cands, key=lambda r: r['cl'])
    (stage/'next_design.json').write_text(json.dumps(best['design'], indent=2)+'\n')
    gain = (best['cl']-h[0]['cl'])/h[0]['cl']
    lo, hi = best['minimum_x_m']-BOX['xmin'], BOX['xmax']-best['maximum_x_m']
    s = {'stage': k, 'exit': rc, 'drag_stop': drag_stop, 'wall_min': wall/60, 'xmax': a.xmax, 'ranks': a.ranks,
         'drag_cap': cap, 'evaluations': len(h), 'failed_evaluations': sum(r['cfd_failed'] for r in h),
         'start_clc_m': h[0]['cl']*.5, 'next_start_evaluation': best['evaluation'],
         'best_clc_m': best['cl']*.5, 'best_cdc_m': best['cd']*.5, 'gain': gain,
         'drag_ratio': best['cd']/cap, 'drag_active': drag_active, 'next_drag_cap': nextcap,
         'x_margin_lead_m': lo, 'x_margin_trail_m': hi, 'chord_m': best['actual_chord_m'],
         'min_clearance_m': best['minimum_clearance_m'],
         'max_thickness_m': best['maximum_section_thickness_m'], 'incidence_deg': best['chord_incidence_deg']}
    nact = (lo <= a.x_tol) + (hi <= a.x_tol)
    if nact == 2: s['stop'] = 'both x limits active'
    elif gain < a.stall_gain and not drag_active: s['stop'] = 'stage stalled before the x limits'
    return s

def watch(stage, cap, alive, stop):
    """Poll a running stage; end it on the drag cap. Returns whether it did."""
    while alive():
        time.sleep(20)
        # End the stage on the drag cap only after that evaluation's gradient is
        # written: best_feasible is saved by then, and the kill lands in a flow
        # solve rather than in a file write.
        h = records(stage)
        g = (stage/'gradients.jsonl').read_text().splitlines() if (stage/'gradients.jsonl').exists() else []
        if drag_hit(h, cap) and g and json.loads(g[-1])['evaluation_count'] == h[-1]['evaluation']+1:
            stop(); return True
    return False

def stage_pids(stage):
    r = subprocess.run(['pgrep', '-f', f'timeout .*{stage}/source/run_slip'], capture_output=True, text=True)
    return [int(x) for x in r.stdout.split()]

def restart_flow(stage, s):
    """The next stage's initial flow: the chosen iterate's own saved flow (--save-iterate-flows),
    else the stage's best feasible flow, which may be several steps away."""
    f = stage/f"iterates/eval_{s.get('next_start_evaluation', -1):03d}/flow_000_vol.cgns"
    return f if f.exists() else stage/'best_feasible/flow_000_vol.cgns'

def finish(stage, s):
    (stage/'stage_summary.json').write_text(json.dumps(s, indent=2)+'\n')
    with log.open('a') as f: f.write(json.dumps(s)+'\n')
    print(json.dumps(s), flush=True)
    if (stage/'start_vol.cgns').exists(): (stage/'start_vol.cgns').unlink()  # the previous stage's flow

design, cap = a.start_design.resolve(), a.drag_cap
flow = None if str(a.start_flow) == 'none' else a.start_flow.resolve()
done = sorted(chain.glob('stage_*'))
if str(a.start_flow) == 'none' and (chain/'stage_00/source/aerodynamic_derivative_check.json').exists():
    a.checks = chain/'stage_00/source/aerodynamic_derivative_check.json'  # the remeshed baseline's own check
if done:
    last = done[-1]
    if not (last/'stage_summary.json').exists():
        prev = float((last/'command.txt').read_text().split('--drag-cap ')[1].split()[0])
        pids, tic = stage_pids(last), time.time()
        # A stage still running under a previous driver is adopted, not restarted.
        def alive(): return bool(stage_pids(last))
        def stop():
            for q in pids: subprocess.run(['kill', str(q)])
            while alive(): time.sleep(2)
        stopped = watch(last, prev, alive, stop) if pids else False
        finish(last, summarise(last, len(done)-1, prev, None, time.time()-tic, stopped))
    s = json.loads((last/'stage_summary.json').read_text())
    if 'stop' in s: raise SystemExit(f'chain stopped at {last.name}: {s["stop"]}')
    design, flow, cap = last/'next_design.json', restart_flow(last, s), s['next_drag_cap']
for k in range(len(done), a.max_stages):
    stage = chain/f'stage_{k:02d}'
    if str(a.start_flow) == 'none' and (chain/'stage_00/source/aerodynamic_derivative_check.json').exists():
        a.checks = chain/'stage_00/source/aerodynamic_derivative_check.json'  # the remeshed baseline's own check
    (stage/'source').mkdir(parents=True)
    for f in a.source_dir.resolve().iterdir():
        shutil.copy(f, stage/'source')
    shutil.copy(P/'scripts/run_slip_optimization.py', stage/'source')
    if flow is not None or reuse_without_flow:
        shutil.copy(a.checks, stage/'source/aerodynamic_derivative_check.json')
    if flow is not None:
        shutil.copy(flow, stage/'start_vol.cgns')
    shutil.copy(design, stage/'start_design.json')
    (stage/'rankfile').write_text(''.join(f'rank {i}=localhost slot={a.first_core+i}\n' for i in range(a.ranks)))
    cmd = ['timeout', a.stage_timeout, 'mpirun', '-np', str(a.ranks), '--map-by', 'rankfile:file=rankfile',
           '--bind-to', 'core', '--mca', 'btl', 'self,sm', '-x', 'GFORTRAN_UNBUFFERED_PRECONNECTED=y',
           'python', '-u', str(stage/'source/run_slip_optimization.py'), '--root', str(stage),
           '--ffd-stations', '5', '--shape-directions', 'xy', '--drag-cap', str(cap), '--min-local-thickness', str(a.min_local_thickness), '--optimizer', 'PSQP',
           '--maxit', str(a.maxit), '--xmax', str(a.xmax),
           '--ncycles', str(a.first_stage_ncycles if flow is None and a.first_stage_ncycles else a.ncycles), '--solver', a.solver,
           '--restart-design', str(stage/'start_design.json')] \
        + (['--reuse-checks', '--restart-flow', str(stage/'start_vol.cgns')] if flow is not None
           else ['--reuse-checks'] if reuse_without_flow else []) \
        + (['--baseline', str(a.baseline.resolve())] if a.baseline else []) \
        + (['--save-iterate-flows'] if a.save_iterate_flows else []) \
        + [x for f in ['nk_switch_tol', 'l2_convergence', 'adjoint_l2_convergence', 'nk_fallback_tol', 'warp_err_tol'] if getattr(a, f) is not None
           for x in ('--'+f.replace('_', '-'), str(getattr(a, f)))]
    (stage/'command.txt').write_text(' '.join(cmd)+'\n')
    tic = time.time()
    with (stage/'run.log').open('w') as out:
        proc = subprocess.Popen(cmd, cwd=stage, stdout=out, stderr=subprocess.STDOUT)
        drag_stop = watch(stage, cap, lambda: proc.poll() is None, lambda: (proc.terminate(), proc.wait()))
        rc = proc.returncode
    s = summarise(stage, k, cap, rc, time.time()-tic, drag_stop)
    finish(stage, s)
    if 'stop' in s: break
    design, flow, cap = stage/'next_design.json', restart_flow(stage, s), s['next_drag_cap']

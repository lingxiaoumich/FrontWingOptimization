"""Run the min-loss Cp-limited runs 036-039 until PSQP's own KKT test is met (user, 2026-11-25: "run 36-39 does not
look to have converged. please try to keep running them until KKT is met").

Every round of these chains was ended by run_chain.py's stall rule (under 3e-5 of loss merit in 15 min), never by
PSQP. Here each run gets new rounds without the stall rule: a round is a remesh at the newest PSQP iterate
(run_chain.prepare, as any chain round: fresh grid at that design, its flow as the initial flow, design variables at
zero) and runs until PSQP exits by itself. PSQP's exit 4 is its KKT test (GMAX, the largest component of the
Lagrangian gradient, <= TOLG 1e-4, with the constraints met to TOLC 1e-6; pyOptSparse's message for 4, "Maximum
constraint value ...", is misleading): the run is done. Any other exit, or a worn-out grid (3 failed evaluations in a
row, 6 of 8, or an hour without an evaluation), starts the next round from the newest iterate.
run_038 (C_L c <= -1.4, Cp >= -4) never had a strictly feasible design: its loss-cap starts sit at Cp -5.2..-5.7. Its
continuation starts instead from run_042_area050MinSuctionClc12CapStart r00 evaluation 23 (C_L c -1.4006, Cp_min
-3.875, loss 0.005428), which is feasible for run_038's limits; run_038's own command (limits) is kept.
Four runs at once, 12 ranks each (cores 12-23, 24-35, 36-47, 48-59); helpers on cores 8-11 and 60-63.
Second start (22:3x): the first rounds (036/037 r02, 039 r03) all ended with PSQP exit -6 after one iteration and a 21-step
line search: the outlet loss jumped by up to 3 % between designs 1e-8 apart. The flow fields of those trial points
differ only by a uniform pressure level (ratio of loss jump to field Cp shift 0.87-0.99), which the low-Mach solve leaves
unconverged by ~1e-6 p while C_L and Cp are immune. The rounds from then on measure the loss against the same solution's
total pressure on an upstream plane, run_slip_optimization.py --loss-reference inlet (x = -4 m), which cancels it
(emulated on the saved fields: jumps 7e-5 -> ~2e-6). The rounds running at the switch were stopped.
Fourth start (11-26 01:3x): the loss planes' points moved off the background cell-centre rows (--plane-points offset) and
each round starts at shape variables of 1e-9 m instead of the grid as made; see RUNS.
11-26 09:xx: restarted (following all four rounds) for the 48 h continuation (user: "keep trying for the next 48 hours or until kkt
is met"), with WATCH.hang_minutes 120.
Fifth start (11-26 14:48): 036 met KKT in r06 (13:59), but 037 and 039 sat with every line search failing: each accepted
point lay in a small pocket where the loss read 5e-7 (037) lower than at any trial 1e-5 m or more away, reproducibly and
whatever the solve's start. The cause is the far-field BC on the top boundary: IDWarp warps the whole overset system with
only the airfoil as a moving surface, so the domain boundary moves with the wing (36-56 mm by 037 r07 eval 172) and the
top picks up local slopes of +-7e-3; the far-field BC treats each face as inflow or outflow by the sign of its normal
velocity (U times that slope), so the converged flow jumps wherever a face changes over. A replay with the far field held
in the warp (run_slip_optimization.py --fixed-farfield) and the top tilted to outflow everywhere (--top-tilt 2e-3) is
smooth to 1e-9 through the same steps. That also changes the model away from a round's start (037 r07 eval 172: loss
-1.3 %, C_L +0.25 %), and removes the far-field drift that made remeshes move C_L c by 1-3 % (038). All four runs
continue on the new model: 036 from its KKT point, 037/039 from their last accepted points, 038 from r07 eval 49.
Sixth start (11-26 19:4x): 039 met KKT in r08 (19:32). 038 r08 cycled for 13 iterations (C 0.06-0.2, F 1.47-1.57, every
step accepted at full length): each step's linearization predicted Cp_min to fall 0.06-0.16 onto its limit, and the
solved Cp_min fell by about nothing, the same sign every step (039's early steps: errors <1e-3). The suction peak is that
curved (derivative check along its gradient: slope 49.8 forward vs 41.5 backward at 1e-4 m), and 1.2-4.3 mm steps
(XMAX 1.0 = 10 mm) leave a curvature error that the next step only removes by making another. 038 now runs with
--xmax 0.05 (0.5 mm, the driver's default; error ~(0.5/2.5)^2 x 0.1) on 24 ranks (cores 36-59, 039's freed), from r08's
newest iterate; 036/037 keep their rounds (queue restarted with arguments 036 037 038).
11-27 00:2x: 036 (r07b, 21:57) and 037 (r08, 00:09) met KKT too. 038 r09 is past its cycling (C 0.05 -> 5e-4) and runs
in a straight line at the 0.5 mm cap (cos 1.000 between steps, linearization errors ~1e-7), some 56 mm from its next
constraint (section area). A hedge on the freed cores 12-35: the branch 'run_038_...@branch_xmax020' (rounds under
run_038_.../branch_xmax020/, template from the main chain) from r09's newest iterate with --xmax 0.2; a second queue
process runs it alone (python -u kkt_queue.py xmax020). The main lane is untouched.
11-27 09:2x: the branch (F 0.907 at r10 NIT 198) left the main lane behind (F 1.199, crawling along the section-area
bound at the 0.5 mm cap), and its 2 mm steps have stayed straight with linearization errors ~9e-7 for 100+ steps, the
objective slope falling only 4.6e-4 -> 3.7e-4 per mm in 108 mm: travel-limited. The main lane (r09) is stopped and its
cores 36-59 run 'run_038_...@branch_xmax050' from the branch's newest iterate with --xmax 0.5 (python -u kkt_queue.py
xmax050); the 2 mm branch goes on as the safe lane.
11-27 15:2x: both branches reached the solver's noise floor: G 2-4e-4, C ~5e-8, line searches halving to <1e-5 mm on F
noise of 1e-8..6e-7, and branch_xmax020 r10 ended with PSQP -6. A replay from its newest iterate (scratchpad modetest
N12/N13, 0.36 mm toward eval 385) reads F to +-3e-7 with --l2-convergence 1e-12 and +-1e-8 with 1e-13 (solves 10-25 %
longer). Not through this queue (a round here always remeshes, which moves Cp_min ~0.006 and the KKT point with it):
same-grid continuations launched by hand, each round's command with --l2-convergence 1e-13, --reuse-checks and the
newest iterate's design and flow as restart: branch_xmax020/r10t (r10 eval 432, cores 12-35), branch_xmax050/r11t (r11
eval 226, cores 36-59). See each command.txt.
Run with the MACH-Aero env, detached:
  setsid nohup taskset -c 8-11,60-63 python -u kkt_queue.py > kkt.log 2>&1 &
"""
import json
import shutil
import sys
import threading
import time
import types
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # FQ.L lanes.py, FQ.RC run_chain, FQ.PY, FQ.S, FQ.subprocess
L, RC = FQ.L, FQ.RC
O = L.O
# run: (round whose command is the template, cores, start (round folder, evaluation or None for its newest iterate), first new round)
# First start: 036/037 from r01's newest, 038 from (run_042 r00, 23), 039 from r02's newest. Second start (inlet reference):
# Second start (22:17, inlet reference): 036/037 from r02's newest -> r04, 038 from r02's newest -> r03, 039 from r03's -> r05.
# Third start (23:0x): with the inlet reference the loss still jumped by 8.5e-6 between PSQP's base point and the limit of
# its line search (036/037 r04): the pressure-level offset is not quite uniform (downstream it moves 10-15 % more than
# upstream). Solves stop at ~4e-11 of the freestream residual, a 1e-6..1e-7 reduction from the warm start, which leaves
# the slow pressure-level mode where the previous design had it; Newton-Krylov gains a factor 5-7 an iteration with no
# plateau, so --l2-convergence 1e-12 (100x tighter) costs ~3 more NK iterations a solve. 038's r03 first solve (from
# r02 eval 25) cycled between 1e-6 and 6e-5 for 29,000 iterations; it restarts from r02 eval 24, which solved normally.
# Fourth start (01:3x, 11-26): the loss "noise" was the loss planes themselves. Every other y point of the outlet and inlet
# planes (1.25 mm points over 2.5 mm cells) lies on a row of background cell centres, the faces of ADflow's donor cells, so
# points there are lost or kept by round-off and by the rank split: run_039 r06 eval 38's saved flow reads loss 0.00420 on
# 12 ranks and 0.00435 on 4 (outlet mass flow 16.706, 16.745 kg/s); with the points moved to the midpoints
# (run_slip_optimization.py --plane-points offset) both read 0.004764 and 16.969 kg/s, and a replayed line search is smooth
# down to steps of 2.4e-4 of its first trial step. Separately, the grid as made (all design variables exactly zero, each
# round's start) is not the limit of the warped grids: C_L c 2.4e-5 and the loss 1.1e-5 off, reproducibly. So each round now
# starts PSQP at shape variables of 1e-9 m (START), not at zero. Each run continues from its newest iterate.
RUNS = {
    # fifth start (see the docstring): every run on --fixed-farfield --top-tilt 2e-3, from these points
    # 036 r07's first solve stalled in NK at ~2e-9 of the freestream residual (steps 0.03-0.18, 50 min), as 038's did with
    # --nk-switch-tol 1e-6: r07b, the same start, without it (DROP)
    'run_036_area050MinLossClc10CpMin4CapStart': ('r06', range(12, 24), (O/'run_036_area050MinLossClc10CpMin4CapStart/r06', 215), 'r07b'),  # its r06 KKT point
    'run_037_area050MinLossClc08CpMin4CapStart': ('r07', range(24, 36), (O/'run_037_area050MinLossClc08CpMin4CapStart/r07', 222), 'r08'),
    # 038 r07 never became feasible (C_L c -1.36..-1.40, Cp -3.9..-5.2) in 50 evaluations of 10-50 min; eval 49 is its
    # nearest to the C_L limit (-1.3978, Cp -4.21)
    # 'run_038_area050MinLossClc14CpMin4CapStart': ('r07', range(36, 48), (O/'run_038_area050MinLossClc14CpMin4CapStart/r07', 49), 'r08'),
    # sixth start: from r08's newest iterate, with --xmax 0.05 (SET), on 24 ranks
    'run_038_area050MinLossClc14CpMin4CapStart': ('r08', range(36, 60), (O/'run_038_area050MinLossClc14CpMin4CapStart/r08', None), 'r09'),
    'run_039_area050MinLossClc12CpMin4': ('r07', range(48, 60), (O/'run_039_area050MinLossClc12CpMin4/r07', 188), 'r08'),
    # 'run@branch': a lane whose rounds go under O/run/branch, its template from O/run/<round>
    'run_038_area050MinLossClc14CpMin4CapStart@branch_xmax020': ('r09', range(12, 36), (O/'run_038_area050MinLossClc14CpMin4CapStart/r09', None), 'r10'),
    'run_038_area050MinLossClc14CpMin4CapStart@branch_xmax050': ('branch_xmax020/r10', range(36, 60), (O/'run_038_area050MinLossClc14CpMin4CapStart/branch_xmax020/r10', None), 'r11'),
}
EXTRA = ['--loss-reference', 'inlet', '--l2-convergence', '1e-12', '--plane-points', 'offset', '--top-tilt', '2e-3']  # (flag, value) pairs added to every round's command
FLAGS = ['--fixed-farfield']  # flags without a value added to every round's command
START = 1e-8  # m: a round's first design has shape variables +-START, alternating (see the fourth start). Not all +START:
# that moves every FFD point alike, a rigid 1e-9 m translation, which is the grid as made again (037 r06's first line search
# converged to a loss 1.6e-5 above its base, as the zero design's does)
DROP = {'run_038_area050MinLossClc14CpMin4CapStart': ['--nk-switch-tol'], 'run_036_area050MinLossClc10CpMin4CapStart': ['--nk-switch-tol']}  # flags (with their value) a run's rounds go without
SET = {'run_038_area050MinLossClc14CpMin4CapStart': ['--xmax', '0.05'], 'run_038_area050MinLossClc14CpMin4CapStart@branch_xmax020': ['--xmax', '0.2'],
       'run_038_area050MinLossClc14CpMin4CapStart@branch_xmax050': ['--xmax', '0.5']}  # (flag, value) pairs a run's rounds get, overriding the template (sixth start)
NP = 12  # ranks of a run unless its cores say otherwise (command(..., np=len(cores)))
MAX_ROUNDS = 12
WATCH = types.SimpleNamespace(max_failures=3, stall_minutes=1e9, min_gain=.005, min_loss_gain=3e-5, min_cp_gain=.01,
                              hang_minutes=120., poll=30.)  # 120, not 60: run_038's evaluations take up to 54 min (walks of 3 solves) (11-26)
lock = threading.Lock()
def log(*w):
    with lock: L.log(*w)

def command(template, root, mesh, np=NP):
    """The round's command: `template` (a command.txt first line) on `root` and the grid `mesh`, as run_chain.prepare
    rewrites it, with NP ranks."""
    cmd = template.split()
    for flag, values in [('--baseline', 1), ('--check-design', 1), ('--restart-design', 1), ('--restart-flow', 1), ('--reuse-checks', 0),
                         ('--hotstart', 1), ('--solver', 1), ('--thickness-measure', 1)]:
        if flag in cmd: del cmd[cmd.index(flag):cmd.index(flag)+1+values]
    cmd[cmd.index('--root')+1] = str(root)
    cmd[[i for i, w in enumerate(cmd) if w.endswith('/source/run_slip_optimization.py')][0]] = str(root/'source/run_slip_optimization.py')
    cmd[cmd.index('-np')+1] = str(np)
    if '--first-ncycles' not in cmd: cmd += ['--first-ncycles', '60000']
    for k in range(0, len(EXTRA), 2):  # set, overriding a value the template has
        if EXTRA[k] in cmd: cmd[cmd.index(EXTRA[k])+1] = EXTRA[k+1]
        else: cmd += EXTRA[k:k+2]
    cmd += [f for f in FLAGS if f not in cmd]
    start = root/'source/start_design.json'  # written by new_round: zero would be the grid as made
    return ' '.join(cmd + ['--baseline', str(mesh/'baseline.json'), '--check-design', str(mesh/'offset_design.json'),
                           '--restart-design', str(start), '--restart-flow', str(root/'source/restart_vol.cgns')])

def back_of(rnd, evaluation):
    saved = sorted(d for d in (rnd/'iterates').glob('eval_*') if (d/'design.json').exists())
    return len(saved)-1-saved.index(rnd/'iterates'/f'eval_{evaluation:03d}')

def new_round(link, source, back, name, template, cores):
    """Remesh `source`'s iterate `back` before its newest as round `name` of `link`, start it; (root, why-not)."""
    skip = []
    while True:
        made = RC.prepare(source, link, name, skip, back)
        if made is None: return None, 'no pyHyp setting gave a usable grid'
        root, _, hyp = made
        mesh = next(p for p in (link/'grids').iterdir() if p.name.startswith(f'{name}_from_'))
        for f in (root/'source').glob('*.py'):
            if (L.S/f.name).exists(): shutil.copy(L.S/f.name, f)  # the current scripts
        design = json.loads((source/'iterates'/sorted(d.name for d in (source/'iterates').glob('eval_*') if (d/'design.json').exists())[0]/'design.json').read_text())
        (root/'source/start_design.json').write_text(json.dumps({k: ([START*(-1)**i for i in range(len(v))] if k != 'incidence' else [-3.]*len(v)) for k, v in design.items()}) + '\n')
        (root/'rankfile').write_text(''.join(f'rank {i}=localhost slot={c}\n' for i, c in enumerate(cores)))
        note = (root/'command.txt').read_text().split('\n', 1)[1]
        (root/'command.txt').write_text(command(template, root, mesh, len(cores)) + '\n' + note +
            f"# KKT continuation, comparisons/kkt_20261125/kkt_queue.py: no stall rule, PSQP runs until it exits; {len(cores)} ranks on cores "
            f"{cores[0]}-{cores[-1]}, {time.strftime('%c')}\n")
        FQ.subprocess.Popen((root/'command.txt').read_text().split('\n')[0], shell=True, cwd=root, stdout=(root/'run.log').open('w'),
                            stderr=FQ.subprocess.STDOUT, start_new_session=True)
        t = time.time()
        while not L.ranks(root) and time.time()-t < 600: time.sleep(10)
        if not L.ranks(root): return root, 'the ranks never came up'
        n = RC.orphans(root)
        if not n: return root, None
        RC.stop(root, f'{n} orphan cells in the overset grid; the adjoint cannot converge with them')
        log(f'{link.name}: {name}: {n} orphan cells; next pyHyp setting')
        skip.append(hyp); name = RC.retry_name(name)

def chain_rounds(link):
    """Rounds of the chain so far, in order (those with a history)."""
    rs = [r for r in sorted(link.glob('r[0-9]*')) if r.is_dir() and (r/'history.jsonl').exists() and (r/'history.jsonl').stat().st_size]
    return sorted(rs, key=lambda r: (int(r.name[1:3]), r.name))

def kkt_state(root):
    """PSQP's exit of a finished round: (value, message) or None."""
    f = root/'run_complete.json'
    if not f.exists(): return None
    inf = json.loads(f.read_text()).get('inform')
    return (inf.get('value'), inf.get('message')) if isinstance(inf, dict) else (inf, '')

def last_record(root):
    rows = [json.loads(l) for l in (root/'history.jsonl').read_text().splitlines() if l.strip()]
    ok = [r for r in rows if not r['cfd_failed']]
    return ok[-1] if ok else None

def drive(name):
    last, cores, start, rnd = RUNS[name]
    run, _, branch = name.partition('@')
    link = O/run/branch if branch else O/run
    link.mkdir(exist_ok=True)
    template = (O/run/last/'command.txt').read_text().split('\n')[0]
    for flag in DROP.get(name, []):
        w = template.split()
        if flag in w: del w[w.index(flag):w.index(flag)+2]
        template = ' '.join(w)
    for k in range(0, len(SET.get(name, [])), 2):
        w, (flag, value) = template.split(), SET[name][k:k+2]
        if flag in w: w[w.index(flag)+1] = value
        else: w += [flag, value]
        template = ' '.join(w)
    source, back = start[0], (0 if start[1] is None else back_of(*start))
    (link/'summary').mkdir(exist_ok=True)
    for k in range(MAX_ROUNDS):
        if k == 0 and (link/rnd).is_dir() and L.ranks(link/rnd):  # a restart of this queue: follow the round already running
            root = link/rnd
            log(f'{name}: following {rnd}, already running')
        else:
            root, err = new_round(link, source, back, rnd, template, cores)
            if err:
                log(f'{name}: {rnd}: {err}' + (f'; see {root.name}/run.log' if root else '')); return
            log(f'{name}: started {root.name} from {source.parent.name if source.parent != link else ""} {source.name} '
                f'(iterate {back} before the newest), cores {cores[0]}-{cores[-1]}')
            FQ.subprocess.Popen([FQ.PY, '-u', str(FQ.S/'live_plots.py'), '--out', str(link/'summary'), *map(str, chain_rounds(link))],
                                stdout=(link/'summary/live_plots.log').open('a'), stderr=FQ.subprocess.STDOUT, start_new_session=True)
        why, go_on = RC.watch(root, WATCH)
        st, r = kkt_state(root), last_record(root)
        its = sorted(d for d in (root/'iterates').glob('eval_*') if (d/'design.json').exists())
        summary = (f"{len(its)} iterates; last: loss {r['cpt_loss_outlet']:.6f}, C_L c {.5*r['cl']:.4f}, Cp_min {r.get('cp_min_ks', float('nan')):.3f}, "
                   f"feasible {r['feasible']}") if r else 'no converged evaluation'
        if st and st[0] == 4:
            log(f'{name}: KKT MET in {root.name}: PSQP exit 4 ({st[1]}); {summary}'); break
        log(f'{name}: {root.name} ended: {why}' + (f' (PSQP exit {st[0]}: {st[1]})' if st else '') + f'; {summary}')
        if not its:
            if not go_on or source is None: log(f'{name}: {root.name} saved no iterate; stopping, needs a look'); return
            rnd = RC.retry_name(root.name); continue  # same start again on the next name
        source, back, rnd = root, 0, RC.next_name(root.name)
    else:
        log(f'{name}: {MAX_ROUNDS} rounds without KKT')
    with (link/'summary/final_plots.log').open('w') as f:
        FQ.subprocess.run([FQ.PY, str(FQ.S/'plot_run_chain.py'), '--out', str(link/'summary'), *map(str, chain_rounds(link))],
                          stdout=f, stderr=FQ.subprocess.STDOUT)
    log(f'{name}: done')

def worker(name):
    try: drive(name)
    except Exception as e:
        import traceback; log(f'{name}: driver error {type(e).__name__}: {e}\n{traceback.format_exc()}')

if __name__ == '__main__':
    names = [n for n in RUNS if not sys.argv[1:] or any(a in n for a in sys.argv[1:])]
    log('KKT queue:', ', '.join(names))
    threads = [threading.Thread(target=worker, args=(n,)) for n in names]
    for t in threads: t.start(); time.sleep(20)  # stagger the remeshes
    for t in threads: t.join()

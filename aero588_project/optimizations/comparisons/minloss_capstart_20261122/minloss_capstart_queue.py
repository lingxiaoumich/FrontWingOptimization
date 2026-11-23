"""The min-loss runs started from a loss-cap design instead of the level start (user, 2026-11-22: "start from loss-cap
design"; and for 036/037 "change target to downforce 1.0 and 0.8"):
  run_031_area050MinLossClc18CapStart          C_L c <= -1.8
  run_032_area050MinLossClc16CapStart          C_L c <= -1.6
  run_033_area050MinLossClc14CapStart          C_L c <= -1.4
  run_036_area050MinLossClc10CpMin4CapStart    C_L c <= -1.0, Cp >= -4
  run_037_area050MinLossClc08CpMin4CapStart    C_L c <= -0.8, Cp >= -4
  run_038_area050MinLossClc14CpMin4CapStart    C_L c <= -1.4, Cp >= -4
Why: from the level start, run_034 (C_L c <= -1.2) ended at loss 0.00404 and run_039 (+ Cp >= -4) at 0.00308, both at
C_L c -1.199, while run_030 (maximum downforce, loss <= 0.003067) had passed a strictly feasible design at C_L c -1.204
with loss 0.00295: the min-loss runs settle in poorer designs. The level-start runs 031/032 were stopped for these
(their folders keep the old names, without CapStart).
Start: among the saved PSQP iterates (each has its flow) of the loss-cap runs on the same limits (run_019/020/021/030,
all run_016's setup: area >= 50 %, the reduced box, the 50 mm floor), the strictly feasible one with C_L c <= the
limit and the lowest outlet loss. It is remeshed as a chain round is (run_chain.prepare: the grid's zero point is
that design, the initial flow is its flow mapped onto the new grid), and the command is the one of runs 031-039
(lanes.COMMON with --drag-cap 10 --max-clc T [--cp-min -4]). All limits are absolute, so they are the same as there.
The Cp runs start from designs whose Cp is not known (the loss-cap runs had no Cp output); PSQP brings Cp to the limit.
Run by comparisons/priority_20261121/priority_queue.py (job kind 'minlosscap', gated on READY here, written once
build() made all six grids in the scratchpad).
"""
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
for d in ('floor_20261121', 'minloss_20261121', 'minloss_cpmin4_20261121'):
    sys.path.insert(0, str(HERE.parent/d))
import floor_queue as FQ   # chain(), log(), RANKFILE; FQ.L is lanes.py, FQ.RC run_chain
import minloss_queue as M
import minloss_cpmin4_queue as MC
L, RC = FQ.L, FQ.RC
POOL = ['run_019_area050CptLoss075', 'run_020_area050CptLoss050', 'run_021_area050CptLoss025', 'run_030_area050CptLoss015']
CP = MC.CP

def iterates():
    """(round, back, iterate folder, result) for every saved iterate of the pool; `back` as in run_chain.last_iterate."""
    for R in POOL:
        for r in sorted((L.O/R).glob('r[0-9]*')):
            saved = sorted(d for d in (r/'iterates').glob('eval_*') if (d/'design.json').exists())
            for k, d in enumerate(saved):
                yield r, len(saved)-1-k, d, json.loads((d/'result.json').read_text())

def start(t):
    ok = [(q['cpt_loss_outlet'], str(d), r, back, d, q) for r, back, d, q in iterates()
          if q.get('feasible') and not q.get('cfd_failed') and 0.5*q['cl'] <= t and (d/'flow_000_vol.cgns').exists()]
    return min(ok)[2:] if ok else None

def build(name, t, cp, lane, test=None, rnd='r00', l2=None):
    """Make the run folder (or `test`, a scratch folder) with round `rnd` (r00, or a retry of it in the existing folder, e.g.
    r00b) on the remeshed start, with --l2-convergence `l2` if given; returns (link, root, what)."""
    s = start(t)
    if s is None: raise RuntimeError(f'no strictly feasible loss-cap iterate with C_L c <= {t:g}')
    run, back, it, q = s
    what = (f"{run.parent.name} {run.name} {it.name} (C_L c {0.5*q['cl']:.4f}, loss {q['cpt_loss_outlet']:.6f}, "
            f"C_D c {0.5*q['cd']:.4f}, area {q['section_area_m2']/L.A0:.1%})")
    if test: link = Path(test); link.mkdir(parents=True)
    elif rnd != 'r00': link = L.O/name
    else:
        link, real = L.O/name, L.BIG/name
        if link.exists(): raise RuntimeError(f'{link} exists already')
        real.mkdir(parents=True); link.symlink_to(real)
    (link/'summary').mkdir(exist_ok=rnd != 'r00')
    new = RC.prepare(run, link, rnd, [], back)
    if new is None: raise RuntimeError(f'the remesh of {what} failed')
    root, _, hyp = new
    mesh = link/'grids'/f"{rnd}_from_{run.name}_{it.name.replace('_', '')}"
    if (mesh/'source.txt').read_text().split()[-1] != str(q['evaluation']): raise RuntimeError(f'prepare() took another iterate than {it}')
    (mesh/'source.txt').write_text(f"{run.parent.name} {run.name} evaluation {q['evaluation']} ({it})\n")
    for f in (root/'source').glob('*.py'):
        if (L.S/f.name).exists(): shutil.copy(L.S/f.name, f)  # the current scripts
    shutil.copy(FQ.RANKFILE[lane], root/'rankfile')
    cap = f'--drag-cap 10 --max-clc {t:g}' + (f' --cp-min {CP}' if cp else '')
    cmd = (f'mpirun -np 24 --map-by rankfile:file=rankfile --bind-to core --mca btl self,sm -x GFORTRAN_UNBUFFERED_PRECONNECTED=y '
           f'python -u {root}/source/run_slip_optimization.py --root {root} ' + L.COMMON.format(cap=cap, amin=L.AMIN) +
           f' --baseline {mesh}/baseline.json --check-design {mesh}/offset_design.json' + (f' --l2-convergence {l2}' if l2 else ''))
    cores = '12-35' if lane == 'A' else '36-59'
    (root/'command.txt').write_text(
        cmd + f"\n# {name}: as runs 031-039 ({L.SETUP.split(', section')[0].replace('level NACA 0012 stretched to the box x-ends, ', '')}, "
        f"section area >= 50 % of the level start, no drag cap; the outlet Cpt loss minimized with C_L c <= {t:g}"
        + (f', Cp >= {CP}' if cp else '') + f"), but started from the loss-cap design {what} instead of the level start: "
        f"this grid is its remesh (pyHyp {hyp}), source/restart_vol.cgns its flow. 24 ranks on cores {cores}, lane {lane} of "
        f"comparisons/priority_20261121/priority_queue.py via comparisons/minloss_capstart_20261122/minloss_capstart_queue.py, "
        f"{time.strftime('%c')}." + ('' if test else f' Folder on the 4T HD ({L.BIG/name}), linked from optimizations/') + '\n')
    return link, root, what

def run(name, t, cp, lane):
    link, root, what = build(name, t, cp, lane)
    FQ.log(f'{name}: started r00 (--max-clc {t:g}' + (f' --cp-min {CP}' if cp else '') + f') from {what}')
    FQ.log(f'{name}: {FQ.chain(link, [root], lane)}')
    (MC if cp else M).finish(link, name, t)

def finish(link, name, t, cp, rounds):
    """As minloss_queue/minloss_cpmin4_queue.finish(), on the given rounds plus those run_chain.py started."""
    text = (link/'summary/run_chain.log').read_text()
    rounds = rounds + [link/l.split()[-1] for l in text.splitlines() if ' started r' in l and link/l.split()[-1] not in rounds]
    rounds = [r for r in dict.fromkeys(rounds) if (r/'history.jsonl').exists()]
    time.sleep(120)  # live_plots.py makes its last refresh when the run stops
    with (link/'summary/final_plots.log').open('w') as f:
        FQ.subprocess.run([FQ.PY, str(FQ.S/'plot_run_chain.py'), '--out', str(link/'summary'), *map(str, rounds)], stdout=f, stderr=FQ.subprocess.STDOUT)
    rows = [json.loads((r/'best_feasible/result.json').read_text()) | {'round': r.name} for r in rounds if (r/'best_feasible/result.json').exists()]
    b = min(rows, key=lambda q: q['cpt_loss_outlet']) if rows else None
    if b: FQ.log(f"{name}: best loss {b['cpt_loss_outlet']:.6f} at C_L c {.5*b['cl']:.3f} (limit {t:g})"
                 + (f", Cp min (KS) {b.get('cp_min_ks', float('nan')):.2f} (limit {CP})" if cp else '')
                 + f", C_D c {.5*b['cd']:.4f}, area {b['section_area_m2']/L.A0:.1%} at {b['round']} evaluation {b['evaluation']}")
    else: FQ.log(f'{name}: no strictly feasible design')

"""Suction-peak run (user, 2026-11-25: "optimize for reduction of suction peak Cp, while keeping Clc <= 1.2"):
  run_042_area050MinSuctionClc12CapStart: minimize the suction peak -Cp_min (ADflow's KS cpmin, rho 100, over the wing)
  with C_L c <= -1.2 (downforce >= 1.2), otherwise run_016's limits (area >= 50 % of the level start, the reduced box,
  the 50 mm floor, no drag cap): run_slip_optimization.py --drag-cap 10 --max-clc -1.2 --min-suction.
Start (user's choice): run_027_area050CpMin4 r00 evaluation 41, its best strictly feasible design (maximum downforce with
Cp >= -4: C_L c -1.393, Cp_min -3.98), remeshed as a chain round is (run_chain.prepare: the grid's zero point is that
design, the initial flow is its flow). It is past the lift limit, so PSQP can approach C_L c -1.2 from the feasible
side: from the level start, runs 034/039 (C_L c <= -1.2) ended at -1.199 and never became strictly feasible.
If the solver reports orphan cells, the start is remeshed with the next pyHyp setting as r00b, ...
It waits for the first lane that comparisons/rotstart_20261125/rotstart_queue.py frees ('lane X done' in its log).
Run with the MACH-Aero env (conda activate mdolab), detached, helpers on cores 12-63:
  setsid nohup taskset -c 12-63 python -u minsuction_queue.py > queue.log 2>&1 &
  python minsuction_queue.py --dry-run DIR   builds the start round under DIR only (no lane, no run)
"""
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent/'floor_20261121'))
import floor_queue as FQ   # log(), RANKFILE, PY, S; FQ.L is lanes.py, FQ.RC run_chain
L, RC = FQ.L, FQ.RC
NAME = 'run_042_area050MinSuctionClc12CapStart'
T = -1.2
SOURCE = L.O/'run_027_area050CpMin4/r00'
SOURCE_EVAL = 41
ROT_LOG = HERE.parent/'rotstart_20261125/queue.log'

def source_iterate():
    """(back, iterate folder, result) of SOURCE_EVAL; `back` as in run_chain.last_iterate."""
    saved = sorted(d for d in (SOURCE/'iterates').glob('eval_*') if (d/'design.json').exists())
    it = SOURCE/'iterates'/f'eval_{SOURCE_EVAL:03d}'
    q = json.loads((it/'result.json').read_text())
    assert it in saved and q['feasible'] and not q['cfd_failed'] and .5*q['cl'] <= T and (it/'flow_000_vol.cgns').exists()
    return len(saved)-1-saved.index(it), it, q

def build(link, lane, rnd='r00', skip=()):
    """Round `rnd` of the run in `link` on the remeshed start; returns (root, pyHyp setting, what) or None."""
    back, it, q = source_iterate()
    what = (f"run_027_area050CpMin4 r00 evaluation {SOURCE_EVAL} (C_L c {.5*q['cl']:.4f}, Cp_min (KS) {q['cp_min_ks']:.3f}, "
            f"C_D c {.5*q['cd']:.4f}, area {q['section_area_m2']/L.A0:.1%})")
    new = RC.prepare(SOURCE, link, rnd, list(skip), back)
    if new is None: return None
    root, _, hyp = new
    mesh = link/'grids'/f"{rnd}_from_{SOURCE.name}_{it.name.replace('_', '')}"
    if (mesh/'source.txt').read_text().split()[-1] != str(SOURCE_EVAL): raise RuntimeError(f'prepare() took another iterate than {it}')
    (mesh/'source.txt').write_text(f"run_027_area050CpMin4 r00 evaluation {SOURCE_EVAL} ({it})\n")
    for f in (root/'source').glob('*.py'):
        if (L.S/f.name).exists(): shutil.copy(L.S/f.name, f)  # the current scripts
    shutil.copy(FQ.RANKFILE[lane], root/'rankfile')
    cap = f'--drag-cap 10 --max-clc {T:g} --min-suction'
    cmd = (f'mpirun -np 24 --map-by rankfile:file=rankfile --bind-to core --mca btl self,sm -x GFORTRAN_UNBUFFERED_PRECONNECTED=y '
           f'python -u {root}/source/run_slip_optimization.py --root {root} ' + L.COMMON.format(cap=cap, amin=L.AMIN) +
           f' --baseline {mesh}/baseline.json --check-design {mesh}/offset_design.json')
    cores = '12-35' if lane == 'A' else '36-59'
    (root/'command.txt').write_text(
        cmd + f"\n# {NAME}: the suction peak -Cp_min (ADflow KS cpmin, rho 100) minimized with C_L c <= {T:g}; otherwise "
        f"run_016_area050's limits (20 FFD pairs moving only up and down, tilt locked, the reduced box, section area >= "
        f"{L.AMIN:.9f} m^2 = 0.50 x the level start's, no drag cap). Started from {what}, the best design of maximum downforce "
        f"with Cp >= -4: this grid is its remesh (pyHyp {hyp}), source/restart_vol.cgns its flow. 24 ranks on cores {cores}, "
        f"lane {lane} of comparisons/minsuction_20261125/minsuction_queue.py, {time.strftime('%c')}. Folder {link}\n")
    return root, hyp, what

def free_lane():
    """The first lane rotstart_queue.py has finished with."""
    while True:
        text = ROT_LOG.read_text() if ROT_LOG.exists() else ''
        for lane in 'AB':
            if f'lane {lane} done' in text: return lane
        time.sleep(60)

def finish(link, rounds):
    text = (link/'summary/run_chain.log').read_text()
    rounds = rounds + [link/l.split()[-1] for l in text.splitlines() if ' started r' in l and link/l.split()[-1] not in rounds]
    rounds = [r for r in dict.fromkeys(rounds) if (r/'history.jsonl').exists()]
    time.sleep(120)  # live_plots.py makes its last refresh when the run stops
    with (link/'summary/final_plots.log').open('w') as f:
        FQ.subprocess.run([FQ.PY, str(FQ.S/'plot_run_chain.py'), '--out', str(link/'summary'), *map(str, rounds)],
                          stdout=f, stderr=FQ.subprocess.STDOUT)
    rows = [json.loads((r/'best_feasible/result.json').read_text()) | {'round': r.name} for r in rounds if (r/'best_feasible/result.json').exists()]
    b = max(rows, key=lambda q: q['cp_min_ks']) if rows else None
    if b: FQ.log(f"{NAME}: best Cp_min (KS) {b['cp_min_ks']:.3f} at C_L c {.5*b['cl']:.3f} (limit {T:g}), C_D c {.5*b['cd']:.4f}, "
                 f"loss {b['cpt_loss_outlet']:.6f}, area {b['section_area_m2']/L.A0:.1%} at {b['round']} evaluation {b['evaluation']}")
    else: FQ.log(f'{NAME}: no strictly feasible design')

def main():
    lane = free_lane()
    FQ.log(f'{NAME}: lane {lane} is free')
    link = L.O/NAME
    if link.exists(): raise RuntimeError(f'{link} exists already')
    (link/'summary').mkdir(parents=True)
    skip, rnd = [], 'r00'
    while True:
        made = build(link, lane, rnd, skip)
        if made is None: FQ.log(f'{NAME}: no pyHyp setting gave a usable start grid'); return
        root, hyp, what = made
        FQ.subprocess.Popen((root/'command.txt').read_text().split('\n')[0], shell=True, cwd=root,
                            stdout=(root/'run.log').open('w'), stderr=FQ.subprocess.STDOUT, start_new_session=True)
        FQ.log(f'{NAME}: started {rnd} (pyHyp {hyp}) from {what}')
        t = time.time()
        while not L.ranks(root) and time.time()-t < 600: time.sleep(10)
        if not L.ranks(root): FQ.log(f'{NAME}: the ranks never came up; see {rnd}/run.log'); return
        n = RC.orphans(root)
        if not n: break
        RC.stop(root, f'{n} orphan cells in the overset grid; the adjoint cannot converge with them')
        FQ.log(f'{NAME}: {n} orphan cells on the start grid; trying the next pyHyp setting')
        skip.append(hyp); rnd = RC.retry_name(rnd)
    time.sleep(30)
    FQ.subprocess.Popen([FQ.PY, '-u', str(FQ.S/'live_plots.py'), '--out', str(link/'summary'), str(root)],
                        stdout=(link/'summary/live_plots.log').open('a'), stderr=FQ.subprocess.STDOUT, start_new_session=True)
    with (link/'summary/run_chain.log').open('w') as f:
        FQ.subprocess.run([FQ.PY, '-u', str(FQ.S/'run_chain.py'), '--live-plots', str(link/'summary'), '--stall-minutes', '15', str(root)],
                          stdout=f, stderr=FQ.subprocess.STDOUT)
    FQ.log(f"{NAME}: {(link/'summary/run_chain.log').read_text().strip().splitlines()[-1]}")
    finish(link, [root])

if __name__ == '__main__':
    if sys.argv[1:2] == ['--dry-run']:
        d = Path(sys.argv[2]); (d/'summary').mkdir(parents=True)
        print(build(d, 'A'))
    else:
        try: main()
        except Exception as e: FQ.log(f'{NAME}: queue error {type(e).__name__}: {e}')

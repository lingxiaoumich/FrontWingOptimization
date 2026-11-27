"""Collect per-run facts (setup, best feasible design, meshes, times) for the report: python report_gather.py runs.json
(MACH-Aero env: the KKT check reads the optimizer databases with pyoptsparse).

Since 2026-11-27 (report_20261127): runs 036-039 were run to PSQP's KKT test (comparisons/kkt_20261125). Their best
design is the KKT point named in KKT, not the lowest loss over all rounds: the rounds before 11-26 measured the loss
with the old planes (9-12 % low) and a moving far field, so their values are not comparable. Each KKT entry also gets
PSQP's exit and final G, the independent residual of kkt_residual.py, and the cost of all its rounds (branches too).
A run with the --min-suction objective (042) keeps its least-suction feasible design."""
import json, re, glob, os
from pathlib import Path
import numpy as np
import kkt_residual as KR
O = Path('/home/lingxiao77/AER588/aero588_project/optimizations')
A0 = 0.021549464275835817
KKT = {'run_036_area050MinLossClc10CpMin4CapStart': ['r07b'], 'run_037_area050MinLossClc08CpMin4CapStart': ['r08'],
       'run_038_area050MinLossClc14CpMin4CapStart': ['branch_xmax020/r10t', 'branch_xmax050/r11t'],
       'run_039_area050MinLossClc12CpMin4': ['r08']}  # the first round of each list is the run's reported design
out = {}

def best_dict(best, rnd):
    return dict(round=rnd, evaluation=best['evaluation'], clc=0.5*best['cl'], cdc=0.5*best['cd'],
              loss=best.get('cpt_loss_outlet'), loss_flux=best.get('cpt_loss_outlet_flux_c_m'), area=best.get('section_area_m2'), area_pct=100*best['section_area_m2']/A0 if best.get('section_area_m2') else None,
              clearance_mm=1e3*best['minimum_clearance_m'], cp_min_ks=best.get('cp_min_ks'), incidence=best.get('chord_incidence_deg'),
              te_base_mm=1e3*best['te_base_thickness_m'] if best.get('te_base_thickness_m') is not None else None,
              thinnest_mm=1e3*best['minimum_surface_gap_m'] if best.get('minimum_surface_gap_m') is not None else None,
              max_height_mm=1e3*best['maximum_height_m'] if best.get('maximum_height_m') else None)

def kkt_point(link, rnd):
    """A KKT round's reported design (its best feasible evaluation), PSQP's exit and G, and the independent residual."""
    root = link/rnd
    q = json.loads((root/'best_feasible/result.json').read_text())
    inform = json.loads((root/'run_complete.json').read_text())['inform']
    psqp = (root/'PSQP.out').read_text()  # the final line (with ITERM) follows 'EXIT FROM PSQP'
    call, obj, res, mult, n_act = KR.last_points(root, 1)[0]
    st = json.loads((root/'settings.json').read_text())
    return best_dict(q, rnd) | dict(psqp_exit=inform.get('value'), psqp_g=float(re.findall(r'\sG=\s*(\S+)', psqp)[-1]),
                                    psqp_iterations=int(re.findall(r'NIT=\s*(\d+)', psqp)[-1]), kkt_residual=res, n_active=n_act,
                                    multipliers=[(a, round(b, 4)) for a, b in mult], l2=st.get('l2_convergence'))

def cost(link):
    """Evaluations and wall hours of every round under the run (branches too), and of the KKT push alone (rounds that
    kkt_queue.py started, or hand-launched same-grid continuations)."""
    c = dict(evals=0, wall_h=0., core_h=0., kkt_evals=0, kkt_wall_h=0., kkt_core_h=0., kkt_rounds=0)
    for hj in sorted(link.rglob('history.jsonl')):
        H = [json.loads(l) for l in open(hj) if l.strip()]
        if not H: continue
        cmd = (hj.parent/'command.txt').read_text() if (hj.parent/'command.txt').exists() else ''
        w = cmd.split('\n')[0].split(); np_ = int(w[w.index('-np')+1]) if '-np' in w else 24; h = H[-1]['elapsed_wall_time_s']/3600
        c['evals'] += len(H); c['wall_h'] += h; c['core_h'] += h*np_
        if 'KKT continuation' in cmd or 'same-grid continuation' in cmd:
            c['kkt_evals'] += len(H); c['kkt_wall_h'] += h; c['kkt_core_h'] += h*np_; c['kkt_rounds'] += 1
    return c
for link in sorted(O.glob('run_0*')):
    R = link.name
    cs = link/'summary/chain_summary.json'
    rounds = [link/r['run'] for r in json.loads(cs.read_text())['runs']] if cs.exists() else sorted(p for p in link.glob('r0*') if (p/'history.jsonl').exists())
    rounds = [r for r in rounds if (r/'history.jsonl').exists()]
    if not rounds: out[R] = dict(rounds=[]); continue
    st = json.loads((rounds[0]/'settings.json').read_text())
    cmd = (rounds[0]/'command.txt').read_text().split('\n')[0].split()
    opt = lambda k: cmd[cmd.index(k)+1] if k in cmd else None
    minloss = st.get('max_clc_m') is not None and not st.get('min_suction')
    minsuction = bool(st.get('min_suction'))
    rows, evals, ft, adj, wall, solves, failed = [], 0, [], [], 0., 0, 0
    for r in rounds:
        H = [json.loads(l) for l in open(r/'history.jsonl')]
        evals += len(H); wall += H[-1]['elapsed_wall_time_s']
        ft += [h['evaluation_wall_time_s'] for h in H if not h['cfd_failed']]
        solves += sum(h.get('solves', 1) for h in H); failed += sum(h['cfd_failed'] for h in H)
        if (r/'run.log').exists():
            adj += [float(m) for m in re.findall(r'Solving ADjoint Transpose with PETSc time \(s\) =\s*([\d.]+)', (r/'run.log').read_text(errors='ignore'))]
        b = r/'best_feasible/result.json'
        if b.exists():
            q = json.loads(b.read_text()); q['round'] = r.name; rows.append(q)
    best = None
    if rows:
        if minsuction: best = max(rows, key=lambda q: q['cp_min_ks'])
        elif minloss:
            rows2 = [q for q in rows if not q.get('zero_deformation')] or rows
            best = min(rows2, key=lambda q: q['cpt_loss_outlet'])
        else: best = min(rows, key=lambda q: q['cl'])
    start = json.loads(open(rounds[0]/'history.jsonl').readline())
    its = json.loads(cs.read_text())['iterates'] if cs.exists() else []
    out[R] = dict(
        rounds=[r.name for r in rounds], evals=evals, wall_h=wall/3600, solves=solves, failed=failed,
        flow_s_median=float(np.median(ft)) if ft else None, adjoint_s_median=float(np.median(adj)) if adj else None, n_adjoint=len(adj),
        iterates=len(its),
        box=st['box_m'], drag_cap_cdc=(0.5*st['drag_cap'] if st.get('drag_cap', 10) < 5 else None), loss_cap=st.get('loss_cap'),
        cp_min=st.get('cp_min'), max_clc=st.get('max_clc_m'), min_area=st.get('min_section_area_m2'),
        shape_directions=st.get('shape_directions'), ffd=st.get('ffd_chordwise_stations'), incidence_fixed=st.get('incidence_fixed'),
        ground=st.get('ground_boundary'), cells=st.get('discretization_cells_before_overset_blanking'), ranks=st.get('mpi_ranks'),
        l2=opt('--l2-convergence'),
        start=dict(clc=0.5*start['cl'], cdc=0.5*start['cd'], loss=start.get('cpt_loss_outlet'), area=start.get('section_area_m2'),
                   incidence=start.get('chord_incidence_deg')),
        best=None if best is None else best_dict(best, best['round']), min_suction=minsuction,
        last=dict(clc=its[-1]['cl_c'], cdc=its[-1]['cd_c']) if its else None,
        running=bool(os.popen(f"ps -eo args | grep -c '[p]ython -u {link}/'").read().strip() not in ('0', '')),
    )
    if R in KKT:
        pts = [kkt_point(link, r) for r in KKT[R]]
        out[R].update(best=pts[0], kkt_alternatives=pts[1:], **cost(link))
json.dump(out, open(__import__('sys').argv[1], 'w'), indent=1)
for R, d in out.items():
    if not d.get('rounds'): print(R, 'no rounds'); continue
    b = d['best']
    print(f"{R:44s} {len(d['rounds'])}r {d['evals']:4d}ev {d['wall_h']:4.1f}h flow {d['flow_s_median'] or 0:5.1f}s adj {d['adjoint_s_median'] or 0:5.1f}s "
          f"dirs {d['shape_directions']} ffd {d['ffd']} fixinc {d['incidence_fixed']} cells {d['cells']} ranks {d['ranks']} "
          + (f"best {b['clc']:.3f}/{b['cdc']:.4f} loss {b['loss'] if b['loss'] is None else round(b['loss'],5)} area {b['area_pct'] or 0:.0f}% clr {b['clearance_mm']:.1f} cp {b['cp_min_ks'] if b['cp_min_ks'] is None else round(b['cp_min_ks'],2)} {b['round']}" if b else 'no best')
          + (' RUNNING' if d['running'] else '')
          + (f" KKT exit {b['psqp_exit']} G {b['psqp_g']:.1e} residual {b['kkt_residual']:.2e}; push {d['kkt_evals']} ev {d['kkt_wall_h']:.1f} h" if b and 'psqp_exit' in b else ''))

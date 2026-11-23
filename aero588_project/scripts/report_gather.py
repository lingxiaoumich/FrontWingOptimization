"""Collect per-run facts (setup, best feasible design, meshes, times) for the report: python report_gather.py runs.json"""
import json, re, glob, os
from pathlib import Path
import numpy as np
O = Path('/home/lingxiao77/AER588/aero588_project/optimizations')
A0 = 0.021549464275835817
out = {}
for link in sorted(O.glob('run_0*')):
    R = link.name
    cs = link/'summary/chain_summary.json'
    rounds = [link/r['run'] for r in json.loads(cs.read_text())['runs']] if cs.exists() else sorted(p for p in link.glob('r0*') if (p/'history.jsonl').exists())
    rounds = [r for r in rounds if (r/'history.jsonl').exists()]
    if not rounds: out[R] = dict(rounds=[]); continue
    st = json.loads((rounds[0]/'settings.json').read_text())
    cmd = (rounds[0]/'command.txt').read_text().split('\n')[0].split()
    opt = lambda k: cmd[cmd.index(k)+1] if k in cmd else None
    minloss = st.get('max_clc_m') is not None
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
        if minloss:
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
        best=None if best is None else dict(round=best['round'], evaluation=best['evaluation'], clc=0.5*best['cl'], cdc=0.5*best['cd'],
              loss=best.get('cpt_loss_outlet'), loss_flux=best.get('cpt_loss_outlet_flux_c_m'), area=best.get('section_area_m2'), area_pct=100*best['section_area_m2']/A0 if best.get('section_area_m2') else None,
              clearance_mm=1e3*best['minimum_clearance_m'], cp_min_ks=best.get('cp_min_ks'), incidence=best.get('chord_incidence_deg'),
              te_base_mm=1e3*best['te_base_thickness_m'] if best.get('te_base_thickness_m') is not None else None,
              thinnest_mm=1e3*best['minimum_surface_gap_m'] if best.get('minimum_surface_gap_m') is not None else None,
              max_height_mm=1e3*best['maximum_height_m'] if best.get('maximum_height_m') else None),
        last=dict(clc=its[-1]['cl_c'], cdc=its[-1]['cd_c']) if its else None,
        running=bool(os.popen(f"ps -eo args | grep -c '[p]ython -u {link}/'").read().strip() not in ('0', '')),
    )
json.dump(out, open(__import__('sys').argv[1], 'w'), indent=1)
for R, d in out.items():
    if not d.get('rounds'): print(R, 'no rounds'); continue
    b = d['best']
    print(f"{R:44s} {len(d['rounds'])}r {d['evals']:4d}ev {d['wall_h']:4.1f}h flow {d['flow_s_median'] or 0:5.1f}s adj {d['adjoint_s_median'] or 0:5.1f}s "
          f"dirs {d['shape_directions']} ffd {d['ffd']} fixinc {d['incidence_fixed']} cells {d['cells']} ranks {d['ranks']} "
          + (f"best {b['clc']:.3f}/{b['cdc']:.4f} loss {b['loss'] if b['loss'] is None else round(b['loss'],5)} area {b['area_pct'] or 0:.0f}% clr {b['clearance_mm']:.1f} cp {b['cp_min_ks'] if b['cp_min_ks'] is None else round(b['cp_min_ks'],2)} {b['round']}" if b else 'no best')
          + (' RUNNING' if d['running'] else ''))

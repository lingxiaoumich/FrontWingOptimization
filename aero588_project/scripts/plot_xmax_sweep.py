"""Compare short PSQP restarts that differ only in the step limit XMAX.

Reads every <sweep>/xmax_*/ run made by run_xmax_sweep.sh and writes
<sweep>/xmax_sweep.{png,svg} and <sweep>/xmax_sweep_summary.{json,txt}.
Forces are reported as C_L*c = L'/q and C_D*c = D'/q (m), independent of the
reference chord. Accepted iterates come from the NIT/NFV lines of PSQP.out;
evaluations between them are line-search trials.
"""
import argparse
import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('--sweep',type=Path,required=True)
a=p.parse_args();sweep=a.sweep
ramp=['#86b6ef','#3987e5','#1c5cab','#0d366b','#082447']
ink,muted,grid='#2b2b2b','#6b6b6b','#e4e4e4'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'figure.facecolor':'white','savefig.facecolor':'white','font.family':'DejaVu Sans'})

def scaled(design):
    return np.r_[design['incidence'],100*np.asarray(design['shape_x']),100*np.asarray(design['shape_y'])]

def load(run):
    settings=json.loads((run/'settings.json').read_text())
    c0=settings['reference_chord_m'];xmax=settings['optimizer_options']['XMAX']
    history=[json.loads(x) for x in (run/'history.jsonl').read_text().splitlines() if x.strip()]
    lines=re.findall(r'NIT=\s*(\d+)\s+NFV=\s*(\d+)\s+NFG=\s*(\d+)\s+F=\s*(\S+)',(run/'PSQP.out').read_text(errors='replace'))
    nit={int(i):int(v) for i,v,_,_ in lines}  # last NFV reported for each major iteration
    complete=(run/'run_complete.json').exists()
    nfv=max(nit.values())
    # optimizer evaluations are the last nfv records (plus the final re-solve when complete)
    end=len(history)-1 if complete else len(history)
    opt=history[end-nfv:end]
    acc=[opt[nit[i]-1] for i in sorted(nit)]
    for i,r in zip(sorted(nit),acc):
        f=float(dict((int(n),float(F)) for n,_,_,F in lines)[i])
        assert abs(r['cl']-f)<5e-5,(run.name,i,r['cl'],f)
    X=np.array([scaled(r['design']) for r in acc])
    step=np.linalg.norm(np.diff(X,axis=0),axis=1)
    exit_=json.loads((run/'optimizer_exit.json').read_text())['inform'] if (run/'optimizer_exit.json').exists() else None
    t0=opt[0]['elapsed_wall_time_s']-opt[0]['evaluation_wall_time_s']
    return dict(name=run.name,xmax=xmax,c0=c0,complete=complete,exit=exit_,
        it=np.array(sorted(nit)),clc=c0*np.array([r['cl'] for r in acc]),cdc=c0*np.array([r['cd'] for r in acc]),
        t=np.array([(r['elapsed_wall_time_s']-t0)/60 for r in acc]),step=step,
        evals=len(opt),failed=sum(r['cfd_failed'] for r in opt),
        trials=len(opt)-len(acc),last=acc[-1],first=acc[0])

runs=sorted((load(r) for r in sweep.glob('xmax_*') if (r/'PSQP.out').exists()),key=lambda r:r['xmax'])
cap=.05*runs[0]['c0']
fig,axs=plt.subplots(2,2,figsize=(12,8.4),constrained_layout=True)
shade=[ramp[int(round(i))] for i in np.linspace(0,len(ramp)-1,len(runs))] if len(runs)>1 else ramp[-1:]
for r,c in zip(runs,shade):
    lab=f"XMAX = {r['xmax']:g}"+('' if r['complete'] else f" (incomplete: {int(r['it'][-1])} iterations)")
    axs[0,0].plot(r['it'],-r['clc'],'-o',c=c,ms=4,lw=1.6,label=lab)
    axs[0,1].plot(r['t'],-r['clc'],'-o',c=c,ms=4,lw=1.6,label=lab)
    axs[1,0].plot(-r['clc'],r['cdc'],'-o',c=c,ms=4,lw=1.6,label=lab)
    if len(r['step']):
        axs[1,1].plot(r['it'][1:],r['step'],'-o',c=c,ms=4,lw=1.6,label=lab)
        axs[1,1].axhline(r['xmax'],c=c,ls=':',lw=1)
axs[0,0].set(xlabel='PSQP major iteration after restart',ylabel="Downforce $-C_L c = -L'/q_\\infty$ (m)",
             title='(a) Downforce per major iteration')
axs[0,1].set(xlabel='Wall time since restart (min, 4 MPI ranks)',ylabel="$-C_L c$ (m)",
             title='(b) Downforce per unit wall time')
axs[1,0].axhline(cap,c=muted,ls='--',lw=1)
axs[1,0].text(axs[1,0].get_xlim()[0],cap,f'  $C_D c$ cap {cap:.3f} m',va='bottom',color=muted,fontsize=9)
axs[1,0].set(xlabel="$-C_L c$ (m)",ylabel="$C_D c = D'/q_\\infty$ (m)",title='(c) Drag polar of accepted iterates')
axs[1,1].set(ylim=(0,1.15*max(r['xmax'] for r in runs)),xlabel='PSQP major iteration after restart',ylabel='Step length, scaled variables',
             title='(d) Step taken (dotted: XMAX limit)')
for ax in axs.flat:ax.grid(c=grid,lw=.6);ax.set_axisbelow(True)
axs[0,0].legend(frameon=False,fontsize=9)
fig.suptitle('PSQP step-limit sweep: restarts from the October 29 final design (slip ground, 88,548-cell mesh)',fontsize=11)
fig.savefig(sweep/'xmax_sweep.png',dpi=190,bbox_inches='tight');fig.savefig(sweep/'xmax_sweep.svg',bbox_inches='tight')

rows=[];summary={}
for r in runs:
    L=r['last'];F=r['first']
    s=dict(xmax=r['xmax'],complete=r['complete'],optimizer_exit=r['exit'],major_iterations=int(r['it'][-1]),
        cfd_evaluations=r['evals'],line_search_trials=r['trials'],failed_cfd=int(r['failed']),
        start_cl_c_m=float(r['clc'][0]),end_cl_c_m=float(r['clc'][-1]),
        downforce_gain_pct=float(100*(r['clc'][-1]/r['clc'][0]-1)),
        end_cd_c_m=float(r['cdc'][-1]),end_cd_c_over_cap=float(r['cdc'][-1]/cap),
        steps_at_limit=int(np.sum(r['step']>.999*r['xmax'])),steps=len(r['step']),
        wall_time_min=float(r['t'][-1]),
        end_actual_chord_m=L['actual_chord_m'],end_min_clearance_m=L['minimum_clearance_m'],
        end_max_thickness_m=L['maximum_section_thickness_m'],end_box_margin_m=L['box_minimum_margin_m'],
        end_min_x_m=L['minimum_x_m'],end_max_x_m=L['maximum_x_m'],end_incidence_deg=L['chord_incidence_deg'])
    summary[r['name']]=s
    rows.append(f"{r['xmax']:>6g} {s['major_iterations']:>4d} {s['cfd_evaluations']:>5d} {s['line_search_trials']:>6d} {s['failed_cfd']:>6d}"
                f" {s['end_cl_c_m']:>9.4f} {s['downforce_gain_pct']:>+7.2f}% {s['end_cd_c_m']:>8.5f} {s['steps_at_limit']:>3d}/{s['steps']:<3d}"
                f" {s['end_actual_chord_m']:>7.4f} {s['end_min_clearance_m']:>8.4f} {s['end_max_thickness_m']:>8.4f} {s['wall_time_min']:>6.1f}"
                f"  {'' if r['complete'] else 'incomplete'}")
head=(f"{'XMAX':>6} {'NIT':>4} {'evals':>5} {'trials':>6} {'failed':>6} {'C_L c(m)':>9} {'gain':>8} {'C_D c(m)':>8} {'lim':>7}"
      f" {'chord':>7} {'clear.':>8} {'thick.':>8} {'min':>6}")
(sweep/'xmax_sweep_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(sweep/'xmax_sweep_summary.txt').write_text(head+'\n'+'\n'.join(rows)+'\n')
print(head);print('\n'.join(rows))

"""Create before/after CFD figures and a human-readable optimization report."""
import argparse
import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
p.add_argument('--after',default='best_feasible');a=p.parse_args();root=a.root
out=root/'plots';out.mkdir(exist_ok=True)
labels=['baseline',a.after]
results=[json.loads((root/l/'result.json').read_text()) for l in labels]
fields=[np.load(root/l/'fields.npz') for l in labels]
surfaces=[np.loadtxt(root/l/'surface_geometry.csv',delimiter=',',skiprows=1) for l in labels]
cp=[np.loadtxt(root/l/'surface_cp.csv',delimiter=',',skiprows=1) for l in labels]
settings=json.loads((root/'settings.json').read_text())
history=[json.loads(x) for x in (root/'history.jsonl').read_text().splitlines() if x.strip()]
colors=['#175676','#cf642b'];titles=['Initial NACA 0012','Optimized airfoil']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'figure.facecolor':'white','savefig.facecolor':'white','font.family':'DejaVu Sans'})
def save(fig,name):
    fig.savefig(out/(name+'.png'),dpi=190,bbox_inches='tight')
    fig.savefig(out/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def body(ax,k,box=False):
    s=surfaces[k]
    ax.fill(s[:,0],s[:,1],color='#263947',ec='#173446',lw=.8,zorder=9)
    ax.axhline(0,color='#273e48',lw=2,zorder=10)
    if box:ax.add_patch(Rectangle((-.125,.05),.75,.30,fill=False,ls='--',ec='#888',lw=1.1,zorder=10))
    ax.set(xlim=(-.15,.8),ylim=(-.015,.375),xlabel='x (m)',ylabel='Height above ground (m)')
    ax.set_aspect('equal')
    ax.set_title(titles[k]+f"\n$C_L$ = {results[k]['cl']:.4f}, $C_D$ = {results[k]['cd']:.5f}")
def centers(x):return .25*(x[:-1,:-1]+x[1:,:-1]+x[1:,1:]+x[:-1,1:])
def tri(f,tag):
    x,y=[centers(f[tag+'_'+v]) for v in ['CoordinateX','CoordinateY']]
    ids=np.arange(x.size).reshape(x.shape)
    aa,bb,cc,dd=[v.ravel() for v in [ids[:-1,:-1],ids[1:,:-1],ids[1:,1:],ids[:-1,1:]]]
    t=np.r_[np.c_[aa,bb,cc],np.c_[aa,cc,dd]]
    if tag=='near':
        aa,bb,cc,dd=[v.ravel() for v in [ids[-1,:-1],ids[0,:-1],ids[0,1:],ids[-1,1:]]]
        t=np.r_[t,np.c_[aa,bb,cc],np.c_[aa,cc,dd]]
    valid=np.isin(f[tag+'_Iblank'],[1,-1]).ravel()
    t=t[valid[t].all(axis=1)]
    tr=mtri.Triangulation(x.ravel(),y.ravel(),t)
    return tr
tris=[{tag:tri(f,tag) for tag in ['background','near']} for f in fields]

fig,ax=plt.subplots(figsize=(10,4))
for k,s in enumerate(surfaces):ax.plot(s[:,0],s[:,1],c=colors[k],lw=2,label=titles[k])
ax.add_patch(Rectangle((-.125,.05),.75,.30,fill=False,ls='--',ec='#777',label='Allowed surface box'))
ax.axhline(0,c='#273e48',lw=2,label='Slip ground')
ax.set(xlim=(-.15,.8),ylim=(-.015,.37),xlabel='x (m)',ylabel='Height above ground (m)',
       title='Shape and position: before / after optimization')
ax.set_aspect('equal');ax.legend(loc='upper right',fontsize=9);ax.grid(alpha=.15)
save(fig,'geometry_comparison')


fig,axes=plt.subplots(1,2,figsize=(11,3.8),constrained_layout=True)
for k,s in enumerate(surfaces):
    axes[0].plot(s[:,0],s[:,1],c=colors[k],lw=1.6,label=titles[k])
base=surfaces[0];ch=(base[0,:2]+base[272,:2])/2-base[136,:2]
unit=ch/np.linalg.norm(ch)
material_x=(base[:273,:2]-base[136,:2])@unit/.5
delta=(surfaces[1][:273,1]-base[:273,1])*1000
for sl,ls,label in [(slice(0,137),'-','Upper surface'),(slice(136,273),'--','Lower surface')]:
    order=np.argsort(material_x[sl]);axes[1].plot(material_x[sl][order],delta[sl][order],ls,lw=1.6,label=label)
axes[0].set(xlabel='x (m)',ylabel='Height (m)',title='Airfoil outline, close view');axes[0].set_aspect('equal');axes[0].legend(loc='upper left',bbox_to_anchor=(0,1.75),ncol=2,frameon=False,fontsize=9)
axes[1].axhline(0,c='#999',lw=.7)
axes[1].set(xlabel='Initial chord station $x/c_0$',ylabel='Change in surface height (mm)',title='Deformation, including position and incidence')
axes[1].legend()
for ax in axes:ax.grid(alpha=.2)
save(fig,'shape_change')

fig,axes=plt.subplots(1,2,figsize=(13,4.2),constrained_layout=True)
for k,(f,ax) in enumerate(zip(fields,axes)):
    for tag,stride,col in [('background',4,'#adc1cf'),('near',3,'#316a88')]:
        x,y=[f[tag+'_'+v] for v in ['CoordinateX','CoordinateY']]
        seg=[np.c_[x[i,:],y[i,:]] for i in range(0,x.shape[0],stride)]
        seg += [np.c_[x[:,j],y[:,j]] for j in range(0,x.shape[1],stride)]
        ax.add_collection(LineCollection(seg,colors=col,linewidths=.3,rasterized=True))
    body(ax,k,True)
fig.suptitle('Structured overset grid — selected grid lines; dashed box is a surface constraint')
save(fig,'mesh_comparison')

fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
for k,s in enumerate(cp):
    # Slice points come as ordered surface segments, with duplicated endpoints.
    # Material upper/lower side is identified in chord coordinates.
    for side,ls in [(s[:,4]>=0,'-'),(s[:,4]<0,'--')]:
        section=s[side];order=np.argsort(section[:,3])
        axes[0].plot(section[order,3],section[order,6],ls,c=colors[k],lw=1.6)
    axes[0].plot([],[],c=colors[k],label=titles[k])
    axes[1].scatter(s[:,0],s[:,6],s=5,c=colors[k],label=titles[k],alpha=.8)
axes[0].plot([],[],'k-',lw=1.2,label='Upper surface')
axes[0].plot([],[],'k--',lw=1.2,label='Lower surface')
axes[0].set(xlabel='$x_{chord}/c_0$',ylabel='$C_p$',title='Surface pressure in chord coordinates')
axes[1].set(xlabel='Global x (m)',ylabel='$C_p$',title='Surface pressure at actual position')
for ax in axes:ax.invert_yaxis();ax.grid(alpha=.2);ax.legend(fontsize=9)
save(fig,'cp_comparison')

for name,kind,lo,hi,cmap,title in [
 ('cp_contour_comparison','CoefPressure',-5,1,'coolwarm','Static pressure coefficient $C_p$'),
 ('total_pressure_comparison','TotalPressure_Pa',.65,1.02,'viridis','Total pressure $p_0/p_{0,\\infty}$'),
 ('pressure_loss_comparison','loss',-.03,1.,'magma','Total pressure loss $(p_{0,\\infty}-p_0)/q_\\infty$')]:
    if kind=='CoefPressure':
        lows=[float(f[tag+'_CoefPressure'][np.isin(f[tag+'_Iblank'],[1,-1])].min()) for f in fields for tag in ['near','background']]
        lo=np.floor(min(lows)*2)/2
    elif kind=='TotalPressure_Pa':
        ratios=[f[tag+'_TotalPressure_Pa'][np.isin(f[tag+'_Iblank'],[1,-1])]/float(f['total_pressure_inf_Pa']) for f in fields for tag in ['near','background']]
        lo=np.floor(min(float(v.min()) for v in ratios)*100)/100
        hi=np.ceil(max(float(v.max()) for v in ratios)*1000)/1000
    fig,axes=plt.subplots(1,2,figsize=(13,4.1),constrained_layout=True)
    norm=Normalize(lo,hi);levels=np.linspace(lo,hi,65)
    for k,(f,ax) in enumerate(zip(fields,axes)):
        for tag in ['background','near']:
            if kind=='loss':values=(float(f['total_pressure_inf_Pa'])-f[tag+'_TotalPressure_Pa'])/float(f['dynamic_pressure_inf_Pa'])
            elif kind=='TotalPressure_Pa':values=f[tag+'_'+kind]/float(f['total_pressure_inf_Pa'])
            else:values=f[tag+'_'+kind]
            ax.tricontourf(tris[k][tag],values.ravel(),levels=levels,cmap=cmap,norm=norm,extend='both',zorder=2)
        body(ax,k)
    fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes,shrink=.85,pad=.015,label=title,extend='both')
    fig.suptitle(title+' — same color scale before / after')
    save(fig,name)

fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
e=np.array([r['evaluation'] for r in history]);cl=np.array([r['cl'] for r in history]);cd=np.array([r['cd'] for r in history])
feas=np.array([r['feasible'] for r in history]);failed=np.array([r['cfd_failed'] for r in history])
axes[0,0].plot(e,-cl,'.-',ms=3,lw=.8,c=colors[0]);axes[0,0].scatter(e[feas],-cl[feas],s=18,c=colors[1],label='Feasible CFD evaluation')
axes[0,0].scatter(e[failed],-cl[failed],s=35,c='#b42525',marker='x',label='CFD failed; rejected')
axes[0,0].set(ylabel='Downforce coefficient $-C_L$');axes[0,0].legend(fontsize=8)
axes[0,1].plot(e,cd,'.-',ms=3,lw=.8);axes[0,1].scatter(e[failed],cd[failed],s=35,c='#b42525',marker='x');axes[0,1].axhline(settings['drag_cap'],ls='--',c='r',label='Drag limit');axes[0,1].set(ylabel='$C_D$');axes[0,1].legend()
axes[1,0].plot(e,[r['chord_incidence_deg'] for r in history],'.-');axes[1,0].set(ylabel='Actual chord incidence (deg)')
axes[1,1].semilogy(e,[r['residual_ratio_initial'] for r in history],'.-');axes[1,1].set(ylabel='Final / reference initial CFD residual')
for ax in axes.flat:ax.set_xlabel('CFD evaluation (baseline and all trial points)');ax.grid(alpha=.2)
fig.suptitle('Optimization evaluation history')
save(fig,'optimization_history')

fig,axes=plt.subplots(1,3,figsize=(11,3.6),constrained_layout=True)
for ax,key,title in zip(axes,['downforce_N_per_m','drag_N_per_m','downforce_to_drag'],['Downforce (N/m span)','Drag (N/m span)','Downforce / drag']):
    vals=[r[key] for r in results]
    ax.bar(['Before','After'],vals,color=colors,width=.6)
    for i,v in enumerate(vals):ax.text(i,v,f'{v:.2f}',ha='center',va='bottom')
    ax.set_ylim(0,max(vals)*1.16);ax.set_title(title)
fig.suptitle('Wing loads at 200 km/h; fixed reference chord 0.5 m')
save(fig,'performance_comparison')

psqp_log=root/'PSQP.out'
if psqp_log.exists():
    pattern=r'NIT=\s*(\d+)\s+NFV=\s*(\d+)\s+NFG=\s*(\d+)\s+F=\s*([+\-.\dEeDd]+)\s+C=\s*([+\-.\dEeDd]+)\s+G=\s*([+\-.\dEeDd]+)'
    iterations=np.array([[float(v.replace('D','E')) for v in row] for row in re.findall(pattern,psqp_log.read_text(errors='replace'))])
    if iterations.size:
        np.savetxt(root/'optimizer_iterations.csv',iterations,delimiter=',',header='major_iteration,function_evaluations,gradient_evaluations,CL,max_constraint_violation,lagrangian_gradient_log_precision',comments='')
        fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
        axes[0].plot(iterations[:,0],-iterations[:,3],'.-',c=colors[0]);axes[0].set(ylabel='Downforce coefficient $-C_L$')
        valid=iterations[:,0]>0
        axes[1].semilogy(iterations[valid,0],iterations[valid,5],'.-');axes[1].axhline(settings['optimizer_options'].get('TOLG',1e-4),c='r',ls='--',label='Gradient target');axes[1].set(ylabel='Lagrangian gradient norm');axes[1].legend(fontsize=8)
        axes[2].semilogy(iterations[:,0],np.maximum(iterations[:,4],1e-12),'.-');axes[2].axhline(settings['optimizer_options'].get('TOLC',1e-6),c='r',ls='--',label='Feasibility target');axes[2].set(ylabel='Maximum constraint violation');axes[2].legend(fontsize=8)
        for ax in axes:ax.set_xlabel('PSQP major iteration');ax.grid(alpha=.2)
        fig.suptitle('Optimizer convergence (norms at printed log precision)')
        save(fig,'optimizer_convergence')

keys=[('cl','CL, upward positive'),('cd','CD'),('cmz','CMz about fixed initial quarter chord'),
 ('downforce_N_per_m','Downforce, N/m span'),('drag_N_per_m','Drag, N/m span'),
 ('downforce_to_drag','Downforce / drag'),('force_down_N_extruded_span','Downforce on 0.05 m extrusion, N'),
 ('chord_incidence_deg','Actual chord incidence, deg'),('incidence_variable_deg','Incidence design variable, deg'),
 ('actual_chord_m','Actual LE-to-TE chord, m'),('minimum_clearance_m','Minimum surface ground clearance, m'),
 ('maximum_section_thickness_m','Maximum section thickness, m'),('maximum_height_m','Maximum surface height, m'),
 ('minimum_x_m','Minimum surface x, m'),('maximum_x_m','Maximum surface x, m'),
 ('box_minimum_margin_m','Minimum surface box margin, m'),('surface_yplus_max','Maximum airfoil y+'),
 ('surface_cp_min','Minimum surface Cp'),('surface_cp_max','Maximum surface Cp'),
 ('maximum_active_mach','Maximum active-cell Mach number'),
 ('minimum_total_pressure_ratio','Minimum active-cell total pressure / freestream'),
 ('maximum_total_pressure_ratio','Maximum active-cell total pressure / freestream'),
 ('residual_ratio_initial','Final/reference initial residual ratio')]
lines=['OPTIMIZATION RESULTS: SLIP GROUND',
       'PROVISIONAL: optimizer still running.' if not (root/'run_complete.json').exists() else 'Optimizer exit: '+str(json.loads((root/'run_complete.json').read_text())['inform']),'',f'After design folder: {a.after}',
       'CL is upward-positive; more negative CL means more downforce.',
       'All force coefficients use the fixed initial reference chord 0.5 m and area 0.025 m^2.',
       'Ground is an ideal slip/symmetry boundary; road boundary-layer effects are excluded.','',
       f"{'Quantity':49s} {'Before':>15s} {'After':>15s} {'Change':>15s}",'-'*100]
for key,label in keys:
    x,y=[r.get(key,float('nan')) for r in results]
    lines.append(f'{label:49s} {x:15.8g} {y:15.8g} {y-x:+15.8g}')
improvement=(results[1]['downforce_N_per_m']/results[0]['downforce_N_per_m']-1)*100
lines+=['',f'Downforce increase: {improvement:.4f}%',f"Drag cap: {settings['drag_cap']:.10g}",
 'Box: -0.125 <= x <= 0.625 m; 0.05 <= y <= 0.35 m.',
 'Maximum section thickness must be >= 0.03 m (6% of fixed initial chord).',
 'Incidence has no explicit lower/upper bound. No fixed leading/trailing-edge equalities.',
 f"After feasible according to stored constraint checks: {results[1]['feasible']}",
 f"After CFD failure flag: {results[1]['cfd_failed']}",
 f'CFD evaluations stored: {len(history)} (baseline and all evaluated trial/recheck points).']
exitfile=root/'run_complete.json'
if exitfile.exists():
    run=json.loads(exitfile.read_text());lines += [f"Optimizer termination: {run['inform']}",
      f"Optimizer wall time: {run.get('optimizer_wall_time_s',run['total_wall_time_s']):.1f} s ({run.get('optimizer_wall_time_s',run['total_wall_time_s'])/60:.2f} min), excluding baseline setup and earlier attempts"]
else:lines+=['Status: optimization in progress; these figures are provisional.']
if (root/'optimizer_iterations.csv').exists():
    it=np.loadtxt(root/'optimizer_iterations.csv',delimiter=',',skiprows=1,ndmin=2)
    lines += [f'PSQP major iterations in printed log: {int(it[:,0].max())}',
              f'Last printed Lagrangian gradient norm: {it[-1,5]:.6g} (target {settings["optimizer_options"].get("TOLG",1e-4):.6g})']
lines += ['', 'MODEL AND NUMERICAL LIMITATIONS',
 'Single quasi-2D section, fully turbulent SA RANS, 88,548 cells before overset blanking.',
 'Shape variables move paired span planes together in global y; 20 shape variables + 1 incidence variable.',
 'The geometric incidence variable rotates the initial FFD; local shape changes can also change chord incidence.',
 'Surface box applies to every CFD boundary vertex; linear connecting segments stay inside this convex box.',
 'Thickness uses maximum chord-normal separation at paired material stations and is piecewise differentiable.',
 'The entire overset grid is warped together with frozen interpolation connectivity/weights; adjoints match this model.',
 'No grid-independence study has yet been completed. These are preliminary optimization results.',
 'A returned optimizer status is not evidence of a global optimum.','']
comparison={'before':results[0],'after':results[1],'downforce_increase_percent':improvement,'after_folder':a.after}
check=root/'full_connectivity_validation/final/result.json'
if check.exists():
    v=json.loads(check.read_text());comparison['independent_full_connectivity_validation']=v
    if v.get('accepted_convergence'):
        validated_gain=(v['downforce_N_per_m']/results[0]['downforce_N_per_m']-1)*100
        comparison['full_connectivity_downforce_increase_percent']=validated_gain
        lines += ['INDEPENDENT FULL-CONNECTIVITY CFD CHECK',
          'Same selected geometry; full hole cutting and interpolation rebuilt; saved optimized flow used as restart.',
          f"CL={v['cl']:.10f}; CD={v['cd']:.10f}; downforce={v['downforce_N_per_m']:.6f} N/m span",
          f"Downforce increase relative to initial NACA: {validated_gain:.4f}%",
          f"Converged: {v['accepted_convergence']}; feasible: {v['feasible']}; residual ratio={v['residual_ratio_initial']:.6g}",
          f"Downforce difference from frozen-connectivity prediction: {(v['downforce_N_per_m']/results[1]['downforce_N_per_m']-1)*100:.5f}%",
          'This is an interpolation/connectivity check on the same deformed mesh, not grid independence.','']
    else:lines += ['INDEPENDENT FULL-CONNECTIVITY CHECK DID NOT CONVERGE; see full_connectivity_validation/.','']
(root/'comparison_results.txt').write_text('\n'.join(lines))
(root/'comparison_results.json').write_text(json.dumps(comparison,indent=2)+'\n')
np.savetxt(root/'history_summary.csv',np.array([[r['evaluation'],r['cl'],r['cd'],r['chord_incidence_deg'],
    r['minimum_clearance_m'],r['maximum_section_thickness_m'],r['feasible'],r['cfd_failed'],r['elapsed_wall_time_s']] for r in history]),
    delimiter=',',header='evaluation,CL,CD,chord_incidence_deg,clearance_m,max_thickness_m,feasible,cfd_failed,elapsed_s',comments='')
print('\n'.join(lines[:35]))

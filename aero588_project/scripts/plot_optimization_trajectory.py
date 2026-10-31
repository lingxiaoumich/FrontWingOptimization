"""Plot the iteration history of an FFD optimization run.

Figure 1 (control_point_history): every FFD control point at every PSQP major
iteration, light (early) to dark (late), with the allowed surface box.
Figure 2 (cl_cd_history): CL*c and CD*c histories (L'/q and D'/q in m, so the
choice of reference chord drops out), the drag polar against the drag cap, and
the signed suction-surface skin friction of the saved baseline/final flows,
used to judge whether the run moved into separated flow.

The lattice is rebuilt without pyGeo: rotate the baseline lattice about the
incidence pivot by -(incidence + 3) deg, then add the shape displacements
(the same order pyGeo uses; checked against ffd_geometry/control_points.npz).
Only the baseline and final flow fields are archived, so the separation check
covers those two designs; intermediate designs are judged from CL/CD alone.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root
out=root/'plots';out.mkdir(exist_ok=True)
settings=json.loads((root/'settings.json').read_text())
history=[json.loads(x) for x in (root/'history.jsonl').read_text().splitlines() if x.strip()]
iters=np.loadtxt(root/'optimizer_iterations.csv',delimiter=',',skiprows=1)
box=settings['box_m'];c0=settings['reference_chord_m'];cap=settings['drag_cap']

# Optimizer iterates: one evaluation per major iteration, after the pre-optimization
# baseline and incidence finite-difference checks. Matched to the PSQP log by CL.
nmajor=int(iters[:,0].max())+1
start=len(history)-1-nmajor
run=history[start:start+nmajor]
assert np.allclose([h['cl'] for h in run],iters[:nmajor,3],atol=1e-5), 'history/PSQP iteration mismatch'
k=np.arange(nmajor)
# L'/q and D'/q (m): coefficients times the chord they were referenced to
clc=c0*np.array([h['cl'] for h in run]);cdc=c0*np.array([h['cd'] for h in run]);capc=c0*cap

ramp=['#86b6ef','#6da7ec','#5598e7','#3987e5','#2a78d6','#256abf','#1c5cab','#184f95','#104281','#0d366b']
cmap=LinearSegmentedColormap.from_list('early_to_late',ramp)
norm=Normalize(0,nmajor-1);colors=cmap(norm(k))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'figure.facecolor':'white','savefig.facecolor':'white','font.family':'DejaVu Sans'})
ink,muted,grid='#2b2b2b','#6b6b6b','#e4e4e4'
def save(fig,name):
    fig.savefig(out/(name+'.png'),dpi=190,bbox_inches='tight')
    fig.savefig(out/(name+'.svg'),bbox_inches='tight');plt.close(fig)

# ---------------------------------------------------------------- control points
lattice0=np.load(root/'ffd_geometry/control_points.npz')['baseline'][:,:,0,:2]
xy=np.loadtxt(root/'source/airfoil_coordinates.csv',delimiter=',',skiprows=1)
pivot=np.array([.125,xy[136,1]+.125*np.sin(np.deg2rad(3.))])
def lattice(design):
    t=np.deg2rad(-(design['incidence'][0]+3.))
    R=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    shape=np.stack([np.reshape(design.get('shape_x',np.zeros(lattice0[...,0].size)),lattice0.shape[:2]),
                    np.reshape(design['shape_y'],lattice0.shape[:2])],-1)
    return (lattice0-pivot)@R.T+pivot+shape
P=np.array([lattice(h['design']) for h in run])  # (iteration, station, lower/upper, xy)
saved=np.load(root/'ffd_geometry/control_points.npz')['final'][:,:,0,:2]
assert np.abs(lattice(json.loads((root/'final/design.json').read_text()))-saved).max()<1e-12
nst=P.shape[1]
surf=[np.loadtxt(root/l/'surface_geometry.csv',delimiter=',',skiprows=1) for l in ['baseline','final']]

fig=plt.figure(figsize=(13,9.6))
gs=fig.add_gridspec(3,nst,height_ratios=[2.3,1,1],hspace=.42,wspace=.12)
ax=fig.add_subplot(gs[0,:])
ax.add_patch(Rectangle((box['xmin'],box['ymin']),box['xmax']-box['xmin'],box['ymax']-box['ymin'],
    fill=False,ls='--',ec=muted,lw=1.1))
ax.text(box['xmax']-.005,box['ymax']-.008,'Allowed surface box',ha='right',va='top',color=muted,fontsize=9)
ax.axhline(0,c=ink,lw=2);ax.text(box['xmax']-.005,.004,'Slip ground',ha='right',va='bottom',color=muted,fontsize=9)
ax.fill(surf[0][:,0],surf[0][:,1],color='#dcdcdc',lw=0,label='Baseline airfoil')
ax.plot(surf[1][:,0],surf[1][:,1],c=ink,lw=1,label='Final airfoil')
for j,lat in [(0,P[0]),(1,P[-1])]:
    st=dict(c='#b5b5b5' if j==0 else '#555555',lw=.7,ls=':' if j==0 else '-')
    for s in range(2):ax.plot(lat[:,s,0],lat[:,s,1],**st)
    for i in range(nst):ax.plot(lat[i,:,0],lat[i,:,1],**st)
for i in range(nst):
    for s in range(2):ax.scatter(P[:,i,s,0],P[:,i,s,1],c=colors,s=9,lw=0,zorder=5)
names=[[f'L{i+1}',f'U{i+1}'] for i in range(nst)]
for i in range(nst):
    ax.annotate(names[i][1],P[-1,i,1],xytext=(0,7),textcoords='offset points',ha='center',fontsize=8,color=ink)
    ax.annotate(names[i][0],P[-1,i,0],xytext=(0,-12),textcoords='offset points',ha='center',fontsize=8,color=ink)
ax.plot([],[],c='#b5b5b5',lw=.7,ls=':',label='Baseline FFD lattice')
ax.plot([],[],c='#555555',lw=.7,label='Final FFD lattice')
ax.set(xlim=(box['xmin']-.02,box['xmax']+.02),ylim=(-.012,box['ymax']+.012),xlabel='x (m)',ylabel='Height above ground (m)')
ax.set_aspect('equal');ax.legend(loc='upper left',bbox_to_anchor=(.03,.9),fontsize=8.5,frameon=False)
ax.set_title(f'FFD control-point history, {nmajor} PSQP major iterations '
             f'(design constraints act on the airfoil surface, not on the control points)',fontsize=10.5)

d=(P-P[:1])*1000  # displacement from the baseline position, mm (shape + incidence rotation)
lim=np.abs(d).max()*1.12
for row,s in [(1,1),(2,0)]:
    for i in range(nst):
        z=fig.add_subplot(gs[row,i])
        z.axhline(0,c=grid,lw=.8,zorder=0);z.axvline(0,c=grid,lw=.8,zorder=0)
        z.plot(d[:,i,s,0],d[:,i,s,1],c='#c8c8c8',lw=.6,zorder=1)
        z.scatter(d[:,i,s,0],d[:,i,s,1],c=colors,s=12,lw=0,zorder=2)
        z.set(xlim=(-lim,lim),ylim=(-lim,lim));z.set_aspect('equal')
        z.set_title(f'{names[i][s]}  ({P[0,i,s,0]:.2f}, {P[0,i,s,1]:.3f}) m',fontsize=8.5)
        z.tick_params(labelsize=7.5)
        if i:z.set_yticklabels([])
        else:z.set_ylabel(('Upper' if s else 'Lower')+' points\nΔy (mm)',fontsize=8.5)
        if row==2:z.set_xlabel('Δx (mm)',fontsize=8.5)
        else:z.set_xticklabels([])
cb=fig.colorbar(ScalarMappable(norm,cmap),ax=fig.axes[1:],shrink=.85,pad=.015,aspect=30)
cb.set_label('PSQP major iteration (light = early)')
fig.text(.125,.035,'Lower panels: displacement of each control point from its baseline position, same scale in every panel; '
         'it includes the rigid incidence rotation\nabout the pivot (0.125, 0.080) m, which alone moves the trailing-edge points by about '
         f'{abs(np.deg2rad(run[-1]["incidence_variable_deg"]+3))*.375*1000:.1f} mm. Box: {box["xmin"]} ≤ x ≤ {box["xmax"]} m, '
         f'{box["ymin"]} ≤ y ≤ {box["ymax"]} m on every surface vertex.',fontsize=8,color=muted)
save(fig,'control_point_history')

# ---------------------------------------------------------------- separation check
def suction_side(label):
    """Signed skin friction on each surface face: |Cf| with the sign of the
    near-wall tangential velocity relative to the boundary-layer edge velocity."""
    f=np.load(root/label/'fields.npz')
    s=np.loadtxt(root/label/'surface_cp.csv',delimiter=',',skiprows=1)
    chord=json.loads((root/label/'result.json').read_text())['actual_chord_m']
    X,Y=f['near_CoordinateX'],f['near_CoordinateY']
    assert np.allclose(X[:-1,0],s[:,0],atol=1e-6)
    t=np.c_[np.diff(X[:,0]),np.diff(Y[:,0])];t/=np.linalg.norm(t,axis=1)[:,None]
    uinf=float(f['mach_inf'])*np.sqrt(float(f['gamma']))
    ut=(f['near_VelocityX']*t[:,:1]+f['near_VelocityY']*t[:,1:])/uinf
    cx=.25*(X[:-1,:-1]+X[1:,:-1]+X[:-1,1:]+X[1:,1:]);cy=.25*(Y[:-1,:-1]+Y[1:,:-1]+Y[:-1,1:]+Y[1:,1:])
    wall=.5*(np.c_[X[:-1,0],Y[:-1,0]]+np.c_[X[1:,0],Y[1:,0]])
    dist=np.hypot(cx-wall[:,:1],cy-wall[:,1:])
    edge=ut[np.arange(len(ut)),np.argmin(abs(dist-.012),axis=1)]  # 12 mm: outside the boundary layer
    reversed_=np.sign(ut)!=np.sign(edge)[:,None]
    cf=.5*(s[:,7]+np.roll(s[:,7],-1))*np.where(reversed_[:,0],-1,1)
    xc=.5*(s[:,3]+np.roll(s[:,3],-1))*c0/chord
    lower=slice(136,272)  # leading edge to lower trailing-edge corner; 272+ is the blunt base
    xl,cfl=xc[lower],cf[lower]
    rev=np.where((cfl<0)&(xl>.1))[0]
    first=rev[0];onset=np.interp(0,[cfl[first],cfl[first-1]],[xl[first],xl[first-1]])
    contiguous=bool(np.all(np.diff(rev)==1))
    # height of the reversed layer: first cell off the wall whose flow agrees with the edge
    height=[dist[c,np.argmin(reversed_[c])] for c in range(136,272) if reversed_[c,0]]
    return dict(x=xl,cf=cfl,separation_onset_x_over_c=float(onset),separated_fraction_of_chord=float(1-onset),
                single_reversed_region=contiguous,reversed_until_x_over_c=float(xl[rev[-1]]),
                max_reversed_layer_height_mm=float(max(height)*1000),
                min_signed_cf=float(cfl.min()))
sep={l:suction_side(l) for l in ['baseline','final']}

fig,axs=plt.subplots(2,2,figsize=(12,8.4),constrained_layout=True)
ax=axs[0,0];ax.plot(k,clc,c='#c8c8c8',lw=.8,zorder=1);ax.scatter(k,clc,c=colors,s=16,lw=0,zorder=2)
ax.set(xlabel='PSQP major iteration',ylabel="$C_L c = L'/q_\\infty$ (m, upward positive)",
       title='(a) Lift: more negative = more downforce')
ax=axs[0,1];ax.plot(k,cdc,c='#c8c8c8',lw=.8,zorder=1);ax.scatter(k,cdc,c=colors,s=16,lw=0,zorder=2)
ax.set(xlabel='PSQP major iteration',ylabel="$C_D c = D'/q_\\infty$ (m)",title=f'(b) Drag (cap {capc:.3f} m is off this scale)')
ax=axs[1,0]
ax.axhline(capc,c=muted,ls='--',lw=1)
ax.text(-clc.max()+.003,capc-.0003,f'$C_D c$ cap = {capc:.3f} m (never active)',va='top',color=muted,fontsize=9)
for ld in [40,50,60]:
    xs=np.array([.55,.775]);ax.plot(xs,xs/ld,c=grid,lw=.8,zorder=0)
    ax.text(xs[1],xs[1]/ld,f'  $-L/D$ = {ld}',va='center',fontsize=8,color=muted)
ax.plot(-clc,cdc,c='#c8c8c8',lw=.8,zorder=1);ax.scatter(-clc,cdc,c=colors,s=16,lw=0,zorder=2)
ax.set(xlim=(.56,.84),ylim=(.009,.0265),xlabel='Downforce $-C_L c$ (m)',ylabel='$C_D c$ (m)',
       title='(c) Drag polar of the iterates')
ax=axs[1,1];ax.axhline(0,c=ink,lw=.8)
for l,c,lab in [('baseline',ramp[0],'Baseline (iteration 0)'),('final',ramp[-1],f'Final (iteration {nmajor-1})')]:
    r=sep[l];m=r['x']>.3
    ax.plot(r['x'][m],r['cf'][m]*1000,c=c,lw=2,label=lab+f": reversed from x/c = {r['separation_onset_x_over_c']:.3f}")
    ax.axvspan(r['separation_onset_x_over_c'],1,color=c,alpha=.12,lw=0)
ax.set(xlim=(.3,1.005),xlabel='x/c along the deformed chord',ylabel='Signed $C_f$ × 1000',
       title='(d) Suction-side skin friction ($C_f$ < 0: reversed flow)\n'
             'only the baseline and final flows are archived')
ax.legend(loc='upper right',fontsize=8.5,frameon=False)
inset=ax.inset_axes([.3,.3,.45,.38])
for l,c in [('baseline',ramp[0]),('final',ramp[-1])]:
    r=sep[l];inset.plot(r['x'],r['cf']*1000,c=c,lw=1.6)
inset.axhline(0,c=ink,lw=.6);inset.set(xlim=(.85,1.0),ylim=(-.5,.6));inset.tick_params(labelsize=7.5)
inset.set_title('Trailing-edge zoom',fontsize=8);inset.grid(c=grid,lw=.6)
for ax in axs.flat:ax.grid(c=grid,lw=.6);ax.set_axisbelow(True)
cb=fig.colorbar(ScalarMappable(norm,cmap),ax=axs,shrink=.6,pad=.01,aspect=30)
cb.set_label('PSQP major iteration (light = early)')
save(fig,'cl_cd_history')

slope=np.diff(cdc)/np.diff(-clc)
summary=dict(major_iterations=nmajor,history_evaluations_used=[start,start+nmajor-1],
    cl_c_m=[float(clc[0]),float(clc[-1])],cd_c_m=[float(cdc[0]),float(cdc[-1])],
    cd_c_cap_m=float(capc),max_cd_c_over_cap=float(cdc.max()/capc),
    monotonic_cl_c=bool(np.all(np.diff(clc)<0)),monotonic_cd_c=bool(np.all(np.diff(cdc)>0)),
    dcd_c_d_downforce_first_last_5=[float(slope[:5].mean()),float(slope[-5:].mean())],
    max_control_point_displacement_mm=float(np.linalg.norm(d,axis=-1).max()),
    separation={l:{q:v for q,v in r.items() if q not in ('x','cf')} for l,r in sep.items()},
    separation_method='sign of first-cell tangential velocity vs. velocity 12 mm off the wall; '
                      'lower surface only, blunt trailing-edge base excluded')
(out/'history_separation_check.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))

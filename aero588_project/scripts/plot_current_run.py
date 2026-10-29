"""Plot the final baseline from the extracted CGNS fields and surface CSV.

Run from any directory: python3 scripts/plot_current_run.py
Requires numpy and matplotlib. The cached arrays in reports/current_run_plots
come from baseline_a0_000_vol.cgns; ghost/rind cells have been removed.
"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable

ROOT = Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--case-directory',type=Path,help='Plot a case made by run_ground_case.py')
args=parser.parse_args()
IS_NEW_CASE=args.case_directory is not None
if IS_NEW_CASE:
    ROOT=args.case_directory.resolve()
OUT = ROOT / ('reports/plots' if IS_NEW_CASE else 'reports/current_run_plots')
FIELDS = np.load(OUT / 'baseline_fields.npz')
RESULT_REL='results/cfd_result.json' if IS_NEW_CASE else 'results/baseline_a0.json'
SURFACE_REL='results/surface.csv' if IS_NEW_CASE else 'results/baseline_surface.csv'
REPORT = json.loads((ROOT / RESULT_REL).read_text())
GROUND_LABEL = {'slip':'Ideal slip ground', 'moving_no_slip':'Moving no-slip ground',
                'stationary_no_slip':'Stationary no-slip ground'}.get(
                    REPORT.get('ground_model','slip'),'Ground')
SURFACE = np.genfromtxt(ROOT / SURFACE_REL, delimiter=',', skip_header=1)
C = REPORT['chord_m']
H = REPORT['clearance_m']
INCIDENCE = REPORT['alpha_deg']
SPEED = REPORT['speed_ms']*3.6
HEIGHT = float(SURFACE[:,1].max())
GAMMA = float(FIELDS['gamma'])
PINF = float(FIELDS['pressure_reference_Pa'])
PTINF = float(FIELDS['total_pressure_inf_Pa'])
QINF = float(FIELDS['dynamic_pressure_inf_Pa'])
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11,
                     'axes.titlesize':12, 'axes.labelsize':11,
                     'axes.spines.top':False, 'axes.spines.right':False,
                     'savefig.facecolor':'white', 'figure.facecolor':'white'})
NAVY, ORANGE = '#175676', '#c65c25'
SUBTITLE = f'NACA 0012 | c = {C:.2f} m | clearance = {H:.2f} m | {SPEED:g} km/h | incidence = {INCIDENCE:g}°'

def save(fig, name):
    fig.savefig(OUT / (name+'.png'), dpi=230, bbox_inches='tight')
    fig.savefig(OUT / (name+'.svg'), bbox_inches='tight')
    plt.close(fig)

def cell_centers(a):
    return .25*(a[:-1,:-1]+a[1:,:-1]+a[1:,1:]+a[:-1,1:])

def axes_geometry(ax, bounds, label=True):
    ax.fill(SURFACE[:,0], SURFACE[:,1], facecolor='#263947', edgecolor='#172732', lw=.8, zorder=6)
    ax.axhline(0,color='#263947',lw=2,zorder=7)
    ax.set(xlim=bounds[:2], ylim=bounds[2:], xlabel='x (m)', ylabel='Height above ground (m)')
    ax.set_aspect('equal')
    if label:
        ax.annotate(GROUND_LABEL, xy=(bounds[0]+.025,0), xytext=(0,7),
                    textcoords='offset points', fontsize=9, color='#263947', zorder=8)

def draw_grid(ax, tag, stride, color, lw):
    x,y=(FIELDS[tag+'_'+k] for k in ['CoordinateX','CoordinateY'])
    indices_i=sorted(set(range(0,x.shape[0],stride))|{x.shape[0]-1})
    indices_j=sorted(set(range(0,x.shape[1],stride))|{x.shape[1]-1})
    segments=[np.column_stack((x[i,:],y[i,:])) for i in indices_i]
    segments += [np.column_stack((x[:,j],y[:,j])) for j in indices_j]
    ax.add_collection(LineCollection(segments,colors=color,linewidths=lw,rasterized=True))

def triangulation(tag):
    x,y=(cell_centers(FIELDS[tag+'_'+k]) for k in ['CoordinateX','CoordinateY'])
    ids=np.arange(x.size).reshape(x.shape)
    a,b,c,d=[v.ravel() for v in [ids[:-1,:-1],ids[1:,:-1],ids[1:,1:],ids[:-1,1:]]]
    triangles=np.concatenate([np.column_stack([a,b,c]),np.column_stack([a,c,d])])
    if tag=='near':
        # Connect the periodic circumferential seam of the closed near grid.
        a,b,c,d=[v.ravel() for v in [ids[-1,:-1],ids[0,:-1],ids[0,1:],ids[-1,1:]]]
        triangles=np.concatenate([triangles,np.column_stack([a,b,c]),np.column_stack([a,c,d])])
    valid=np.isin(FIELDS[tag+'_Iblank'],[1,-1]).ravel()
    tri=mtri.Triangulation(x.ravel(),y.ravel(),triangles)
    tri.set_mask(~valid[triangles].all(axis=1))
    return tri

TRIS={tag:triangulation(tag) for tag in ['background','near']}

def contour(ax, values, levels, cmap, norm):
    for tag in ['background','near']:
        cs=ax.tricontourf(TRIS[tag],values[tag].ravel(),levels=levels,cmap=cmap,
                         norm=norm,extend='both',zorder=2)
        cs.set_rasterized(True)

def colorbar(fig, axes, cmap, norm, label, ticks=None):
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes,shrink=.9,pad=.02,extend='both')
    cb.set_label(label)
    if ticks is not None: cb.set_ticks(ticks)

def mesh_plot():
    fig,axs=plt.subplots(2,1,figsize=(11,9),layout='constrained',height_ratios=[1.15,1])
    fig.suptitle('Structured overset mesh\n'+SUBTITLE,fontsize=15)
    for ax,near_stride,bg_stride,bounds,title in [
        (axs[0],3,4,(-.3*C,1.45*C,-.025*C,HEIGHT+.5*C),'Airfoil and ground: selected grid lines'),
        (axs[1],1,1,(-.05*C,.28*C,H-.025*C,HEIGHT+.04*C),'Leading-edge detail: every grid line')]:
        draw_grid(ax,'background',bg_stride,'#a6b7c0',.4)
        draw_grid(ax,'near',near_stride,NAVY,.38)
        axes_geometry(ax,bounds)
        ax.set_title(title,loc='left')
    xlow=float(SURFACE[np.argmin(SURFACE[:,1]),0])
    axs[0].annotate('',xy=(xlow,H),xytext=(xlow,0),arrowprops={'arrowstyle':'<->','color':ORANGE,'lw':1.4})
    axs[0].text(xlow+.02*C,.48*H,f'{H:.2f} m',fontsize=10,color=ORANGE,rotation=90,va='center')
    cells=REPORT.get('total_cells_before_blanking',77093)
    spacing=REPORT.get('first_wall_spacing_m',5e-6)*1e6
    fig.text(.5,-.07,f'{cells:,} total cells before blanking • one spanwise cell • first wall spacing {spacing:g} µm\n'
             'Near-field grid: blue; background: grey. Raw overlapping grids shown.',ha='center',fontsize=10,color='#455660')
    save(fig,'mesh')

def pressure_plot():
    fig,axs=plt.subplots(2,1,figsize=(12,8.5),layout='constrained',height_ratios=[1.8,1])
    fig.suptitle('Total pressure and wake loss\n'+SUBTITLE,fontsize=15)
    vals={tag:FIELDS[tag+'_TotalPressure_Pa']/1000 for tag in TRIS}
    values=np.concatenate([vals[t][np.isin(FIELDS[t+'_Iblank'],[1,-1])] for t in vals])
    levels=np.linspace(np.floor(values.min()*20)/20,np.ceil(values.max()*20)/20,81)
    norm=Normalize(levels[0],levels[-1])
    contour(axs[0],vals,levels,'viridis',norm)
    axes_geometry(axs[0],(-.3*C,2.4*C,-.02*C,HEIGHT+.55*C))
    axs[0].set_title('Absolute local total pressure',loc='left')
    colorbar(fig,axs[0],'viridis',norm,'Total pressure, $p_0$ (kPa)',np.linspace(levels[0],levels[-1],6))
    loss={tag:(PTINF-FIELDS[tag+'_TotalPressure_Pa'])/QINF for tag in TRIS}
    loss_min,loss_max=(-.08,1.1) if IS_NEW_CASE else (-.05,.5)
    loss_levels=np.linspace(loss_min,loss_max,81)
    loss_norm=TwoSlopeNorm(vmin=loss_min,vcenter=0,vmax=loss_max)
    contour(axs[1],loss,loss_levels,'RdBu_r',loss_norm)
    te_y=float(SURFACE[np.argmax(SURFACE[:,3]),1])
    axes_geometry(axs[1],(.7*C,2.45*C,max(0,te_y-.13*C),te_y+.22*C),label=False)
    axs[1].set_title('Wake close-up: positive values indicate total-pressure loss',loc='left')
    loss_ticks=[-.08,0,.2,.4,.6,.8,1.] if IS_NEW_CASE else [-.05,0,.1,.2,.3,.4,.5]
    colorbar(fig,axs[1],'RdBu_r',loss_norm,'$(p_{0,∞}-p_0)/q_∞$',loss_ticks)
    fig.text(.5,-.07,rf'$p_0=p[1+(\gamma-1)M^2/2]^{{\gamma/(\gamma-1)}}$, '
             rf'$\gamma=1.4$; $p_{{0,∞}}={PTINF/1000:.3f}$ kPa; $q_∞={QINF:.1f}$ Pa.'
             f'\nCell-center fields; blanked/flooded cells excluded. Wake scale saturates above {loss_max:.2f}; pressure overshoots retained.',
             ha='center',fontsize=10,color='#455660')
    save(fig,'total_pressure')

def cp_plot():
    le=int(np.argmin(SURFACE[:,3]))
    lower_te=next(i for i in range(le+1,len(SURFACE)) if np.isclose(SURFACE[i,3],1.,atol=2e-6))
    upper=SURFACE[:le+1][::-1]
    lower=SURFACE[le:lower_te+1]
    fig,axs=plt.subplots(1,2,figsize=(12,5.2),layout='constrained',width_ratios=[1.6,1])
    fig.suptitle('Surface pressure coefficient\n'+SUBTITLE,fontsize=15)
    for ax in axs:
        ax.plot(upper[:,3],upper[:,6],color=NAVY,lw=1.8,label='Upper surface')
        ax.plot(lower[:,3],lower[:,6],color=ORANGE,lw=1.8,label='Lower surface (ground-facing)')
        ax.axhline(0,color='#8d9aa1',lw=.8,ls='--')
        ax.set(xlabel='Chordwise position / c',ylabel='$C_p$',
               ylim=(max(1.12,SURFACE[:,6].max()+.1),min(-.65,SURFACE[:,6].min()-.15)))
        ax.grid(alpha=.17)
    axs[0].set(xlim=(0,1),title='Full chord')
    axs[0].legend(loc='upper right',fontsize=9)
    axs[1].set(xlim=(0,.12),title='Leading-edge close-up')
    axs[0].text(.035,.035,rf'$C_L={REPORT["cl_upward_positive"]:.5f}$   '
                rf'$C_D={REPORT["cd"]:.5f}$',transform=axs[0].transAxes,
                bbox={'facecolor':'white','edgecolor':'#d9e0e4','pad':5},fontsize=10)
    fig.text(.5,-.12,'Raw saved surface data, without smoothing. Cp axis inverted by convention.\n'
             'The finite trailing-edge face is excluded from the upper/lower curves.',ha='center',fontsize=10,color='#455660')
    save(fig,'surface_cp')

def static_cp_plot():
    fig,ax=plt.subplots(figsize=(12,5),layout='constrained')
    fig.suptitle('Static-pressure coefficient in the flow field\n'+SUBTITLE,fontsize=15)
    vals={tag:FIELDS[tag+'_CoefPressure'] for tag in TRIS}
    values=np.concatenate([vals[t][np.isin(FIELDS[t+'_Iblank'],[1,-1])] for t in vals])
    cpmin=float(np.floor(values.min()*5)/5)
    cpmax=float(np.ceil(values.max()*5)/5)
    levels=np.linspace(cpmin,cpmax,81)
    norm=TwoSlopeNorm(vmin=cpmin,vcenter=0,vmax=cpmax)
    contour(ax,vals,levels,'RdBu_r',norm)
    axes_geometry(ax,(-.3*C,2.1*C,-.02*C,HEIGHT+.55*C))
    colorbar(fig,ax,'RdBu_r',norm,'Static-pressure coefficient, $C_p$')
    fig.text(.5,-.09,f'Negative Cp: pressure below freestream. {GROUND_LABEL}; saved CFD volume solution.',ha='center',fontsize=10,color='#455660')
    save(fig,'static_pressure_cp')

if __name__=='__main__':
    mesh_plot()
    pressure_plot()
    cp_plot()
    static_cp_plot()
    provenance={}
    volume_rel=f'runs/{REPORT["case_id"]}_000_vol.cgns' if IS_NEW_CASE else 'runs/baseline_a0/baseline_a0_000_vol.cgns'
    for relative in [volume_rel,SURFACE_REL]:
        p=ROOT/relative
        provenance[relative]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
    (OUT/'source_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Plots saved to',OUT)

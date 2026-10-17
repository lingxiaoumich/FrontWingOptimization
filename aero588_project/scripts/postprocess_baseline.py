import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path('/project')
rows=[]
history=root/'reports/baseline_a0_converged_history.log'
if not history.exists():
    history=root/'reports/baseline_a0.log'
for line in history.read_text().splitlines():
    cols=line.split()
    if len(cols)==12 and all(v.isdigit() for v in cols[:3]):
        try:
            rows.append([float(cols[2]),float(cols[-3]),float(cols[-2]),float(cols[-1])])
        except ValueError:
            pass
data=np.array(rows)
report=json.loads((root/'reports/baseline_a0.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
axes[0].semilogy(data[:,0],data[:,3]/report['residual_initial'],color='#245a81')
axes[0].axhline(1e-8,ls='--',color='#b14c40',label='Target')
axes[0].set(xlabel='Solver work iterations (continued run)',ylabel='Residual / original initial residual',title='Baseline convergence')
axes[0].legend()
axes[1].plot(data[:,0],-data[:,1],label='Downforce coefficient')
axes[1].plot(data[:,0],data[:,2],label='Drag coefficient')
axes[1].set(xlabel='Solver work iterations (continued run)',ylabel='Coefficient',title='Airfoil forces')
axes[1].legend()
fig.savefig(root/'reports/baseline_convergence.png',dpi=180)
tail=data[-min(5,len(data)):]
report['last_five_samples_cl_range']=float(np.ptp(tail[:,1]))
report['last_five_samples_cd_range']=float(np.ptp(tail[:,2]))
report['overset_orphan_counts']=[int(v) for v in re.findall(r'Total number of orphans:\s*(\d+)',(root/'reports/baseline_a0.log').read_text())]
slice_path=root/'runs/baseline_a0/baseline_a0_000_slices.dat'
if slice_path.exists():
    points=[]
    for line in slice_path.read_text().splitlines():
        values=line.split()
        if len(values)==9:
            try:
                points.append([float(v) for v in values])
            except ValueError:
                pass
    points=np.array(points)
    report['surface_yplus_max']=float(points[:,8].max())
    report['surface_yplus_median']=float(np.median(points[:,8]))
    np.savetxt(root/'reports/baseline_surface.csv',points,delimiter=',',header='x_m,y_m,z_m,x_over_c,y_over_c,z_over_c,Cp,Cf_magnitude,yplus',comments='')
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    mid=(points[:,1].max()+points[:,1].min())/2
    for mask,label in [(points[:,1]>=mid,'Upper surface'),(points[:,1]<mid,'Lower surface')]:
        values=points[mask]
        values=values[np.argsort(values[:,0])]
        ax.plot(values[:,0]/0.7,values[:,6],label=label)
    ax.invert_yaxis()
    ax.set(xlabel='x / chord',ylabel='Pressure coefficient',title='Zero-incidence NACA 0012: ideal slip ground')
    ax.legend()
    fig.savefig(root/'reports/baseline_pressure.png',dpi=180)
(root/'reports/baseline_a0.json').write_text(json.dumps(report,indent=2)+'\n')
status='CONVERGED' if report['accepted_convergence'] else 'NOT CONVERGED'
text=f'''# NACA 0012 ground-effect baseline

Status: {status}. Zero-incidence baseline, not an optimized design.

- Chord: 0.70 m; speed: 100 km/h; lowest-point clearance: 0.20 m.
- Ground: ideal slip/symmetry plane. Road shear is omitted; this is not a moving no-slip road simulation.
- Air: 288.15 K and 101325 Pa; Reynolds number {report['reynolds']:.0f}; Mach {report['mach']:.5f}.
- Fully turbulent Spalart-Allmaras RANS. One spanwise cell with symmetry boundaries models a two-dimensional section.
- Standard NACA 0012 finite trailing edge, 1.764 mm thick. The earlier sharp trailing-edge coordinate file is not used for this CFD mesh.
- 77093 cells before overset blanking; no nonpositive corner Jacobians or orphan overset cells found.
- Maximum wall y+: {report.get('surface_yplus_max',float('nan')):.3f}; median wall y+: {report.get('surface_yplus_median',float('nan')):.3f}.

## Computed baseline quantities

Use these only within the stated model and convergence status.

| Quantity | Value |
|---|---:|
| Lift coefficient, upward positive | {report['cl_upward_positive']:.6f} |
| Drag coefficient | {report['cd']:.6f} |
| Downforce per unit span | {report['downforce_N_per_m']:.3f} N/m |
| Drag per unit span | {report['drag_N_per_m']:.3f} N/m |
| Downforce / drag | {report['downforce_to_drag']:.3f} |
| Residual / original initial residual | {report['residual_reduction']:.3e} |

## Limitations and next checks

Mesh independence is not established. The ground has no viscous road boundary layer. Transition, finite-span effects, wheels, and the vehicle body are excluded. These are preliminary isolated-airfoil results. Verify grid sensitivity and ground treatment before using them for racing-car design decisions or optimization claims.

Angle-of-attack bounds are -4 to +4 degrees. For a ground-effect angle study, rotate the airfoil geometry while keeping the freestream horizontal and enforcing the 0.20 m minimum clearance; changing freestream angle alone would create a physically inconsistent road setup.

## Reproduction

Use mdolab/public:u24-gcc-ompi-stable with /home/mdolabuser/.bashrc_mdolab sourced. Mount a local working project at /project. Run scripts/build_ground_mesh.py, then run scripts/run_baseline.py with four MPI ranks. Add --restart to continue from runs/baseline_a0/baseline_a0_000_vol.cgns. Run scripts/plot_mesh.py and scripts/postprocess_baseline.py for the plots. Copy deliverables to G:\\My Drive\\aero588_project.
'''
(root/'reports/baseline_report.md').write_text(text)

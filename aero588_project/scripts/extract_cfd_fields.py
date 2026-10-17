"""Extract plot arrays and surface Cp from one completed ground-effect case."""
from pathlib import Path
import argparse
import json
import re
import numpy as np
from adf_reader import ADFReader

parser=argparse.ArgumentParser()
parser.add_argument('--case-directory',type=Path,required=True)
parser.add_argument('--adf-library',type=Path,required=True)
args=parser.parse_args()
root=args.case_directory.resolve()
cfg=json.loads((root/'case_config.json').read_text())
prefix=cfg['case_id']
volume=root/'runs'/f'{prefix}_000_vol.cgns'
tree=ADFReader(args.adf_library).load(volume)
base=next(t for t in tree['children'] if t['label']=='CGNSBase_t')
def child(t,name):return next(c for c in t['children'] if c['name']==name)
ref=child(base,'ReferenceState')
pref=float(child(ref,'Pressure')['data'][0])
minf=float(child(ref,'Mach')['data'][0])
gamma=float(child(child(child(base,'FlowEquationSet'),'GasModel'),'SpecificHeatRatio')['data'][0])
ptinf=pref*(1+(gamma-1)/2*minf**2)**(gamma/(gamma-1))
qinf=.5*gamma*pref*minf**2
output={'pressure_reference_Pa':np.array(pref),'mach_inf':np.array(minf),
        'gamma':np.array(gamma),'total_pressure_inf_Pa':np.array(ptinf),
        'dynamic_pressure_inf_Pa':np.array(qinf)}
stats={}
for zone in (c for c in base['children'] if c['label']=='Zone_t'):
    tag='near' if 'airfoil' in zone['name'] else 'background'
    coords=child(zone,'GridCoordinates')
    sol=next(c for c in zone['children'] if c['label']=='FlowSolution_t')
    rind=child(sol,'Rind')['data']
    # CGNS Rind_t specifies low/high rind for each index dimension.
    core=tuple(slice(int(rind[2*i]),-int(rind[2*i+1]) or None) for i in range(3))
    project=(lambda a:a[:,0,:]) if tag=='near' else (lambda a:a[:,:,0])
    for name in ['CoordinateX','CoordinateY']:
        output[tag+'_'+name]=project(child(coords,name)['data'])
    for name in ['Pressure','Density','VelocityX','VelocityY','Mach','CoefPressure','Iblank']:
        output[tag+'_'+name]=project(child(sol,name)['data'][core]).astype(float)
    p,m,cp,ib=[output[tag+'_'+n] for n in ['Pressure','Mach','CoefPressure','Iblank']]
    valid=np.isin(ib,[1,-1])
    ptotal=pref*p*(1+(gamma-1)/2*m**2)**(gamma/(gamma-1))
    output[tag+'_TotalPressure_Pa']=ptotal
    cp_error=float(np.max(np.abs((p[valid]-1)*pref/qinf-cp[valid])))
    calculated_mach=np.hypot(output[tag+'_VelocityX'],output[tag+'_VelocityY'])/np.sqrt(gamma*p/output[tag+'_Density'])
    mach_error=float(np.max(np.abs(calculated_mach[valid]-m[valid])))
    assert cp_error<5e-5,(tag,cp_error)
    assert mach_error<1e-6,(tag,mach_error)
    assert np.isfinite(ptotal[valid]).all()
    assert (p[valid]>0).all()
    stats[tag]={'total_pressure_Pa_range':[float(ptotal[valid].min()),float(ptotal[valid].max())],
                'Cp_reconstruction_max_error':cp_error,'Mach_reconstruction_max_error':mach_error,
                'iblank_counts':{str(int(k)):int(v) for k,v in zip(*np.unique(ib,return_counts=True))}}
out=root/'reports/plots'
out.mkdir(exist_ok=True)
np.savez_compressed(out/'baseline_fields.npz',**output)
(out/'field_validation.json').write_text(json.dumps({'pressure_reference_Pa':pref,
     'total_pressure_inf_Pa':ptinf,'qinf_Pa':qinf,'zones':stats},indent=2)+'\n')

rows=[]
for line in (root/'runs'/f'{prefix}_000_slices.dat').read_text().splitlines():
    values=line.split()
    if len(values)==9:
        try: rows.append([float(v) for v in values])
        except ValueError: pass
points=np.array(rows)
assert points.shape[1]==9 and len(points)>200
mesh_info=json.loads((root/'reports/mesh_check.json').read_text())
le=np.array(mesh_info['leading_edge_m'])
angle=np.deg2rad(cfg['geometry_rotation_deg'])
delta=points[:,:2]-le
points[:,3]=(delta[:,0]*np.cos(angle)+delta[:,1]*np.sin(angle))/cfg['chord_m']
points[:,4]=(-delta[:,0]*np.sin(angle)+delta[:,1]*np.cos(angle))/cfg['chord_m']
np.savetxt(root/'results/surface.csv',points,delimiter=',',
           header='x_m,y_m,z_m,chord_x_over_c,chord_y_over_c,z_over_c,Cp,Cf_magnitude,yplus',comments='')
result_path=root/'results/cfd_result.json'
result=json.loads(result_path.read_text())
result['surface_yplus_max']=float(points[:,8].max())
result['surface_yplus_median']=float(np.median(points[:,8]))
result['overset_orphan_counts']=[int(v) for v in re.findall(r'Total number of orphans:\s*(\d+)',(root/'reports/cfd.log').read_text())]
result['total_cells_before_blanking']=mesh_info['total_cells_before_blanking']
assert result['overset_orphan_counts'] and max(result['overset_orphan_counts'])==0
result_path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'field_validation':stats,'surface_yplus_max':result['surface_yplus_max']},indent=2))

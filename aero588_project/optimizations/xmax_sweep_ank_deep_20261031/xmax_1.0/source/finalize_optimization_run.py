"""Save resolved settings/provenance and verify the selected optimization design."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from airfoil_geometry import maximum_chord_normal_thickness
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root
settings=json.loads((root/'settings.json').read_text())
settings['sensitivity_overset_mode']='frozen, consistent with primal USMesh warping'
meshcheck=json.loads((root/'source/mesh_check.json').read_text())
le=meshcheck['leading_edge_m'];te=meshcheck['trailing_edge_midpoint_m']
fields_path=root/'baseline/fields.npz'
if fields_path.exists():
    fields=np.load(fields_path)
    settings['mesh_definition']={
      'units':'m','nearfield_vertices_ijk':[289,2,81],'nearfield_cells':23040,
      'background_vertices_ijk':[319,207,2],'background_cells':65508,
      'total_cells_before_blanking':88548,'initial_compute_cells':84781,'initial_fringe_cells':1036,
      'initial_domain_x_m':[float(fields['background_CoordinateX'].min()),float(fields['background_CoordinateX'].max())],
      'initial_domain_y_m':[float(fields['background_CoordinateY'].min()),float(fields['background_CoordinateY'].max())],
      'span_m':[0.,.05],'spanwise_cells':1,'ground_first_cell_height_m':.0025,
      'airfoil_first_wall_spacing_m':2e-6,'nearfield_march_distance_m':.025,
      'airfoil_boundary':'stationary viscous no-slip wall',
      'ground_boundary':'symmetry / ideal slip at y=0',
      'span_boundaries':'z symmetry at z=0 and z=0.05 m',
      'outer_boundaries':'farfield on inlet, outlet and top',
      'deformation':'IDWarp USMesh applied to the whole overset grid; airfoil is the design surface'}
shape_groups=settings.get('shape_groups',['shape'])
settings['design_variable_units']={**{name:'m' for name in shape_groups},'incidence':'deg'}
settings['optimizer_variable_scales']={**{name:100. for name in shape_groups},'incidence':1.}
settings['geometry_constraints_on']='289 vertices of the actual CFD section surface, including closure; tied span planes have identical x-y coordinates.'
settings['moment_reference_point_m']=[.125,le[1]+.125*np.sin(np.deg2rad(3.)),.025]
if settings['optimizer'].endswith('PSQP'):
    settings['optimizer_all_options']={'XMAX':1e16,'TOLX':1e-16,'TOLC':1e-6,'TOLG':1e-6,'RPF':.0001,'MIT':1000,'MFV':2000,'MET':2,'MEC':2,'IPRINT':2,'IOUT':6,'IFILE':'PSQP.out'}
    settings['optimizer_all_options'].update(settings['optimizer_options'])
    settings['numerical_step_limit_note']='XMAX applies to optimizer steps in scaled coordinates, not to final design bounds. Incidence remains unbounded.'
settings['force_family']='airfoil_only: all stationary no-slip wing surface faces'
settings['ground_force_in_objective']=False
settings['stored_feasibility_tolerance']={'surface_box_nondimensional':1e-6,
    'thickness_nondimensional':1e-6,'relative_drag_cap':1e-5,
    'note':'Physical inequalities unchanged; tolerances only classify finite precision results.'}
(root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
log=(root/'run.log').read_text(errors='replace')
# The complete resolved defaults/overrides are printed by ADflow and IDWarp.
start=log.index("{'ADPC':")
end=log.index("  -> Reading CGNS File:",start)
resolved=log[start:end]
software=(root/'source/software_versions.json').read_text()
(root/'settings.txt').write_text('CFD AND OPTIMIZER SETTINGS\n\n'+json.dumps(settings,indent=2)
    +'\n\nSOFTWARE VERSIONS AND SOURCE COMMITS\n'+software
    +'\nALL RESOLVED ADFLOW AND IDWARP OPTIONS (FROM RUN LOG)\n'+resolved
    +'\nPrimal and adjoint both use frozen overset connectivity while USMesh warps the entire grid.\n')
orphans=[int(v) for v in re.findall(r'Total number of orphans:\s*(\d+)',log)]
assert orphans and max(orphans)==0,orphans
assert 'Full overset update derivatives not implemented' not in log
assert 'ERROR STOP' not in log
summary={}
for label in ['baseline','final','best_feasible']:
    path=root/label/'result.json'
    if not path.exists():continue
    r=json.loads(path.read_text());xy=np.loadtxt(root/label/'surface_geometry.csv',delimiter=',',skiprows=1)[:,:2]
    # Detect proper intersections of nonadjacent boundary segments.
    aa=xy[:-1];bb=xy[1:]
    def cross(u,v):return u[...,0]*v[...,1]-u[...,1]*v[...,0]
    ab=(bb-aa)[:,None,:];cd=(bb-aa)[None,:,:]
    o1=cross(ab,aa[None,:,:]-aa[:,None,:]);o2=cross(ab,bb[None,:,:]-aa[:,None,:])
    o3=cross(cd,aa[:,None,:]-aa[None,:,:]);o4=cross(cd,bb[:,None,:]-aa[None,:,:])
    intersections=int(np.count_nonzero(np.triu((o1*o2 < -1e-24)&(o3*o4 < -1e-24),k=1)))
    assert intersections==0,(label,intersections)
    tol=5e-7
    if label in ['baseline','best_feasible']:
        assert xy[:,0].min()>=-.125-tol and xy[:,0].max()<=.625+tol
        assert xy[:,1].min()>=.05-tol and xy[:,1].max()<=.35+tol
        if settings.get('shape_directions')=='xy':
            thickness=maximum_chord_normal_thickness(xy)
        else:
            chord=(xy[0]+xy[272])/2-xy[136]
            normal=np.array([-chord[1],chord[0]])/np.linalg.norm(chord)
            thickness=np.max((xy[:137]-xy[272:135:-1])@normal)
        assert thickness>=.03-tol
        assert r['cd']/settings['drag_cap']<=1+1e-5
    assert not r['cfd_failed'],(label,'CFD failure')
    assert r['residual_ratio_initial']<=1e-10,(label,'CFD residual')
    grids=r.get('grid_validation',{})
    assert grids and all(v['nonpositive_corner_areas']==0 for v in grids.values())
    assert abs(grids['background']['minimum_grid_y_m'])<1e-12
    summary[label]={'self_intersections':intersections,'cfd_converged':True,
        'feasible':r['feasible'],'grid_corner_checks_passed':True,
        'minimum_surface_box_margin_m':r['box_minimum_margin_m'],
        'max_yplus':r.get('surface_yplus_max')}
validation={'solutions':summary,'overset_updates_checked':len(orphans),'max_orphans':max(orphans),
    'geometry_derivative_check':json.loads((root/'source/geometry_derivative_check.json').read_text()),
    'aerodynamic_derivative_check':json.loads((root/'source/aerodynamic_derivative_check.json').read_text())}
extended=root/'source/geometry_derivative_check_extended.json'
if extended.exists():validation['extended_geometry_derivative_check']=json.loads(extended.read_text())
validation['optimizer_terminated']=(root/'run_complete.json').exists()
if validation['optimizer_terminated']:validation['optimizer_exit']=json.loads((root/'run_complete.json').read_text())['inform']
(root/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
manifest={}
for path in (root/'source').iterdir():
    if path.is_file() and path.name!='sha256_manifest.json':manifest[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
(root/'source/sha256_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(validation,indent=2))

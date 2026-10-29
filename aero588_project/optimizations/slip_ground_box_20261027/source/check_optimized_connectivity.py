"""Independent full-connectivity CFD check, starting from saved optimized states.

Stage the original source grid/geometry and selected design.json under source/;
use the selected solution as source/restart_vol.cgns. No optimizer or adjoint is
run. The executed optimization driver defines the identical geometry and flow.
"""
from pathlib import Path
# Reuse the exact setup rather than maintaining a second geometry model.
setup = (Path(__file__).with_name('run_slip_optimization.py')).read_text()
setup = setup.split('# Baseline data only uses')[0]
setup = setup.replace("'oversetUpdateMode':'frozen'", "'oversetUpdateMode':'full'")
exec(compile(setup, 'shared_optimization_setup', 'exec'))
settings.update({'description':'Independent CFD validation with fully rebuilt overset connectivity',
    'optimizer':None,'optimizer_options':{},'sensitivity_overset_mode':'not used; no adjoints evaluated',
    'drag_cap':json.loads((root/'source/replay_settings.json').read_text())['drag_cap'],
    'flow_derived':{'density_kg_m3':float(ap.rho),'viscosity_Pa_s':float(ap.mu),
                  'mach':float(ap.mach),'reynolds':float(ap.rho*V*c0/ap.mu)}})
dump(root/'settings.json',settings)
if comm.rank==0:
    (root/'settings.txt').write_text('INDEPENDENT FULL-CONNECTIVITY CFD VALIDATION\n\n'+json.dumps(settings,indent=2)+'\n')
target={k:np.array(v) for k,v in json.loads((root/'source/design_to_check.json').read_text()).items()}
geo.setDesignVars(target)
pts=geo.update('surface_check')
solver(ap,writeSolution=False)
failed=bool(ap.solveFailed or ap.fatalFail)
ri,rs,rf=map(float,solver.getResNorms())
r={'cl':None,'cd':None,'cmz':None,'cfd_failed':failed,
   'residual_initial':ri,'residual_start':rs,'residual_final':rf,
   'residual_ratio_initial':rf/ri if ri else None,'design':realdict(target),
   'incidence_variable_deg':float(target['incidence'][0]),**metrics(pts),
   'elapsed_wall_time_s':time.perf_counter()-start}
if not ap.fatalFail:
    f={};solver.evalFunctions(ap,f)
    cl,cd,cmz=[float(f['wing_'+name+'_wing']) for name in ['cl','cd','cmz']]
    q=.5*ap.rho*V**2
    r.update({'cl':cl,'cd':cd,'cmz':cmz,'downforce_N_per_m':float(-cl*q*c0),
        'drag_N_per_m':float(cd*q*c0),'downforce_to_drag':float(-cl/cd),
        'force_down_N_extruded_span':float(-cl*q*c0*span)})
    solver.writeSolution(outputDir=str(root/'final'),baseName='flow',number=0)
r['accepted_convergence']=bool(not failed and r['residual_ratio_initial']<=1e-10)
r['feasible']=bool(r['accepted_convergence'] and r['cd']<=settings['drag_cap']*(1+1e-5)
    and r['minimum_x_m']>=-.125-5e-7 and r['maximum_x_m']<=.625+5e-7
    and r['minimum_clearance_m']>=.05-5e-7 and r['maximum_height_m']<=.35+5e-7
    and r['maximum_section_thickness_m']>=.03-5e-7)
if comm.rank==0:
    np.savetxt(root/'final/surface_geometry.csv',pts,delimiter=',',header='x_m,y_m,z_m',comments='')
dump(root/'final/design.json',realdict(target))
dump(root/'final/result.json',r)
dump(root/'validation_complete.json',r)
if comm.rank==0:print('FULL_CONNECTIVITY_VALIDATION',json.dumps(r),flush=True)

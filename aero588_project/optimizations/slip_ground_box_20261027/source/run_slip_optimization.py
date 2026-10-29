"""Adjoint shape/incidence optimization of the project's slip-ground airfoil.

Execute under mpirun in mdolab/public:u24-gcc-ompi-stable. All dimensional
references are fixed to the initial chord. The bounds apply to surface mesh
vertices (and thus every segment of the piecewise-linear CFD surface).
"""
import argparse
import copy
from dataclasses import asdict, is_dataclass
import json
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI
from pygeo import DVGeometry
from pyspline import Curve
from idwarp import USMesh
from adflow import ADFLOW
from baseclasses import AeroProblem
from pyoptsparse import Optimization, OPT

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, required=True)
p.add_argument('--geometry-only', action='store_true')
p.add_argument('--maxit', type=int, default=50)
p.add_argument('--restart-design', type=Path)
p.add_argument('--hotstart', type=Path)
p.add_argument('--reuse-checks', action='store_true')
p.add_argument('--evaluate-design', type=Path)
p.add_argument('--optimizer',choices=['SLSQP','PSQP'],default='PSQP')
a = p.parse_args()
root = a.root.resolve()
comm = MPI.COMM_WORLD
start = time.perf_counter()
for d in ['source', 'baseline', 'final', 'best_feasible', 'evaluations', 'plots']:
    if comm.rank == 0:
        (root/d).mkdir(parents=True, exist_ok=True)
comm.Barrier()

def json_default(value):
    if is_dataclass(value):return asdict(value)
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,Path):return str(value)
    raise TypeError(f'Cannot serialize {type(value).__name__}')

def dump(path, value):
    if comm.rank == 0:
        path.write_text(json.dumps(value, indent=2, default=json_default)+'\n')

def realdict(d):
    return {k: np.real(v).tolist() for k,v in d.items()}

c0, span, V, T, P = .5, .05, 200/3.6, 288.15, 101325.
xy0 = np.loadtxt(root/'source/airfoil_coordinates.csv', delimiter=',', skiprows=1)
angle0 = np.deg2rad(3.)
# Invert the baseline rigid rotation to identify the original NACA coordinates.
xlocal = (1+np.cos(np.linspace(0,2*np.pi,273)))/2*c0
ylocal = .6*(.2969*np.sqrt(xlocal/c0)-.126*xlocal/c0-.3516*(xlocal/c0)**2
             +.2843*(xlocal/c0)**3-.1015*(xlocal/c0)**4)*c0
ylocal = np.where(np.arange(273)<=136, ylocal, -ylocal)
offset = float(xy0[136,1]+.25*c0*np.sin(angle0))
pivot = np.array([.25*c0, offset, 0.])

# Ten chordwise stations; two heights; two tied span planes = 40 FFD points,
# 20 independent shape variables. FFD enclosure is not the regulation box.
ffd = np.empty((10,2,2,3))
for i,x in enumerate(np.linspace(-.01,c0+.01,10)):
    for j,y in enumerate([-.06,.06]):
        for k,z in enumerate([-.005,span+.005]):
            ffd[i,j,k] = [np.cos(angle0)*(x-.25*c0)-np.sin(angle0)*y+.25*c0,
                          np.sin(angle0)*(x-.25*c0)+np.cos(angle0)*y+offset,z]
if comm.rank == 0:
    with (root/'source/ffd.xyz').open('w') as f:
        f.write('1\n10 2 2\n')
        for axis in range(3):
            np.savetxt(f, ffd[...,axis].ravel(order='F'), fmt='%.16e')
comm.Barrier()
geo = DVGeometry(str(root/'source/ffd.xyz'))
geo.addRefAxis('incidence_axis', curve=Curve(k=2,t=[0.,0.,1.,1.],coef=np.array([
    [pivot[0],pivot[1],-.005],[pivot[0],pivot[1],span+.005]])), axis=np.array([0.,0.,1.]))
def incidence(val, g):
    g.rot_z['incidence_axis'].coef[:] = -(val[0]+3.)
geo.addGlobalDV('incidence', [-3.], incidence, scale=1.)
idx = geo.getLocalIndex(0)
shapes = [{int(idx[i,j,k]):np.array([0.,1.,0.]) for k in range(2)}
          for i in range(10) for j in range(2)]
geo.addShapeFunctionDV('shape', shapes, scale=100.)
points = np.c_[xy0,np.full(len(xy0),span/2)]
geo.addPointSet(points, 'surface_check')
xinitial = copy.deepcopy(geo.getValues())

# Geometry functions are evaluated directly on the actual CFD boundary, not FFD.
# Thickness is the maximum chord-normal separation of corresponding upper/lower
# NACA material stations. Piecewise differentiable at a change of peak station.
# Closed TE strip is included in the box but not in paired section thickness.
def geometry_values(pts):
    xy=pts[:,:2]
    le=xy[136]
    te=(xy[0]+xy[272])/2
    chord=te-le
    normal=np.array([-chord[1],chord[0]])/np.sqrt(np.dot(chord,chord))
    separation=np.sum((xy[:137]-xy[272:135:-1])*normal,axis=1)
    peak=int(np.argmax(np.real(separation)))
    return {'surface_x':xy[:,0]/c0, 'surface_y':xy[:,1]/c0,
            'max_thickness':np.array([separation[peak]/c0])}

def metrics(pts):
    vals=geometry_values(pts)
    xy=pts[:,:2]; chord=(xy[0]+xy[272])/2-xy[136]
    return {'minimum_clearance_m':float(np.min(xy[:,1])),
            'maximum_height_m':float(np.max(xy[:,1])),
            'minimum_x_m':float(np.min(xy[:,0])), 'maximum_x_m':float(np.max(xy[:,0])),
            'maximum_section_thickness_m':float(vals['max_thickness'][0]*c0),
            'maximum_section_thickness_over_reference_chord':float(vals['max_thickness'][0]),
            'actual_chord_m':float(np.linalg.norm(chord)),
            'chord_incidence_deg':float(-np.rad2deg(np.arctan2(chord[1],chord[0]))),
            'box_minimum_margin_m':float(min(xy[:,0].min()+.125,.625-xy[:,0].max(),
                                            xy[:,1].min()-.05,.35-xy[:,1].max()))}

def geometry_sens(pts):
    # Analytic pyGeo total Jacobian combined with complex-step derivatives of
    # the inexpensive coordinate functions; no CFD finite-difference gradients.
    out={}
    for name in ['surface_x','surface_y']:
        axis=0 if name=='surface_x' else 1
        jac=np.zeros((len(pts),len(pts),3))
        jac[np.arange(len(pts)),np.arange(len(pts)),axis]=1/c0
        out[name]=geo.totalSensitivity(jac,'surface_check')
    dp=np.zeros_like(pts)
    for i in range(len(pts)):
        for j in range(2):
            z=pts.astype(complex);z[i,j]+=1e-30j
            dp[i,j]=np.imag(geometry_values(z)['max_thickness'][0])/1e-30
    out['max_thickness']=geo.totalSensitivity(dp,'surface_check')
    return out

# Verify geometry derivatives independently at a small nonzero perturbation.
xcheck=copy.deepcopy(xinitial)
xcheck['shape']=np.linspace(-1e-5,1e-5,20)
xcheck['incidence']=np.array([-3.001])
geo.setDesignVars(xcheck);pc=geo.update('surface_check')
gs=geometry_sens(pc)
checks={}
for key,step in [('incidence',1e-4),('shape',1e-6)]:
    direction=np.ones(len(xcheck[key]));direction/=np.linalg.norm(direction)
    xp=copy.deepcopy(xcheck);xm=copy.deepcopy(xcheck)
    xp[key]+=step*direction;xm[key]-=step*direction
    geo.setDesignVars(xp);vp=geometry_values(geo.update('surface_check'))
    geo.setDesignVars(xm);vm=geometry_values(geo.update('surface_check'))
    for name in gs:
        fd=(vp[name]-vm[name])/(2*step)
        exact=np.asarray(gs[name][key])@direction
        checks[name+'_'+key]=float(np.max(np.abs(fd-exact)))
geo.setDesignVars(xinitial)
dump(root/'source/geometry_derivative_check.json',checks)
assert max(checks.values())<1e-6,checks
if a.geometry_only:
    if comm.rank==0: print('GEOMETRY_CHECK',checks,metrics(geo.update('surface_check')),flush=True)
    raise SystemExit

options={
    'gridFile':str(root/'source/ground.cgns'),'outputDirectory':str(root/'evaluations'),
    'equationType':'RANS','equationMode':'steady','liftIndex':2,'turbulenceModel':'SA',
    'turbulenceOrder':'first order','discretization':'central plus scalar dissipation',
    'useWallFunctions':False,'useANKSolver':True,'useNKSolver':True,
    'ANKSwitchTol':10.,'NKSwitchTol':1e-4,'ANKCoupledSwitchTol':1e-2,
    'ANKSecondOrdSwitchTol':1e-2,'ANKADPC':True,'NKADPC':True,'NKViscPC':True,
    'NKSubspaceSize':500,'NKPCILUFill':5,'NKLinearSolveTol':.1,'NKUseEW':False,
    'turbResScale':10000.,'nCycles':4000,'L2Convergence':1e-10,'L2ConvergenceRel':1e-10,
    'MGCycle':'sg','useZipperMesh':False,'nearWallDist':.01,'oversetUpdateMode':'frozen',
    'designSurfaceFamily':'airfoil','meshSurfaceFamily':'airfoil',
    'monitorVariables':['resrho','resturb','cl','cd'],
    'volumeVariables':['cp','mach','blank'],'surfaceVariables':['cp','cf','yplus'],
    'writeTecplotSurfaceSolution':True,'writeVolumeSolution':True,'writeSurfaceSolution':True,
    'adjointL2Convergence':1e-7,'adjointMaxIter':2000,'adjointSubspaceSize':500,
    'ADPC':True,'viscPC':True,'ILUFill':4,'frozenTurbulence':False,
}
if (root/'source/restart_vol.cgns').exists():
    options['restartFile']=str(root/'source/restart_vol.cgns')
meshoptions={'gridFile':options['gridFile'],
    'symmetryPlanes':[[[0.,0.,0.],[0.,1.,0.]],[[0.,0.,0.],[0.,0.,1.]],
                      [[0.,0.,span],[0.,0.,1.]]]}
optimizeroptions=({'ACC':1e-5,'MAXIT':a.maxit,'IPRINT':2,'IFILE':str(root/'SLSQP.out')}
 if a.optimizer=='SLSQP' else {'XMAX':.05,'TOLX':1e-6,'TOLC':1e-6,'TOLG':1e-4,
 'MIT':a.maxit,'MFV':150,'MET':1,'MEC':2,'IPRINT':2,'IFILE':str(root/'PSQP.out')})
settings={'description':'Minimize upward-positive CL at baseline CD cap; slip ground',
 'reference_chord_m':c0,'reference_area_m2':c0*span,'span_m':span,
 'speed_ms':V,'temperature_K':T,'pressure_Pa':P,'initial_incidence_deg':-3.,
 'incidence_explicit_bounds':None,'shape_explicit_bounds':None,
 'box_m':{'xmin':-.125,'xmax':.625,'ymin':.05,'ymax':.35},
 'max_section_thickness_min_over_initial_chord':.06,
 'nearfield_normal_nodes':81,'nearfield_march_distance_m':.025,'first_wall_spacing_m':2e-6,
 'shape_variables':20,'incidence_variables':1,'ffd_control_points':40,
 'ffd_span_planes_tied':True,'endpoint_equalities':False,
 'thickness_definition':'Maximum chord-normal paired upper/lower material-station separation',
 'geometry_constraints_on':'All 289 actual CFD surface vertices; edges linearly connect them',
 'discretization_cells_before_overset_blanking':88548,'mpi_ranks':comm.size,
 'cfd_options':options,'mesh_options':meshoptions,'optimizer':'pyOptSparse '+a.optimizer,
 'optimizer_options':optimizeroptions,'derivatives':'ADflow discrete adjoints + pyGeo analytic geometry',
 'sensitivity_overset_mode':'frozen, consistent with primal USMesh warping',
 'mesh_independent':False,'ground_boundary':'symmetry/slip at y=0; no viscous road boundary layer'}
dump(root/'settings.json',settings)
if comm.rank==0:
    (root/'settings.txt').write_text('CFD AND OPTIMIZER SETTINGS\n\n'+json.dumps(settings,indent=2)+'\n')
solver=ADFLOW(options=options,comm=comm)
solver.setMesh(USMesh(options=meshoptions,comm=comm))
solver.setDVGeo(geo)
solver.addFamilyGroup('airfoil_only',['airfoil'])
for f in ['cl','cd','cmz']:
    solver.addFunction(f,'airfoil_only',name=f+'_wing')
solver.addSlices('z',[span/2],sliceType='absolute',groupName='airfoil_only')
ap=AeroProblem(name='wing',V=V,P=P,T=T,alpha=0.,areaRef=c0*span,chordRef=c0,
               xRef=pivot[0],yRef=pivot[1],zRef=span/2,
               evalFuncs=['cl_wing','cd_wing','cmz_wing'])
# Baseline data only uses the wing force families. The slip wall carries no
# viscous drag, but pressure forces on the infinite ground are not wing forces.
count=0
best=None
cdcap=None
cache={}
anchor_state=None
previous_failed=False

def save_state(folder, record, pts, x):
    folder=root/folder
    solver.writeSolution(outputDir=str(folder),baseName='flow',number=0)
    if comm.rank==0:
        np.savetxt(folder/'surface_geometry.csv',pts,delimiter=',',
                   header='x_m,y_m,z_m',comments='')
        (folder/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        (folder/'design.json').write_text(json.dumps(realdict(x),indent=2)+'\n')

def objfun(x, force=False):
    global count,best,cdcap,cache,anchor_state,previous_failed
    if previous_failed and anchor_state is not None:
        solver.setStates(anchor_state)
        previous_failed=False
    geo.setDesignVars(x)
    pts=geo.update('surface_check')
    vals=geometry_values(pts)
    # Do not invoke CFD for a folded/crossed section or grossly displaced mesh.
    # Invalid topology is a numerical failure, not an extra design constraint.
    xy=pts[:,:2]
    sep=xy[:137]-xy[272:135:-1]
    ch=(xy[0]+xy[272])/2-xy[136]
    signed=(sep[:,1]*ch[0]-sep[:,0]*ch[1])/np.linalg.norm(ch)
    invalid=bool(np.min(signed[1:-1])<0 or np.min(xy[:,1])<.015)
    tic=time.perf_counter()
    dump(root/'incoming_design.json',realdict(x))
    if invalid:
        vals.update({'objective':10.,'drag':10.})
        if comm.rank==0:print('INVALID_GEOMETRY',realdict(x),flush=True)
        return vals,True
    solver(ap,writeSolution=False)
    f={};solver.evalFunctions(ap,f)
    cl=float(f['wing_cl_wing']);cd=float(f['wing_cd_wing']);cmz=float(f['wing_cmz_wing'])
    ri,rs,rf=map(float,solver.getResNorms())
    failed=bool(ap.solveFailed or ap.fatalFail or not np.isfinite([cl,cd,cmz]).all())
    if not failed:anchor_state=solver.getStates().copy()
    previous_failed=failed
    if cdcap is None:
        if failed:raise RuntimeError('Initial CFD baseline did not converge')
        replay=root/'source/replay_settings.json'
        cdcap=json.loads(replay.read_text())['drag_cap'] if replay.exists() else cd
        settings['drag_cap']=cdcap
        settings['flow_derived']={'density_kg_m3':float(ap.rho),'viscosity_Pa_s':float(ap.mu),
                                 'mach':float(ap.mach),'reynolds':float(ap.rho*V*c0/ap.mu)}
        dump(root/'settings.json',settings)
        if comm.rank==0:
            (root/'settings.txt').write_text('CFD AND OPTIMIZER SETTINGS\n\n'+json.dumps(settings,indent=2)+'\n')
    vals.update({'objective':cl,'drag':cd/cdcap})
    feasible=bool(not failed and np.min(vals['surface_x'])>=-.25-1e-6
        and np.max(vals['surface_x'])<=1.25+1e-6
        and np.min(vals['surface_y'])>=.1-1e-6 and np.max(vals['surface_y'])<=.7+1e-6
        and vals['max_thickness'][0]>=.06-1e-6 and cd/cdcap<=1+1e-5)
    q=.5*ap.rho*V**2
    record={'evaluation':count,'cl':cl,'cd':cd,'cmz':cmz,
            'downforce_N_per_m':float(-cl*q*c0),'drag_N_per_m':float(cd*q*c0),
            'downforce_to_drag':float(-cl/cd),'force_down_N_extruded_span':float(-cl*q*c0*span),
            'residual_initial':ri,'residual_start':rs,'residual_final':rf,
            'residual_ratio_initial':rf/ri,'cfd_failed':failed,'feasible':feasible,
            'evaluation_wall_time_s':time.perf_counter()-tic,'elapsed_wall_time_s':time.perf_counter()-start,
            'incidence_variable_deg':float(x['incidence'][0]),**metrics(pts),'design':realdict(x)}
    if comm.rank==0:
        with (root/'history.jsonl').open('a') as f:
            f.write(json.dumps(record)+'\n')
        dump(root/'latest.json',record)
        print('OPT_EVAL',json.dumps(record),flush=True)
    if count==0:save_state('baseline',record,pts,x)
    if feasible and (best is None or cl<best['cl']):
        best=copy.deepcopy(record)
        save_state('best_feasible',record,pts,x)
    if failed:
        # pyOptSparse SLSQP ignores fail flags. Explicitly reject failed CFD
        # with finite penalties instead of optimizing unconverged loads.
        vals.update({'objective':1000.,'drag':1000.})
    cache={'record':record,'values':vals,'points':pts.copy(),'design':copy.deepcopy(x)}
    count+=1
    return vals,failed

def sensfun(x, funcs):
    geo.setDesignVars(x)
    pts=geo.update('surface_check')
    g=geometry_sens(pts)
    # Frozen connectivity is the documented USMesh approach when the entire
    # overset grid is warped together. The adjoint matches this primal model.
    if cache.get('record',{}).get('cfd_failed',False):
        zeros={k:np.zeros((1,len(v))) for k,v in x.items()}
        g.update({'objective':copy.deepcopy(zeros),'drag':copy.deepcopy(zeros)})
        return g,True
    fs={};solver.evalFunctionsSens(ap,fs,evalFuncs=['cl_wing','cd_wing'])
    failed=bool(ap.adjointFailed)
    if failed:raise RuntimeError('Adjoint failed: refusing an unconverged optimization gradient')
    g['objective']={k:np.asarray(v) for k,v in fs['wing_cl_wing'].items()}
    g['drag']={k:np.asarray(v)/cdcap for k,v in fs['wing_cd_wing'].items()}
    if comm.rank==0:
        print('OPT_GRADIENT',count,'failed',failed,flush=True)
        with (root/'gradients.jsonl').open('a') as f:
            f.write(json.dumps({'evaluation_count':count,'failed':failed,
                'objective':realdict(g['objective']),'drag':realdict(g['drag'])})+'\n')
    return g,failed

# Check one independent CFD derivative before starting the optimization.
# Small +/- incidence perturbations use the same consistent frozen connectivity.
# Restore the baseline primal state before SLSQP resumes.
def aerodynamic_derivative_check(x0):
    f0,fail=objfun(x0)
    g,adjfail=sensfun(x0,f0)
    if adjfail:raise RuntimeError('Baseline adjoint did not converge')
    h=1e-2
    xp=copy.deepcopy(x0);xm=copy.deepcopy(x0)
    xp['incidence']+=h;xm['incidence']-=h
    fp,pfail=objfun(xp);fm,mfail=objfun(xm)
    report={}
    for name in ['objective','drag']:
        fd=float((fp[name]-fm[name])/(2*h))
        adj=float(np.asarray(g[name]['incidence']).ravel()[0])
        report[name]={'central_difference_per_deg':fd,'adjoint_per_deg':adj,
                      'relative_error':abs(fd-adj)/max(abs(fd),abs(adj),1e-8)}
    report['primal_failed']=pfail or mfail
    dump(root/'source/aerodynamic_derivative_check.json',report)
    if comm.rank==0:print('AERODYNAMIC_DERIVATIVE_CHECK',report,flush=True)
    if pfail or mfail or max(report[n]['relative_error'] for n in ['objective','drag'])>.05:
        raise RuntimeError('Aerodynamic derivative check failed; see saved report')
    objfun(x0)

# Baseline always starts from the same NACA shape, even for a continuation.
basefuncs,fail=objfun(xinitial)
if fail: raise RuntimeError('Baseline CFD failed')
if a.reuse_checks:
    report=json.loads((root/'source/aerodynamic_derivative_check.json').read_text())
    assert not report['primal_failed'] and max(report[n]['relative_error'] for n in ['objective','drag'])<.05
else:
    aerodynamic_derivative_check(xinitial)
if a.evaluate_design:
    target={k:np.array(v) for k,v in json.loads(a.evaluate_design.read_text()).items()}
    funcs,failed=objfun(target)
    save_state('final',cache['record'],cache['points'],target)
    dump(root/'run_complete.json',{'inform':{'value':-999,'text':'Independent design feasibility evaluation; not optimizer convergence'},
        'final_failed':failed,'final':cache['record'],'best_feasible':best,'total_wall_time_s':time.perf_counter()-start})
    raise SystemExit
xstart=copy.deepcopy(xinitial)
if a.restart_design:
    xstart={k:np.array(v) for k,v in json.loads(a.restart_design.read_text()).items()}
    geo.setDesignVars(xstart)
optprob=Optimization('Slip-ground downforce with surface box',objfun,comm=comm)
geo.addVariablesPyOpt(optprob)
optprob.addObj('objective')
optprob.addCon('drag',upper=1.)
optprob.addConGroup('surface_x',len(points),lower=-.25,upper=1.25)
optprob.addConGroup('surface_y',len(points),lower=.1,upper=.7)
optprob.addCon('max_thickness',lower=.06)
if comm.rank==0:
    (root/'problem.txt').write_text(str(optprob))
optimizer=OPT(a.optimizer,options=optimizeroptions)
solution=optimizer(optprob,sens=sensfun,storeHistory=str(root/'optimization.hst'),
                   hotStart=str(a.hotstart) if a.hotstart else None)
if comm.rank==0:
    (root/'solution.txt').write_text(str(solution))
    print(solution,flush=True)
dump(root/'optimizer_exit.json',{'inform':solution.optInform,'variables':realdict(solution.xStar),
                                'optimizer_wall_time_s':float(solution.optTime),'total_wall_time_s':time.perf_counter()-start,'evaluations':count})
# Re-solve and write the optimizer's returned design, keeping best feasible data
# separately in case SLSQP terminates at an infeasible line-search point.
finalfuncs,finalfail=objfun(solution.xStar)
save_state('final',cache['record'],cache['points'],solution.xStar)
dump(root/'run_complete.json',{'inform':solution.optInform,'final_failed':finalfail,
       'final':cache['record'],'best_feasible':best,'optimizer_wall_time_s':float(solution.optTime),'total_wall_time_s':time.perf_counter()-start})

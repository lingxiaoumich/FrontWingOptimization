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
from airfoil_geometry import maximum_chord_normal_thickness, boundary_self_intersects

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, required=True)
p.add_argument('--geometry-only', action='store_true')
p.add_argument('--maxit', type=int, default=50)
p.add_argument('--restart-design', type=Path)
p.add_argument('--hotstart', type=Path)
p.add_argument('--reuse-checks', action='store_true')
p.add_argument('--evaluate-design', type=Path)
p.add_argument('--optimizer',choices=['SLSQP','PSQP'],default='PSQP')
p.add_argument('--ffd-stations',type=int,default=5)
p.add_argument('--shape-directions',choices=['y','xy'],default='xy')
p.add_argument('--drag-cap',type=float,default=.05)
p.add_argument('--save-iterate-flows',action='store_true',
               help='write the flow, design and surface of every PSQP major iterate to iterates/eval_NNN')
p.add_argument('--iterate-volume-only',action='store_true',
               help='with --save-iterate-flows, write only the volume solution and slices (14 of 35 MB)')
p.add_argument('--min-local-thickness',type=float,default=0.,
               help='m; lower bound on local thickness at every paired surface vertex '
                    'away from the leading and trailing edges (0 disables)')
p.add_argument('--thickness-measure',choices=['pair','gap'],default='pair',
               help='local thickness as the chord-normal distance between vertex k and its mirror 272-k (pair), '
                    'or as the signed distance from each of those vertices to the other surface (gap)')
p.add_argument('--max-evaluations',type=int,default=150,help='PSQP MFV: cap on function evaluations')
p.add_argument('--xmax',type=float,default=.05,help='PSQP maximum step length in scaled variables')
p.add_argument('--ncycles',type=int,default=4000,help='ADflow nCycles: caps total (incl. linear) iterations')
p.add_argument('--restart-flow',type=Path,help='converged volume solution at --restart-design; '
               'the first evaluation is then the restart design instead of the NACA baseline')
p.add_argument('--solver',choices=['main','ank_uncoupled','ank_deep'],default='main',
               help='flow-solver variant from the solver study (study/solver-mesh-adjoint)')
p.add_argument('--nk-switch-tol',type=float,help='override NKSwitchTol of the solver variant')
p.add_argument('--l2-convergence',type=float,help='override L2Convergence (default 1e-10)')
p.add_argument('--adjoint-l2-convergence',type=float,help='override adjointL2Convergence (default 1e-7)')
p.add_argument('--baseline',type=Path,
               help='JSON from remesh_design.py: start from a remeshed design, with its FFD lattice '
                    '(ffd.xyz beside the JSON) and the original incidence pivot; design variables restart at zero')
p.add_argument('--warp-err-tol',type=float,
               help="IDWarp errTol (default 5e-4); its far-field approximation can invert large far-field cells")
p.add_argument('--nk-fallback-tol',type=float,
               help='if a flow solve stalls, re-solve the design once from the last converged state with this NKSwitchTol')
p.add_argument('--split-stalled-solves',type=int,default=0,
               help='if a flow solve stalls, reach the design from the last converged one through '
                    'intermediate designs, halving a stalled step down to 2^-N of the way (replaces '
                    '--nk-fallback-tol); the step length that worked is reused for later evaluations')
a = p.parse_args()
if a.restart_flow and not (a.restart_design and a.reuse_checks):
    p.error('--restart-flow needs --restart-design and --reuse-checks (the checks solve the NACA baseline)')
if a.reuse_checks and not a.restart_flow and not a.baseline:
    p.error('--reuse-checks without --restart-flow is only for a remeshed --baseline (checked off its zero point)')
if a.ffd_stations < 2:
    p.error('--ffd-stations must be at least 2')
if not np.isfinite(a.drag_cap) or a.drag_cap <= 0:
    p.error('--drag-cap must be finite and positive')
root = a.root.resolve()
if (root/'run_complete.json').exists():
    raise RuntimeError('This run is complete. Use a fresh --root to preserve its results.')
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
baseline = json.loads(a.baseline.read_text()) if a.baseline else None
if baseline:pivot = np.array(baseline['pivot_m']+[0.])

# Five upper/lower pairs by default, extruded across two tied span planes.
# Each projected point has independent global x/y motion. The enclosure is
# distinct from the regulation box, which constrains actual surface geometry.
nstation = a.ffd_stations
ffd = np.empty((nstation,2,2,3))
for i,x in enumerate(np.linspace(-.01,c0+.01,nstation)):
    for j,y in enumerate([-.06,.06]):
        for k,z in enumerate([-.005,span+.005]):
            ffd[i,j,k] = [np.cos(angle0)*(x-.25*c0)-np.sin(angle0)*y+.25*c0,
                          np.sin(angle0)*(x-.25*c0)+np.cos(angle0)*y+offset,z]
if comm.rank == 0:
    if baseline:
        # The remeshed design's lattice is the deformed one; its x/y moves restart at zero.
        (root/'source/ffd.xyz').write_text((a.baseline.parent/'ffd.xyz').read_text())
    else:
        with (root/'source/ffd.xyz').open('w') as f:
            f.write(f'1\n{nstation} 2 2\n')
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
shape_groups = []
for axis in ([0,1] if a.shape_directions=='xy' else [1]):
    name = 'shape_'+('x' if axis==0 else 'y') if a.shape_directions=='xy' else 'shape'
    direction = np.eye(3)[axis]
    shapes = [{int(idx[i,j,k]):direction.copy() for k in range(2)}
              for i in range(nstation) for j in range(2)]
    geo.addShapeFunctionDV(name, shapes, scale=100.)
    shape_groups.append(name)
points = np.c_[xy0,np.full(len(xy0),span/2)]
geo.addPointSet(points, 'surface_check')
xinitial = copy.deepcopy(geo.getValues())

def load_design(path):
    design={k:np.asarray(v,dtype=float) for k,v in json.loads(path.read_text()).items()}
    if set(design)!=set(xinitial) or any(design[k].shape!=xinitial[k].shape for k in xinitial):
        raise ValueError('Design variable groups/dimensions do not match this FFD setup; do not reuse the old 10-station y-only design/history.')
    if not all(np.isfinite(v).all() for v in design.values()):
        raise ValueError('Design variables must be finite')
    return design

# Geometry functions are evaluated directly on the actual CFD boundary, not FFD.
# Thickness uses actual intersections at a common chordwise coordinate, so
# independent horizontal movement of upper/lower material points is supported.
# Local thickness pairs vertex k with its mirror 272-k (same initial chordwise
# station) and measures their separation normal to the current chord; it turns
# negative where the surfaces cross. Pairs 15..124 skip the leading edge and the
# last 2-3 % of chord, where the section is thin by design.
# At strong camber and incidence the two surfaces slide along the chord relative to
# each other, and a pair can read 2 mm while vertex k touches a different part of the
# other surface. The gap measure instead takes each of those 220 vertices' signed
# distance to the nearest point of the other surface's polyline: positive on the
# material side, negative once the surfaces cross.
local_pairs=np.arange(15,125)
def surface_gap(xy,vertices,other):
    q=xy[vertices][:,None,:]; A=other[None,:-1]; AB=(other[1:]-other[:-1])[None]
    t=((q-A)*AB).sum(-1)/(AB*AB).sum(-1)
    t=np.where(t.real<0,0,np.where(t.real>1,1,t))
    r=A+t[...,None]*AB-q; d=np.sqrt((r*r).sum(-1))
    side=np.sign((AB[...,0]*(q-A)[...,1]-AB[...,1]*(q-A)[...,0]).real)
    j=np.argmin(d.real,axis=1); i=np.arange(len(vertices))
    # Both surfaces run in vertex order around the closed section, so the material is
    # on the same side of every segment: the left for a counterclockwise boundary.
    x,y=xy[:,0].real,xy[:,1].real
    orientation=np.sign(x@np.roll(y,-1)-y@np.roll(x,-1))
    return orientation*side[i,j]*d[i,j]
def gap_thickness(pts):
    xy=pts[:,:2]
    return np.r_[surface_gap(xy,local_pairs,xy[136:273]),surface_gap(xy,272-local_pairs,xy[:137])]
def local_thickness(pts):
    if a.thickness_measure=='gap':return gap_thickness(pts)
    xy=pts[:,:2]; ch=(xy[0]+xy[272])/2-xy[136]; ch=ch/np.sqrt(ch@ch)
    return -(xy[272-local_pairs]-xy[local_pairs])@np.array([-ch[1],ch[0]])
n_local=len(local_pairs)*(2 if a.thickness_measure=='gap' else 1)

def geometry_values(pts):
    xy=pts[:,:2]
    vals={'surface_x':xy[:,0]/c0, 'surface_y':xy[:,1]/c0,
          'max_thickness':np.array([maximum_chord_normal_thickness(pts)/c0])}
    if a.min_local_thickness>0:vals['local_thickness']=local_thickness(pts)/c0
    return vals

def metrics(pts):
    vals=geometry_values(pts)
    xy=pts[:,:2]; chord=(xy[0]+xy[272])/2-xy[136]
    return {'minimum_clearance_m':float(np.min(xy[:,1])),
            'maximum_height_m':float(np.max(xy[:,1])),
            'minimum_x_m':float(np.min(xy[:,0])), 'maximum_x_m':float(np.max(xy[:,0])),
            'maximum_section_thickness_m':float(vals['max_thickness'][0]*c0),
            'maximum_section_thickness_over_reference_chord':float(vals['max_thickness'][0]),
            'minimum_local_thickness_m':float(local_thickness(pts).min()),
            'minimum_surface_gap_m':float(gap_thickness(pts).min()),
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
    _,dp=maximum_chord_normal_thickness(pts,return_gradient=True)
    out['max_thickness']=geo.totalSensitivity(dp/c0,'surface_check')
    if a.min_local_thickness>0:
        dl=np.zeros((n_local,len(pts),3))
        for i in range(len(pts)):
            for axis in range(2):
                pc=pts.astype(complex);pc[i,axis]+=1e-30j
                dl[:,i,axis]=local_thickness(pc).imag/1e-30/c0
        out['local_thickness']=geo.totalSensitivity(dl,'surface_check')
    return out

# Verify geometry derivatives independently at a small nonzero perturbation.
xcheck=copy.deepcopy(xinitial)
for j,name in enumerate(shape_groups):
    xcheck[name]=1e-4*np.sin(np.arange(2*nstation)*.71+(j+1)*.37)
xcheck['incidence']=np.array([-3.001])
geo.setDesignVars(xcheck);pc=geo.update('surface_check')
gs=geometry_sens(pc)
checks={}
for key,step in [('incidence',1e-4)]+[(name,1e-7) for name in shape_groups]:
    direction=np.cos(np.arange(len(xcheck[key]))*.43+.29);direction/=np.linalg.norm(direction)
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
geometry_setup={'ffd_chordwise_stations':nstation,'ffd_dimensions':[nstation,2,2],
    'ffd_control_points':4*nstation,'ffd_projected_points':2*nstation,
    'ffd_span_planes_tied':True,'shape_directions':a.shape_directions,
    'shape_groups':shape_groups,'shape_variables':2*nstation*len(shape_groups),
    'incidence_variables':1,'total_design_variables':2*nstation*len(shape_groups)+1,
    'shape_units':'m','incidence_units':'deg','drag_cap':a.drag_cap,
    'reference_chord_m':c0,'max_section_thickness_min_m':.06*c0,
    'thickness_definition':'Maximum continuous interior chord-normal segment of the actual piecewise-linear surface',
    'box_m':{'xmin':-.125,'xmax':.625,'ymin':.05,'ymax':.35}}
dump(root/'source/geometry_setup.json',geometry_setup)
dump(root/'source/initial_design.json',realdict(xinitial))
if a.geometry_only:
    initial_points=geo.update('surface_check')
    dump(root/'geometry_validation.json',{'status':'Geometry and derivatives checked; no CFD or optimization run',
        'setup':geometry_setup,'maximum_derivative_errors':checks,'baseline_metrics':metrics(initial_points)})
    if comm.rank==0:
        np.savetxt(root/'baseline/surface_geometry.csv',initial_points,delimiter=',',header='x_m,y_m,z_m',comments='')
        print('GEOMETRY_CHECK',checks,metrics(initial_points),flush=True)
    raise SystemExit

options={
    'gridFile':str(root/'source/ground.cgns'),'outputDirectory':str(root/'evaluations'),
    'equationType':'RANS','equationMode':'steady','liftIndex':2,'turbulenceModel':'SA',
    'turbulenceOrder':'first order','discretization':'central plus scalar dissipation',
    'useWallFunctions':False,'useANKSolver':True,'useNKSolver':True,
    'ANKSwitchTol':10.,'NKSwitchTol':1e-4,'ANKCoupledSwitchTol':1e-2,
    'ANKSecondOrdSwitchTol':1e-2,'ANKADPC':True,'NKADPC':True,'NKViscPC':True,
    'NKSubspaceSize':500,'NKPCILUFill':5,'NKLinearSolveTol':.1,'NKUseEW':False,
    'turbResScale':10000.,'nCycles':a.ncycles,'L2Convergence':1e-10,'L2ConvergenceRel':1e-10,
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
if a.restart_flow:
    options['restartFile']=str(a.restart_flow.resolve())
# ank_uncoupled skips coupled ANK; ank_deep also hands over to NK four decades later
# with a larger Krylov space and stronger ILU, which cured the study's NK stalls.
options.update({'main':{},'ank_uncoupled':{'ANKCoupledSwitchTol':1e-16},
    'ank_deep':{'ANKCoupledSwitchTol':1e-16,'NKSwitchTol':1e-8,'NKSubspaceSize':400,
                'NKPCILUFill':4,'NKASMOverlap':3}}[a.solver])
for flag,key in [('nk_switch_tol','NKSwitchTol'),('l2_convergence','L2Convergence'),
                 ('adjoint_l2_convergence','adjointL2Convergence')]:
    if getattr(a,flag) is not None:options[key]=getattr(a,flag)
meshoptions={'gridFile':options['gridFile'],
    'symmetryPlanes':[[[0.,0.,0.],[0.,1.,0.]],[[0.,0.,0.],[0.,0.,1.]],
                      [[0.,0.,span],[0.,0.,1.]]]}
if a.warp_err_tol is not None:meshoptions['errTol']=a.warp_err_tol
optimizeroptions=({'ACC':1e-5,'MAXIT':a.maxit,'IPRINT':2,'IFILE':str(root/'SLSQP.out')}
 if a.optimizer=='SLSQP' else {'XMAX':a.xmax,'TOLX':1e-6,'TOLC':1e-6,'TOLG':1e-4,
 'MIT':a.maxit,'MFV':a.max_evaluations,'MET':1,'MEC':2,'IPRINT':2,'IFILE':str(root/'PSQP.out')})
settings={'description':'Minimize upward-positive CL with prescribed CD cap; slip ground',
 'reference_chord_m':c0,'reference_area_m2':c0*span,'span_m':span,
 'speed_ms':V,'temperature_K':T,'pressure_Pa':P,'initial_incidence_deg':-3.,
 'incidence_explicit_bounds':None,'shape_explicit_bounds':None,
 'box_m':{'xmin':-.125,'xmax':.625,'ymin':.05,'ymax':.35},
 'max_section_thickness_min_over_initial_chord':.06,
 'min_local_thickness_m':a.min_local_thickness,'local_thickness_vertex_pairs':'k and 272-k, k=15..124',
 'local_thickness_measure':{'pair':'chord-normal distance between vertex k and 272-k',
     'gap':'signed distance from vertex k and from 272-k to the nearest point of the other surface'}[a.thickness_measure],
 'nearfield_normal_nodes':81,'nearfield_march_distance_m':.025,'first_wall_spacing_m':2e-6,
 **geometry_setup,
 'ffd_span_planes_tied':True,'endpoint_equalities':False,
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
cdcap=a.drag_cap
cache={}
anchor_state=None
anchor_x=None
previous_failed=False
# Distance (in PSQP's scaled variables, where XMAX applies) that one solve can bridge from
# a converged flow; unknown until a solve stalls.
reach=np.inf
dv_scale={'incidence':1.,**{name:100. for name in shape_groups}}

def scaled_distance(x0,x1):
    return float(np.sqrt(sum((((np.real(x1[k])-np.real(x0[k]))*dv_scale[k])**2).sum() for k in x1)))

def walk(x):
    """Solve x starting from the last converged design, through intermediate designs on the
    straight line between them, each step no longer than the reach. A stalled step is
    halved and retried from its own converged start. The final mesh is the same as a
    direct solve: IDWarp always warps from the original grid."""
    global anchor_state,anchor_x,reach
    x0=copy.deepcopy(anchor_x);d=scaled_distance(x0,x)
    s,h,hmin=0.,1.,2.**-a.split_stalled_solves
    if d>0:h=max(min(1.,reach/d),hmin)
    solves=stalls=0
    while True:
        t=min(1.,s+h)
        geo.setDesignVars({k:np.real(x0[k])+t*(np.real(x[k])-np.real(x0[k])) for k in x})
        solver(ap,writeSolution=False);solves+=1
        if ap.fatalFail:break  # e.g. negative volumes: a shorter step cannot reach x either
        if not ap.solveFailed:
            # Every converged design is a valid start for the next evaluation.
            s=t;anchor_state=solver.getStates().copy()
            anchor_x={k:np.real(x0[k])+s*(np.real(x[k])-np.real(x0[k])) for k in x}
            if s>=1:break
            continue
        stalls+=1
        if comm.rank==0:print('STALLED_STEP',json.dumps({'from':s,'to':t,'distance':d*(t-s)}),flush=True)
        if h/2<hmin:break
        h/=2;solver.setStates(anchor_state)
    geo.setDesignVars(x)
    if stalls:reach=d*h if s>=1 else d*h/2
    elif s>=1:reach*=1.25  # grow back slowly after a clean evaluation
    if reach>4*a.xmax:reach=np.inf
    return {'solves':solves,'stalled_solves':stalls,'reached_fraction':s,'distance':d,
            'reach_after':reach if np.isfinite(reach) else None}

def save_state(folder, record, pts, x, surfaces=True):
    folder=root/folder
    if not surfaces:
        kept={k:solver.getOption(k) for k in ['writeSurfaceSolution','writeTecplotSurfaceSolution']}
        for k in kept:solver.setOption(k,False)
    solver.writeSolution(outputDir=str(folder),baseName='flow',number=0)
    if not surfaces:
        for k,v in kept.items():solver.setOption(k,v)
    if comm.rank==0:
        np.savetxt(folder/'surface_geometry.csv',pts,delimiter=',',
                   header='x_m,y_m,z_m',comments='')
        (folder/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        (folder/'design.json').write_text(json.dumps(realdict(x),indent=2)+'\n')

def objfun(x, force=False):
    global count,best,cdcap,cache,anchor_state,anchor_x,previous_failed
    if previous_failed and anchor_state is not None:
        solver.setStates(anchor_state)
        previous_failed=False
    geo.setDesignVars(x)
    pts=geo.update('surface_check')
    # Do not invoke CFD for a folded/crossed section or grossly displaced mesh.
    # Invalid topology is a numerical failure, not an extra design constraint.
    xy=pts[:,:2]
    invalid=bool(boundary_self_intersects(pts) or np.min(xy[:,1])<.015)
    try:
        vals=geometry_values(pts) if not invalid else None
    except ValueError:
        invalid=True
    tic=time.perf_counter()
    dump(root/'incoming_design.json',realdict(x))
    if invalid:
        vals={'surface_x':xy[:,0]/c0,'surface_y':xy[:,1]/c0,
              'max_thickness':np.array([-1.]),'objective':10.,'drag':10.}
        if a.min_local_thickness>0:vals['local_thickness']=np.full(n_local,-1.)
        if comm.rank==0:print('INVALID_GEOMETRY',realdict(x),flush=True)
        # PSQP may still ask for a gradient here; without this, sensfun would run an
        # adjoint at the previous valid state and pair it with this geometry.
        cache={'record':{'cfd_failed':True}}
        return vals,True
    path={}
    if a.split_stalled_solves and anchor_state is not None:path=walk(x)
    else:solver(ap,writeSolution=False)
    fallback=False
    if not a.split_stalled_solves and a.nk_fallback_tol and ap.solveFailed and not ap.fatalFail and anchor_state is not None:
        # An early NK handover is 3x faster but can stall on a large step; retry the
        # same design once with the later, robust handover before reporting failure.
        if comm.rank==0:print('NK_FALLBACK',realdict(x),flush=True)
        fallback=True
        solver.setStates(anchor_state)
        solver.setOption('NKSwitchTol',a.nk_fallback_tol)
        solver(ap,writeSolution=False)
        solver.setOption('NKSwitchTol',options['NKSwitchTol'])
    f={};solver.evalFunctions(ap,f)
    cl=float(f['wing_cl_wing']);cd=float(f['wing_cd_wing']);cmz=float(f['wing_cmz_wing'])
    ri,rs,rf=map(float,solver.getResNorms())
    failed=bool(ap.solveFailed or ap.fatalFail or not np.isfinite([cl,cd,cmz]).all())
    if not failed:anchor_state,anchor_x=solver.getStates().copy(),copy.deepcopy(x)
    previous_failed=failed
    if count==0:
        if failed:raise RuntimeError('Initial CFD baseline did not converge')
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
        and vals['max_thickness'][0]>=.06-1e-6 and cd/cdcap<=1+1e-5
        and (a.min_local_thickness<=0 or vals['local_thickness'].min()>=a.min_local_thickness/c0-1e-6))
    q=.5*ap.rho*V**2
    record={'evaluation':count,'cl':cl,'cd':cd,'cmz':cmz,
            'downforce_N_per_m':float(-cl*q*c0),'drag_N_per_m':float(cd*q*c0),
            'downforce_to_drag':float(-cl/cd),'force_down_N_extruded_span':float(-cl*q*c0*span),
            'residual_initial':ri,'residual_start':rs,'residual_final':rf,
            'residual_ratio_initial':rf/ri,'cfd_failed':failed,'nk_fallback':fallback,**path,'feasible':feasible,
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
    if a.save_iterate_flows:
        # A gradient request marks a major iterate; the solver still holds its state.
        folder=f"iterates/eval_{cache['record']['evaluation']:03d}"
        if comm.rank==0:(root/folder).mkdir(parents=True,exist_ok=True)
        comm.Barrier()
        save_state(folder,cache['record'],cache['points'],x,surfaces=not a.iterate_volume_only)
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

# Baseline always starts from the same NACA shape, even for a continuation,
# unless --restart-flow supplies a converged flow at the restart design: solving
# the NACA baseline first would force a large jump back to the restart design.
basefuncs,fail=objfun(load_design(a.restart_design) if a.restart_flow else xinitial)
if fail: raise RuntimeError('Baseline CFD failed')
if a.reuse_checks:
    report=json.loads((root/'source/aerodynamic_derivative_check.json').read_text())
    assert not report['primal_failed'] and max(report[n]['relative_error'] for n in ['objective','drag'])<.05
else:
    aerodynamic_derivative_check(xinitial)
if a.evaluate_design:
    target=load_design(a.evaluate_design)
    funcs,failed=objfun(target)
    save_state('final',cache['record'],cache['points'],target)
    dump(root/'run_complete.json',{'inform':{'value':-999,'text':'Independent design feasibility evaluation; not optimizer convergence'},
        'final_failed':failed,'final':cache['record'],'best_feasible':best,'total_wall_time_s':time.perf_counter()-start})
    raise SystemExit
xstart=copy.deepcopy(xinitial)
if a.restart_design:
    xstart=load_design(a.restart_design)
    geo.setDesignVars(xstart)
optprob=Optimization('Slip-ground downforce with surface box',objfun,comm=comm)
geo.addVariablesPyOpt(optprob)
optprob.addObj('objective')
optprob.addCon('drag',upper=1.)
optprob.addConGroup('surface_x',len(points),lower=-.25,upper=1.25)
optprob.addConGroup('surface_y',len(points),lower=.1,upper=.7)
optprob.addCon('max_thickness',lower=.06)
if a.min_local_thickness>0:
    optprob.addConGroup('local_thickness',n_local,lower=a.min_local_thickness/c0)
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

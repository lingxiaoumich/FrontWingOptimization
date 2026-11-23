"""Adjoint shape/incidence optimization of the project's slip-ground airfoil.

Execute under mpirun in mdolab/public:u24-gcc-ompi-stable. All dimensional
references are fixed to the initial chord. The bounds apply to surface mesh
vertices (and thus every segment of the piecewise-linear CFD surface).
"""
import argparse
import copy
from dataclasses import asdict, is_dataclass
import json
import os
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
p.add_argument('--fix-incidence',action='store_true',help='keep the incidence at its start: only the shape variables move')
p.add_argument('--drag-cap',type=float,default=.05)
p.add_argument('--loss-cap',type=float,help='upper bound on the outlet Cpt loss (see --outlet-x); replaces the drag cap, '
               'which then must be off (>= 1)')
p.add_argument('--cp-min',type=float,help='lower bound on the wing surface Cp (e.g. -6): ADflow cpmin, a KS '
               'aggregate (sharpness --cp-min-rho) just below the true minimum; the lowest Cp of the field sits on the wing')
p.add_argument('--cp-min-rho',type=float,default=100.,help='KS sharpness of --cp-min (ADflow cpMinRho)')
p.add_argument('--max-clc',type=float,help='m; upper bound on C_L c (upward positive: -1.6 asks for a downforce C_L c of at least '
               '1.6). The objective becomes the outlet Cpt loss (see --outlet-x) instead of C_L; drag and loss caps must be off')
p.add_argument('--first-solve-tilt',type=float,help='deg; with no initial flow (no source/restart_vol.cgns): the first solve, '
               'from freestream, is made at the start design with its incidence lowered by this much, then the start design '
               'is reached from that converged flow by the stalled-solve walk (needs --split-stalled-solves). For the level '
               'start in free air, whose symmetric zero-lift flow stalled at ~2e-9 of the initial residual when solved from '
               'freestream (lift flipping sign at 1e-5), 2026-11-21')
p.add_argument('--no-ground',action='store_true',help='free air: the grid (remesh_design.py --free-air) has a far field '
               '10 chords below y = 0 instead of the slip ground; no y = 0 symmetry plane in the mesh warp, and the outlet '
               'loss plane spans the full height, y = -5..5 m')
p.add_argument('--outlet-x',type=float,default=7.3,help='m; plane x=const over the full domain height (outlet at 7.5 m) on '
               'which the mass-averaged Cpt loss (p0_inf - p0)/q_inf is measured')
# Geometry constraints on the actual surface (0 turns each off)
p.add_argument('--box',type=float,nargs=4,default=[-.125,.625,.05,.35],metavar=('XMIN','XMAX','YMIN','YMAX'),
               help='m; regulation box that holds every surface vertex')
p.add_argument('--min-local-thickness',type=float,default=0.,
               help='m; floor on the gap from each vertex k=15..124 and 272-k to the other surface')
p.add_argument('--te-min-angle',type=float,default=0.,
               help='deg; trailing-edge taper of that floor over vertices 1..14 and their mirrors')
p.add_argument('--te-min-thickness',type=float,default=0.,help='m; floor on the trailing-edge base thickness')
p.add_argument('--near-wall-dist',type=float,default=.01,help='m; ADflow nearWallDist, the overset near-wall distance; '
               '0.4 x the near-grid depth (0.01 with the original 25 mm)')
p.add_argument('--min-area',type=float,default=0.,
               help='m^2; floor on the section area (inside all 289 surface vertices). Absolute, so the runs on '
               'later grids of a chain keep the floor of the original start')
# Optimizer and flow solver
p.add_argument('--xmax',type=float,default=.05,help='PSQP XMAX: maximum step in scaled variables')
p.add_argument('--max-evaluations',type=int,default=150,help='PSQP MFV: maximum function evaluations')
p.add_argument('--ncycles',type=int,default=4000,help='ADflow nCycles')
p.add_argument('--first-ncycles',type=int,help='nCycles of the first solve, e.g. from freestream on a new grid')
p.add_argument('--nk-switch-tol',type=float,help='ADflow NKSwitchTol (default 1e-8)')
p.add_argument('--warp-err-tol',type=float,help='IDWarp errTol (default 5e-4)')
p.add_argument('--split-stalled-solves',type=int,default=0,
               help='reach a design whose solve stalls in steps from the last converged one, halving down to 2^-N')
# Start and output
p.add_argument('--baseline',type=Path,help='baseline.json of remesh_design.py: start from that grid and FFD lattice')
p.add_argument('--check-design',type=Path,help='design for the derivative check (default: the baseline)')
p.add_argument('--restart-flow',type=Path,help='converged flow of --restart-design; the first solve is that design')
p.add_argument('--save-iterate-flows',action='store_true',help='save every PSQP major iterate in iterates/eval_NNN')
p.add_argument('--iterate-volume-only',action='store_true',help='with --save-iterate-flows: skip the surface solution files')
a = p.parse_args()
if a.ffd_stations < 2:
    p.error('--ffd-stations must be at least 2')
if not np.isfinite(a.drag_cap) or a.drag_cap <= 0:
    p.error('--drag-cap must be finite and positive')
if a.loss_cap is not None and not (np.isfinite(a.loss_cap) and a.loss_cap > 0 and a.drag_cap >= 1):
    p.error('--loss-cap must be finite and positive, with the drag cap off (--drag-cap 10)')
if a.max_clc is not None and not (np.isfinite(a.max_clc) and a.max_clc < 0 and a.loss_cap is None and a.drag_cap >= 1):
    p.error('--max-clc must be a finite negative C_L c, without --loss-cap and with the drag cap off (--drag-cap 10)')
# The objective: C_L, or with --max-clc the outlet Cpt loss over LOSS_SCALE (order 1, as C_L is)
LOSS_OBJ = a.max_clc is not None
LOSS_SCALE = .005
# The bounded aerodynamic constraint: drag, the outlet Cpt loss, or C_L c when the loss is the objective
CAP = 'lift' if LOSS_OBJ else 'loss' if a.loss_cap is not None else 'drag'
if a.first_solve_tilt is not None and not (np.isfinite(a.first_solve_tilt) and a.first_solve_tilt and a.split_stalled_solves):
    p.error('--first-solve-tilt must be a finite nonzero tilt, with --split-stalled-solves')
if a.cp_min is not None and not (np.isfinite(a.cp_min) and a.cp_min < 0):
    p.error('--cp-min must be a finite negative Cp')
# Every aerodynamic constraint, each scaled so that <= 1 is feasible
AERO = [CAP]+(['cp_min'] if a.cp_min is not None else [])
box = dict(zip(['xmin','xmax','ymin','ymax'], a.box))
# Absolute but not resolved: a run folder can be a link to another disk whose path has a space,
# which the solver's file names and run_chain.py's commands must not see.
root = Path(os.path.abspath(a.root))
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
# --fix-incidence: both bounds at the start (pyGeo's freezeVars would pass them complex, which PSQP rejects)
fixed = -3. if a.fix_incidence else None
geo.addGlobalDV('incidence', [-3.], incidence, lower=fixed, upper=fixed, scale=1.)
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
# Vertices 0..136 are the upper surface, 136..272 the lower one, 272..288 the blunt
# trailing-edge base. All functions below take complex coordinates for complex step.
def surface_gap(xy,vertices,other):
    """Signed distance from xy[vertices] to the polyline `other`: positive inside the section."""
    q=xy[vertices][:,None,:]; A=other[None,:-1]; AB=(other[1:]-other[:-1])[None]
    t=((q-A)*AB).sum(-1)/(AB*AB).sum(-1)
    t=np.where(t.real<0,0,np.where(t.real>1,1,t))  # clip to the segment
    r=A+t[...,None]*AB-q; d=np.sqrt((r*r).sum(-1))
    side=np.sign((AB[...,0]*(q-A)[...,1]-AB[...,1]*(q-A)[...,0]).real)
    j=np.argmin(d.real,axis=1); i=np.arange(len(vertices))
    # The material is left of every segment for a counterclockwise boundary.
    x,y=xy[:,0].real,xy[:,1].real
    orientation=np.sign(x@np.roll(y,-1)-y@np.roll(x,-1))
    return orientation*side[i,j]*d[i,j]

# Local thickness: the gap from each vertex k=15..124 (upper) and 272-k (lower) to the
# other surface. The leading edge and the last 2-3 % of chord are left out.
local_pairs=np.arange(15,125)
def gap_thickness(pts):
    xy=pts[:,:2]
    return np.r_[surface_gap(xy,local_pairs,xy[136:273]),surface_gap(xy,272-local_pairs,xy[:137])]

# Trailing-edge taper over vertices 1..14 and their mirrors: gap minus
# min(local floor, sin(angle) x distance to the vertex's own base corner, 0 or 272),
# a wedge of that angle up to the base. Without it the optimizer folds the edge flat.
te_pairs=np.arange(1,15)
def te_taper(pts):
    xy=pts[:,:2]; lower=272-te_pairs
    gap=np.r_[surface_gap(xy,te_pairs,xy[136:273]),surface_gap(xy,lower,xy[:137])]
    d=np.r_[xy[te_pairs]-xy[0],xy[lower]-xy[272]]
    floor=np.sin(np.deg2rad(a.te_min_angle))*np.sqrt((d*d).sum(-1))
    if a.min_local_thickness>0:floor=np.where(floor.real>a.min_local_thickness,a.min_local_thickness,floor)
    return gap-floor

# Base thickness: the gap from each base corner (0, 272) to the other surface.
def te_base(pts):
    xy=pts[:,:2]
    return np.r_[surface_gap(xy,np.array([0]),xy[136:273]),surface_gap(xy,np.array([272]),xy[:137])]

# Section area: the shoelace sum over the closed boundary (vertex 288 is vertex 0), positive for
# the counterclockwise boundary; a folded one is rejected before this counts.
def section_area(pts):
    x,y=pts[:,0],pts[:,1]
    return np.array([.5*(x@np.roll(y,-1)-y@np.roll(x,-1))])

# As a constraint: area over the reference chord, the mean thickness in m, so that like the others
# it is a length here and over c0 in the optimizer (area/c0^2).
def mean_thickness(pts):
    return section_area(pts)/c0

# The active vertex constraints: name -> (function in m, lower bound in m)
vertex_constraints={}
if a.min_local_thickness>0:vertex_constraints['local_thickness']=(gap_thickness,a.min_local_thickness)
if a.te_min_angle>0:vertex_constraints['te_taper']=(te_taper,0.)
if a.te_min_thickness>0:vertex_constraints['te_base']=(te_base,a.te_min_thickness)
if a.min_area>0:vertex_constraints['section_area']=(mean_thickness,a.min_area/c0)

def geometry_values(pts):
    xy=pts[:,:2]
    vals={'surface_x':xy[:,0]/c0, 'surface_y':xy[:,1]/c0,
          'max_thickness':np.array([maximum_chord_normal_thickness(pts)/c0])}
    for name,(fun,_) in vertex_constraints.items():vals[name]=fun(pts)/c0
    return vals

def metrics(pts):
    vals=geometry_values(pts)
    xy=pts[:,:2]; chord=(xy[0]+xy[272])/2-xy[136]
    return {'minimum_clearance_m':float(np.min(xy[:,1])),
            'maximum_height_m':float(np.max(xy[:,1])),
            'minimum_x_m':float(np.min(xy[:,0])), 'maximum_x_m':float(np.max(xy[:,0])),
            'maximum_section_thickness_m':float(vals['max_thickness'][0]*c0),
            'maximum_section_thickness_over_reference_chord':float(vals['max_thickness'][0]),
            'minimum_surface_gap_m':float(gap_thickness(pts).min()),
            'te_base_thickness_m':float(te_base(pts).min()),
            'section_area_m2':float(section_area(pts)[0]),
            **({'minimum_te_taper_margin_m':float(te_taper(pts).min())} if a.te_min_angle>0 else {}),
            'actual_chord_m':float(np.linalg.norm(chord)),
            'chord_incidence_deg':float(-np.rad2deg(np.arctan2(chord[1],chord[0]))),
            'box_minimum_margin_m':float(min(xy[:,0].min()-box['xmin'],box['xmax']-xy[:,0].max(),
                                            xy[:,1].min()-box['ymin'],box['ymax']-xy[:,1].max()))}

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
    for name,(fun,_) in vertex_constraints.items():
        # complex step in x and y of every vertex: the nearest point can be anywhere
        dl=np.zeros((len(fun(pts)),len(pts),3))
        for i in range(len(pts)):
            for axis in range(2):
                pc=pts.astype(complex);pc[i,axis]+=1e-30j
                dl[:,i,axis]=fun(pc).imag/1e-30/c0
        out[name]=geo.totalSensitivity(dl,'surface_check')
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
    'box_m':box}
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
    'MGCycle':'sg','useZipperMesh':False,'nearWallDist':a.near_wall_dist,'oversetUpdateMode':'frozen',
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
    options['restartFile']=os.path.abspath(a.restart_flow)
# "ank_deep" of the solver study (study/solver-mesh-adjoint): no coupled ANK, and a later
# hand-over to NK with a larger Krylov space and stronger ILU, which cured NK stalls.
options.update({'ANKCoupledSwitchTol':1e-16,'NKSwitchTol':1e-8,'NKSubspaceSize':400,
                'NKPCILUFill':4,'NKASMOverlap':3})
if a.nk_switch_tol is not None:options['NKSwitchTol']=a.nk_switch_tol
if a.cp_min is not None:options.update({'computeCavitation':True,'cpMinRho':a.cp_min_rho})  # turns on ADflow's cpmin
meshoptions={'gridFile':options['gridFile'],
    'symmetryPlanes':([] if a.no_ground else [[[0.,0.,0.],[0.,1.,0.]]])+[[[0.,0.,0.],[0.,0.,1.]],
                      [[0.,0.,span],[0.,0.,1.]]]}
# --no-ground and a free-air grid go together (remesh_design.py records it in baseline.json 'mesh')
if a.baseline and json.loads(Path(a.baseline).read_text()).get('mesh',{}).get('free_air',False)!=a.no_ground:
    p.error('--no-ground needs a free-air grid (remesh_design.py --free-air), and a free-air grid needs --no-ground')
if a.warp_err_tol is not None:meshoptions['errTol']=a.warp_err_tol
optimizeroptions=({'ACC':1e-5,'MAXIT':a.maxit,'IPRINT':2,'IFILE':str(root/'SLSQP.out')}
 if a.optimizer=='SLSQP' else {'XMAX':a.xmax,'TOLX':1e-6,'TOLC':1e-6,'TOLG':1e-4,
 'MIT':a.maxit,'MFV':a.max_evaluations,'MET':1,'MEC':2,'IPRINT':2,'IFILE':str(root/'PSQP.out')})
settings={'description':(f'Minimize the outlet Cpt loss with C_L c <= {a.max_clc:g} m' if LOSS_OBJ else
                'Minimize upward-positive CL with prescribed '+('outlet Cpt loss cap' if CAP=='loss' else 'CD cap'))+
               ('; free air (no ground)' if a.no_ground else '; slip ground'),
 'aerodynamic_constraint':CAP,
 'objective':f'outlet Cpt loss / {LOSS_SCALE:g}' if LOSS_OBJ else 'cl',
 'max_clc_m':a.max_clc,
 'first_solve_tilt_deg':a.first_solve_tilt,
 'lift_constraint':'2 - C_L c / max_clc_m <= 1 (C_L c <= max_clc_m, both negative)' if LOSS_OBJ else None,
 'outlet_loss':{'plane_x_m':a.outlet_x,'height_m':[-5. if a.no_ground else 0.,5.],'definition':'(p0_inf - mass-averaged p0 on the plane)/q_inf; '
                'flux form in drag units: that x mdot/(rho_inf V span), m','adflow_functions':['mavgptot','mdot']},
 'reference_chord_m':c0,'reference_area_m2':c0*span,'span_m':span,
 'speed_ms':V,'temperature_K':T,'pressure_Pa':P,'initial_incidence_deg':-3.,
 'incidence_explicit_bounds':None,'shape_explicit_bounds':None,'incidence_fixed':a.fix_incidence,
 'box_m':box,
 'max_section_thickness_min_over_initial_chord':.06,
 'min_local_thickness_m':a.min_local_thickness,'local_thickness_vertex_pairs':'k and 272-k, k=15..124',
 'local_thickness_measure':'signed distance from vertex k and from 272-k to the nearest point of the other surface',
 'te_min_angle_deg':a.te_min_angle,'te_min_thickness_m':a.te_min_thickness,
 'te_base':'gap from each base corner (vertex 0, 272) to the other surface >= te_min_thickness',
 'min_section_area_m2':a.min_area,
 'section_area':'area inside the 289 surface vertices (shoelace) >= min_section_area_m2; '
                'in the optimizer as area/c0^2 (constraint section_area)',
 'te_taper':'vertex pairs k and 272-k, k=1..14: gap to the other surface >= min(min_local_thickness, '
            'sin(te_min_angle) x distance from the vertex to its base corner, 0 or 272)',
 'nearfield_normal_nodes':81,'nearfield_march_distance_m':.025,'first_wall_spacing_m':2e-6,
 **geometry_setup,
 'ffd_span_planes_tied':True,'endpoint_equalities':False,
 'geometry_constraints_on':'All 289 actual CFD surface vertices; edges linearly connect them',
 'discretization_cells_before_overset_blanking':88548,'mpi_ranks':comm.size,
 'cfd_options':options,'mesh_options':meshoptions,'optimizer':'pyOptSparse '+a.optimizer,
 'optimizer_options':optimizeroptions,'derivatives':'ADflow discrete adjoints + pyGeo analytic geometry',
 'sensitivity_overset_mode':'frozen, consistent with primal USMesh warping',
 'mesh_independent':False,'ground_boundary':'none: free air, far field 10 chords below y=0 and above' if a.no_ground else
                                 'symmetry/slip at y=0; no viscous road boundary layer'}
dump(root/'settings.json',settings)
if comm.rank==0:
    (root/'settings.txt').write_text('CFD AND OPTIMIZER SETTINGS\n\n'+json.dumps(settings,indent=2)+'\n')
solver=ADFLOW(options=options,comm=comm)
# Outlet plane x = --outlet-x across the whole domain height (5 m) and span, normal +x, finer than
# the background grid there (2.5 mm up to 0.25 m). ADflow interpolates the flow onto it and
# integrates mass flow and mass-averaged total pressure, with adjoint derivatives.
outlet_plane=root/'source/outlet_plane.xyz'
if comm.rank==0:
    yp=np.unique(np.r_[np.linspace(0,.5,401),np.linspace(.5,1.,101),np.linspace(1.,5.,161)])
    if a.no_ground:yp=np.unique(np.r_[-yp,yp])  # the same spacing below y = 0
    with outlet_plane.open('w') as f:
        f.write(f'1\n{len(yp)} 2 1\n')
        for coord in [np.full(2*len(yp),a.outlet_x),np.tile(yp,2),np.repeat([0.,span],len(yp))]:
            np.savetxt(f,coord,fmt='%.16e')
comm.Barrier()
solver.addIntegrationSurface(str(outlet_plane),'outlet_plane',isInflow=False)
solver.finalizeUserIntegrationSurfaces()
solver.setMesh(USMesh(options=meshoptions,comm=comm))
solver.setDVGeo(geo)
solver.addFamilyGroup('airfoil_only',['airfoil'])
for f in ['cl','cd','cmz']:
    solver.addFunction(f,'airfoil_only',name=f+'_wing')
solver.addSlices('z',[span/2],sliceType='absolute',groupName='airfoil_only')
solver.addFunction('mavgptot','outlet_plane',name='ptot_outlet')
if a.cp_min is not None:solver.addFunction('cpmin','airfoil_only',name='cpmin_wing')
solver.addFunction('mdot','outlet_plane',name='mdot_outlet')
ap=AeroProblem(name='wing',V=V,P=P,T=T,alpha=0.,areaRef=c0*span,chordRef=c0,
               xRef=pivot[0],yRef=pivot[1],zRef=span/2,
               evalFuncs=['cl_wing','cd_wing','cmz_wing','ptot_outlet','mdot_outlet']+(['cpmin_wing'] if a.cp_min is not None else []))
# Freestream total pressure exactly as ADflow forms it (constant gamma, the AeroProblem's P and Mach)
gamma=ap.gamma
p0inf=P*(1+(gamma-1)/2*ap.mach**2)**(gamma/(gamma-1))
qinf=.5*ap.rho*V**2
def outlet_loss(f):
    """(mass-averaged Cpt loss on the outlet plane, the same loss as a flux in drag units C_D c, m)"""
    loss=(p0inf-float(f['wing_ptot_outlet']))/qinf
    return loss, loss*float(f['wing_mdot_outlet'])/(ap.rho*V*span)
# Baseline data only uses the wing force families. The slip wall carries no
# viscous drag, but pressure forces on the infinite ground are not wing forces.
count=0
best=None
cdcap=a.drag_cap
cache={}
anchor_state=None
anchor_x=None
previous_failed=False
# Step (in PSQP's scaled variables) one solve can bridge from a converged flow; unknown until a stall.
reach=np.inf
dv_scale={'incidence':1.,**{name:100. for name in shape_groups}}

def scaled_distance(x0,x1):
    return float(np.sqrt(sum((((np.real(x1[k])-np.real(x0[k]))*dv_scale[k])**2).sum() for k in x1)))

def walk(x):
    """Solve x from the last converged design through designs on the straight line between
    them, each step at most `reach`; a stalled step is halved and retried from its start.
    IDWarp always warps from the original grid, so the final mesh is that of a direct solve."""
    global anchor_state,anchor_x,reach
    x0=copy.deepcopy(anchor_x);d=scaled_distance(x0,x)
    s,h,hmin=0.,1.,2.**-a.split_stalled_solves
    if d>0:h=max(min(1.,reach/d),hmin)
    solves=stalls=0
    while True:
        t=min(1.,s+h)
        geo.setDesignVars({k:np.real(x0[k])+t*(np.real(x[k])-np.real(x0[k])) for k in x})
        solver(ap,writeSolution=False);solves+=1
        if ap.fatalFail:break  # e.g. negative volumes: a shorter step cannot help
        if not ap.solveFailed:
            # every converged design is a start for the next evaluation
            s=t;anchor_state=solver.getStates().copy()
            anchor_x={k:np.real(x0[k])+s*(np.real(x[k])-np.real(x0[k])) for k in x}
            if s>=1:break
            continue
        stalls+=1
        if comm.rank==0:print('STALLED_STEP',json.dumps({'from':s,'to':t,'distance':d*(t-s)}),flush=True)
        h=(t-s)/2  # halve the step that stalled (t-s, not h: s+h may have been cut at 1)
        if h<hmin:break
        solver.setStates(anchor_state)
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

def unwarped(x):
    """The grid's own design: incidence -3 deg and no shape move, so IDWarp leaves the grid as made."""
    return bool(float(np.real(x['incidence'][0]))==-3. and not any(np.any(np.real(v)) for k,v in x.items() if k!='incidence'))

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
    invalid=bool(boundary_self_intersects(pts) or np.min(xy[:,1])<.3*box['ymin'])  # 15 mm for the 50 mm floor
    try:
        vals=geometry_values(pts) if not invalid else None
    except ValueError:
        invalid=True
    tic=time.perf_counter()
    dump(root/'incoming_design.json',realdict(x))
    if invalid:
        vals={'surface_x':xy[:,0]/c0,'surface_y':xy[:,1]/c0,
              'max_thickness':np.array([-1.]),'objective':10.,**{n:10. for n in AERO}}
        for name,(fun,_) in vertex_constraints.items():vals[name]=np.full(len(fun(points)),-1.)
        if comm.rank==0:print('INVALID_GEOMETRY',realdict(x),flush=True)
        # a gradient request here must not pair the previous state's adjoint with this geometry
        cache={'record':{'cfd_failed':True}}
        return vals,True
    path={}
    tilt=count==0 and a.first_solve_tilt and anchor_state is None and 'restartFile' not in options
    if tilt:
        # a converged lifting flow first; the walk below then takes it to x (the run's own convergence test throughout)
        w={k:np.real(np.asarray(v,dtype=float)).copy() for k,v in x.items()};w['incidence']=w['incidence']-a.first_solve_tilt
        geo.setDesignVars(w)
        if a.first_ncycles:solver.setOption('nCycles',a.first_ncycles)
        solver(ap,writeSolution=False)
        solver.setOption('nCycles',a.ncycles)
        if comm.rank==0:print('FIRST_SOLVE_TILT',json.dumps({'tilt_deg':a.first_solve_tilt,'failed':bool(ap.solveFailed or ap.fatalFail)}),flush=True)
        if ap.solveFailed or ap.fatalFail:
            raise RuntimeError(f'the first solve, at the start tilted by {a.first_solve_tilt:g} deg, did not converge')
        anchor_state=solver.getStates().copy();anchor_x=w
        geo.setDesignVars(x)
    if count==0 and a.first_ncycles and not tilt:solver.setOption('nCycles',a.first_ncycles)
    if a.split_stalled_solves and anchor_state is not None:path=walk(x)
    else:solver(ap,writeSolution=False)
    if count==0 and a.first_ncycles and not tilt:solver.setOption('nCycles',a.ncycles)
    f={};solver.evalFunctions(ap,f)
    cl=float(f['wing_cl_wing']);cd=float(f['wing_cd_wing']);cmz=float(f['wing_cmz_wing'])
    loss,loss_c=outlet_loss(f)
    ri,rs,rf=map(float,solver.getResNorms())
    cpmin=float(f['wing_cpmin_wing']) if a.cp_min is not None else None
    failed=bool(ap.solveFailed or ap.fatalFail or not np.isfinite([cl,cd,cmz,loss]+([cpmin] if cpmin is not None else [])).all())
    if not failed:anchor_state,anchor_x=solver.getStates().copy(),copy.deepcopy(x)
    previous_failed=failed
    if count==0:
        if failed:raise RuntimeError('Initial CFD baseline did not converge')
        settings['drag_cap']=cdcap
        settings['loss_cap']=a.loss_cap
        settings['cp_min']=a.cp_min
        settings['outlet_total_pressure_inf_Pa']=p0inf
        settings['flow_derived']={'density_kg_m3':float(ap.rho),'viscosity_Pa_s':float(ap.mu),
                                 'mach':float(ap.mach),'reynolds':float(ap.rho*V*c0/ap.mu)}
        dump(root/'settings.json',settings)
        if comm.rank==0:
            (root/'settings.txt').write_text('CFD AND OPTIMIZER SETTINGS\n\n'+json.dumps(settings,indent=2)+'\n')
    vals.update({'objective':loss/LOSS_SCALE if LOSS_OBJ else cl,
                 CAP:2-cl*c0/a.max_clc if CAP=='lift' else loss/a.loss_cap if CAP=='loss' else cd/cdcap})
    if cpmin is not None:vals['cp_min']=cpmin/a.cp_min  # both negative: <= 1 when cpmin >= the limit
    feasible=bool(not failed and np.min(vals['surface_x'])>=box['xmin']/c0-1e-6
        and np.max(vals['surface_x'])<=box['xmax']/c0+1e-6
        and np.min(vals['surface_y'])>=box['ymin']/c0-1e-6 and np.max(vals['surface_y'])<=box['ymax']/c0+1e-6
        and vals['max_thickness'][0]>=.06-1e-6 and all(vals[n]<=1+1e-5 for n in AERO)
        and all(vals[name].min()>=low/c0-1e-6 for name,(_,low) in vertex_constraints.items()))
    q=.5*ap.rho*V**2
    record={'evaluation':count,'cl':cl,'cd':cd,'cmz':cmz,
            'downforce_N_per_m':float(-cl*q*c0),'drag_N_per_m':float(cd*q*c0),
            'downforce_to_drag':float(-cl/cd),'force_down_N_extruded_span':float(-cl*q*c0*span),
            'cpt_loss_outlet':loss,'cpt_loss_outlet_flux_c_m':loss_c,
            **({'cp_min_ks':cpmin} if cpmin is not None else {}),
            'outlet_ptot_Pa':float(f['wing_ptot_outlet']),'outlet_mdot_kg_s':float(f['wing_mdot_outlet']),
            'residual_initial':ri,'residual_start':rs,'residual_final':rf,
            'residual_ratio_initial':rf/ri,'cfd_failed':failed,**path,'feasible':feasible,
            'objective_value':float(vals['objective']),'zero_deformation':unwarped(x),
            'evaluation_wall_time_s':time.perf_counter()-tic,'elapsed_wall_time_s':time.perf_counter()-start,
            'incidence_variable_deg':float(x['incidence'][0]),**metrics(pts),'design':realdict(x)}
    if comm.rank==0:
        with (root/'history.jsonl').open('a') as f:
            f.write(json.dumps(record)+'\n')
        dump(root/'latest.json',record)
        print('OPT_EVAL',json.dumps(record),flush=True)
    if count==0:save_state('baseline',record,pts,x)
    # With the loss the objective, the grid's own design (no deformation) does not count: any warp, however small,
    # adds 1-4 % to the outlet loss (every round of runs 016-022 shows it; C_L moves as it should), so the round's
    # start would beat every design the optimizer finds. It is the previous round's last design, counted there.
    if feasible and not (LOSS_OBJ and unwarped(x)) and (best is None or record['objective_value']<best['objective_value']):
        best=copy.deepcopy(record)
        save_state('best_feasible',record,pts,x)
    if failed:
        # pyOptSparse SLSQP ignores fail flags. Explicitly reject failed CFD
        # with finite penalties instead of optimizing unconverged loads.
        vals.update({'objective':1000.,**{n:1000. for n in AERO}})
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
        g.update({'objective':copy.deepcopy(zeros),**{n:copy.deepcopy(zeros) for n in AERO}})
        return g,True
    fs={};solver.evalFunctionsSens(ap,fs,evalFuncs=['cl_wing','ptot_outlet' if CAP in ('loss','lift') else 'cd_wing']
                                   +(['cpmin_wing'] if a.cp_min is not None else []))
    failed=bool(ap.adjointFailed)
    if failed:raise RuntimeError('Adjoint failed: refusing an unconverged optimization gradient')
    if LOSS_OBJ:
        g['objective']={k:-np.asarray(v)/qinf/LOSS_SCALE for k,v in fs['wing_ptot_outlet'].items()}
        g['lift']={k:-np.asarray(v)*c0/a.max_clc for k,v in fs['wing_cl_wing'].items()}
    else:
        g['objective']={k:np.asarray(v) for k,v in fs['wing_cl_wing'].items()}
        if CAP=='loss':g['loss']={k:-np.asarray(v)/qinf/a.loss_cap for k,v in fs['wing_ptot_outlet'].items()}
        else:g['drag']={k:np.asarray(v)/cdcap for k,v in fs['wing_cd_wing'].items()}
    if a.cp_min is not None:g['cp_min']={k:np.asarray(v)/a.cp_min for k,v in fs['wing_cpmin_wing'].items()}
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
                'objective':realdict(g['objective']),**{n:realdict(g[n]) for n in AERO}})+'\n')
    return g,failed

# Check independent CFD derivatives before starting the optimization.
# Small +/- perturbations use the same consistent frozen connectivity.
# Restore the baseline primal state before PSQP starts.
# Incidence: central difference. With --check-design also the shape variables, along each
# function's own adjoint gradient at three steps (at zero deformation IDWarp has a kink).
# A function passes within `limit` if
#   - its central difference agrees, or
#   - incidence: the error is under 1 % of the whole scaled gradient (the derivative is
#     near zero at the best incidence), or
#   - shape: the forward or the backward difference agrees at 2 of the 3 steps (the
#     converged response can jump on one side of the design).
# Passing at 50 % but not 5 % lets the run go on with a warning: the slope still points downhill.
def one_sided_agreement(entry,limit=.05):
    return [side for side in ['forward','backward']
            if sum(s[side+'_relative_error']<=limit for s in entry['steps_m'].values())>=2]
def derivative_check_passed(report,limit=.05,names=('objective',CAP)):
    incidence=all(report[n]['relative_error']<=limit or report[n]['error_over_scaled_gradient_norm']<=.01
                  for n in names)
    shape=all(report[n]['relative_error']<=limit or one_sided_agreement(report[n],limit)
              for n in [n+'_shape' for n in names] if n in report)
    return not report['primal_failed'] and incidence and shape
# The outlet loss is report-only (user, 2026-11-15): at low-loss designs its value jumps by 1e-6..2e-5
# between micro-steps (deterministic, likely the approximate far-field mesh warp near the outlet), so
# finite differences there cannot test its adjoint. That adjoint passed within 1 % at run_006's best
# design (run_006_levelNoTilt/cptloss_eval). Only C_L gates a loss-capped run, and with --max-clc
# (the loss the objective) only the C_L limit 'lift'.
def derivative_check_verdict(report):
    gated=(('lift',) if LOSS_OBJ else ('objective',) if CAP=='loss' else ('objective',CAP))+tuple(AERO[1:])
    verdict=('passed' if derivative_check_passed(report,names=gated) else
             'passed with a warning: within 50 % but not 5 %' if derivative_check_passed(report,.5,gated) else 'failed')
    if verdict!='failed' and CAP in ('loss','lift') and not derivative_check_passed(report,names=('objective' if LOSS_OBJ else 'loss',)):
        verdict+=' (C_L); the outlet loss did not agree, which is report-only'
    return verdict

def aerodynamic_derivative_check(x0):
    f0,fail=objfun(x0)
    base={n:float(f0[n]) for n in ['objective',*AERO]}
    g,adjfail=sensfun(x0,f0)
    if adjfail:raise RuntimeError('Baseline adjoint did not converge')
    h=1e-2
    xp=copy.deepcopy(x0);xm=copy.deepcopy(x0)
    xp['incidence']+=h;xm['incidence']-=h
    fp,pfail=objfun(xp);fm,mfail=objfun(xm)
    report={}
    for name in ['objective',*AERO]:
        fd=float((fp[name]-fm[name])/(2*h))
        adj=float(np.asarray(g[name]['incidence']).ravel()[0])
        scaled_norm=float(np.sqrt(sum(((np.asarray(g[name][k],dtype=float).ravel()/dv_scale[k])**2).sum() for k in x0)))
        report[name]={'central_difference_per_deg':fd,'adjoint_per_deg':adj,
                      'relative_error':abs(fd-adj)/max(abs(fd),abs(adj),1e-8),
                      'scaled_gradient_norm':scaled_norm,
                      'error_over_scaled_gradient_norm':abs(fd-adj)/max(scaled_norm,1e-8)}
    report['primal_failed']=pfail or mfail
    if a.check_design:
        # Separated flow is rough at 0.1 mm and solver noise shows at 0.01 mm: the best step counts.
        for name in ['objective',*AERO]:
            u={k:np.asarray(g[name][k],dtype=float).ravel() for k in shape_groups}
            norm=np.sqrt(sum((v**2).sum() for v in u.values()));u={k:v/norm for k,v in u.items()}
            adj=float(sum(np.asarray(g[name][k],dtype=float).ravel()@u[k] for k in shape_groups))
            err=lambda d:abs(d-adj)/max(abs(d),abs(adj),1e-8)
            steps={}
            for hs in [1e-5,3e-5,1e-4]:  # m of control-point motion along the unit direction
                xp=copy.deepcopy(x0);xm=copy.deepcopy(x0)
                for k in shape_groups:xp[k]=x0[k]+hs*u[k];xm[k]=x0[k]-hs*u[k]
                fp,pfail=objfun(xp);fm,mfail=objfun(xm)
                fd=float((fp[name]-fm[name])/(2*hs))
                fwd=float((fp[name]-base[name])/hs);bwd=float((base[name]-fm[name])/hs)
                steps[f'{hs:.0e}']={'central_difference_per_m':fd,'relative_error':err(fd),
                                    'forward_difference_per_m':fwd,'forward_relative_error':err(fwd),
                                    'backward_difference_per_m':bwd,'backward_relative_error':err(bwd),
                                    'plus':float(fp[name]),'minus':float(fm[name])}
                report['primal_failed']=report['primal_failed'] or pfail or mfail
            entry={'direction':'own adjoint gradient in the shape variables','adjoint_per_m':adj,'base':base[name],
                   'steps_m':steps,'relative_error':min(s['relative_error'] for s in steps.values())}
            entry['agreement']=(['central'] if entry['relative_error']<=.05 else [])+one_sided_agreement(entry)
            report[name+'_shape']=entry
    report['verdict']=derivative_check_verdict(report)
    dump(root/'source/aerodynamic_derivative_check.json',report)
    if comm.rank==0:print('AERODYNAMIC_DERIVATIVE_CHECK',report,flush=True)
    if report['verdict']=='failed':
        raise RuntimeError('Aerodynamic derivative check failed; see saved report')
    if report['verdict']!='passed' and comm.rank==0:
        print('AERODYNAMIC_DERIVATIVE_CHECK_WARNING',report['verdict']+'; the run goes on',flush=True)
    objfun(x0)

# Baseline always starts from the same NACA shape, even for a continuation,
# unless --restart-flow gives the converged flow of the restart design.
basefuncs,fail=objfun(load_design(a.restart_design) if a.restart_flow else xinitial)
if fail: raise RuntimeError('Baseline CFD failed')
if a.reuse_checks:
    report=json.loads((root/'source/aerodynamic_derivative_check.json').read_text())
    assert derivative_check_verdict(report)!='failed',report
else:
    aerodynamic_derivative_check(load_design(a.check_design) if a.check_design else xinitial)
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
for n in AERO:optprob.addCon(n,upper=1.)
optprob.addConGroup('surface_x',len(points),lower=box['xmin']/c0,upper=box['xmax']/c0)
optprob.addConGroup('surface_y',len(points),lower=box['ymin']/c0,upper=box['ymax']/c0)
optprob.addCon('max_thickness',lower=.06)
for name,(fun,low) in vertex_constraints.items():
    optprob.addConGroup(name,len(fun(points)),lower=low/c0)
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

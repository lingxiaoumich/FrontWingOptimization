"""MPI validation of prescribed moving walls using a small empty wind tunnel.

Run in the patched container. A stationary-wall control makes this test sensitive
to a missing/no-op boundary implementation; a second speed checks persistence.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from cgnsutilities.cgnsutilities import Grid, Block, Boco

parser=argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--mode', choices=['moving','stationary'], default='moving')
parser.add_argument('--derivatives', action='store_true')
args=parser.parse_args()
comm=MPI.COMM_WORLD
out=args.output.resolve()
if comm.rank==0:
    out.mkdir(parents=True,exist_ok=True)
    xx=np.linspace(-.1,.9,65)
    yy=np.r_[0., .0001*np.cumsum(1.12**np.arange(0,55))]
    yy=yy/yy[-1]*.3
    coords=np.stack(np.meshgrid(xx,yy,[0.,.02],indexing='ij'),axis=-1)
    dims=list(coords.shape[:3])
    block=Block('tunnel',dims,coords)
    for name,kind,axis,index,family in [
        ('inlet','bcfarfield',0,1,'farfield'),('outlet','bcfarfield',0,dims[0],'farfield'),
        ('top','bcsymmetryplane',1,dims[1],'symmetry'),
        ('ground','bcwallviscous',1,1,'ground'),
        ('span0','bcsymmetryplane',2,1,'symmetry'),
        ('span1','bcsymmetryplane',2,dims[2],'symmetry')]:
        points=np.array([[1,dims[0]],[1,dims[1]],[1,dims[2]]])
        points[axis]=index
        block.addBoco(Boco(name,kind,points,family))
    grid=Grid();grid.blocks=[block];grid.writeToCGNS(str(out/'tunnel.cgns'))
comm.Barrier()
options={'gridFile':str(out/'tunnel.cgns'),'outputDirectory':str(out),
         'equationType':'laminar NS','equationMode':'steady','liftIndex':2,
         'useANKSolver':True,'useNKSolver':True,'ANKSwitchTol':10.,'NKSwitchTol':1e-5,
         'ANKADPC':True,'NKADPC':True,'NKViscPC':True,
         'nCycles':12000 if args.mode=='stationary' else 1,
         'L2Convergence':1e-8,'L2ConvergenceRel':1e-12,'MGCycle':'sg',
         'useZipperMesh':False,'monitorVariables':['resrho','cl','cd'],
         'viscousSurfaceVelocities':True,'surfaceVariables':['cp','cf','yplus'],
         'writeTecplotSurfaceSolution':True,'writeVolumeSolution':True,
         'adjointL2Convergence':1e-10,'ADPC':True}
solver=ADFLOW(options=options,comm=comm)
speed=200/3.6
wall_speed=speed if args.mode=='moving' else 0.
solver.setWallVelocity([wall_speed,0,0],groupName='ground')
ap=AeroProblem(name='flat_'+args.mode,V=speed,P=101325.,T=288.15,
               alpha=0.,areaRef=.02,chordRef=1.,evalFuncs=['cl','cd'])
started=time.perf_counter()
solver(ap)
stats=solver.getWallVelocityStats('ground')
assert np.allclose(stats['min_mps'],[wall_speed,0,0],rtol=0,atol=1e-11),stats
assert np.allclose(stats['max_mps'],[wall_speed,0,0],rtol=0,atol=1e-11),stats
funcs={};solver.evalFunctions(ap,funcs)
# Read the entire physical velocity field independently from the state-vector API.
states=solver.getStates().reshape((-1,5))
reference=solver.adflow.flowvarrefstate
fluid_velocity=states[:,1:4]*float(np.sqrt(reference.pref/reference.rhoref))
velocity_error=comm.allreduce(float(np.max(np.abs(fluid_velocity-[speed,0,0]))),op=MPI.MAX)
result={'mode':args.mode,'speed_mps':speed,'wall_velocity':stats,
        'functions':{k:float(v) for k,v in funcs.items()},
        'residuals':[float(x) for x in solver.getResNorms()],
        'solve_failed':bool(ap.solveFailed),'fatal_fail':bool(ap.fatalFail),
        'seconds':time.perf_counter()-started,'mpi_ranks':comm.size,
        'max_fluid_velocity_error_from_uniform_mps':velocity_error}
assert not ap.fatalFail
if args.mode=='moving':
    assert velocity_error/speed<1e-6,velocity_error
    assert abs(funcs[ap.name+'_cd'])<1e-6,funcs
    # Uniform flow starts at roundoff residual: a relative reduction of eight
    # orders cannot be expected. Use an absolute residual check for this test.
    assert result['residuals'][2]<1e-10,result['residuals']
    result['acceptance']='absolute roundoff residual, uniform velocity, and negligible drag'
else:
    assert velocity_error/speed>.1,velocity_error
    assert funcs[ap.name+'_cd']>1e-5,funcs
    assert not ap.solveFailed
    result['acceptance']='converged stationary control with velocity deficit and positive drag'

if args.derivatives:
    rng=np.random.default_rng(1234+comm.rank)
    xv=rng.standard_normal(solver.getSpatialSize())
    jv=solver.computeJacobianVectorProductFwd(xVDot=xv,residualDeriv=True)
    seed=rng.standard_normal(len(jv))
    jtb=solver.computeJacobianVectorProductBwd(resBar=seed,xVDeriv=True)
    lhs=comm.allreduce(float(np.dot(seed,jv)),op=MPI.SUM)
    rhs=comm.allreduce(float(np.dot(xv,jtb)),op=MPI.SUM)
    error=abs(lhs-rhs)/max(abs(lhs),abs(rhs),1.)
    result['spatial_jacobian_dot_product']={'forward':lhs,'reverse':rhs,'relative_error':error}
    assert error<1e-9,result['spatial_jacobian_dot_product']
    # ADflow's FD mode perturbs volume coordinates and recomputes its residual.
    fd=solver.computeJacobianVectorProductFwd(xVDot=xv,residualDeriv=True,mode='FD',h=1e-8)
    sqerr=comm.allreduce(float(np.dot(fd-jv,fd-jv)),op=MPI.SUM)
    sqnorm=comm.allreduce(float(np.dot(jv,jv)),op=MPI.SUM)
    rel=np.sqrt(sqerr/max(sqnorm,1e-30))
    result['spatial_jacobian_fd_relative_error']=float(rel)
    assert rel<5e-3,rel
    result['wall_velocity_after_derivatives']=solver.getWallVelocityStats('ground')

# A solver-internal geometry/velocity refresh must not reset the prescribed speed.
solver.adflow.preprocessingapi.updategridvelocitiesalllevels()
refreshed=solver.getWallVelocityStats('ground')
assert np.allclose(refreshed['min_mps'],[wall_speed,0,0],atol=1e-11)
result['wall_velocity_after_refresh']=refreshed
solver.setWallVelocity([.5*speed,0,0],groupName='ground')
solver.adflow.preprocessingapi.updategridvelocitiesalllevels()
changed=solver.getWallVelocityStats('ground')
assert np.allclose(changed['min_mps'],[.5*speed,0,0],atol=1e-11)
result['wall_velocity_after_changing_setting']=changed
solver.setWallVelocity([wall_speed,0,0],groupName='ground')
solver.adflow.preprocessingapi.updategridvelocitiesalllevels()
result['validation_passed']=True
if comm.rank==0:
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

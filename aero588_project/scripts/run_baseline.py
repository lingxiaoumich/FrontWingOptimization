"""RANS baseline; forces per metre, ideal slip ground, zero incidence."""
import json
import argparse
from pathlib import Path
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW

root=Path('/project')
parser=argparse.ArgumentParser()
parser.add_argument('--restart', action='store_true')
args=parser.parse_args()
out=root/'runs/baseline_a0'
if MPI.COMM_WORLD.rank == 0:
    out.mkdir(parents=True,exist_ok=True)
MPI.COMM_WORLD.Barrier()
V=100/3.6
T=288.15
P=101325.0
rho=P/(287.058*T)
mu=1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4)
mach=V/np.sqrt(1.4*287.058*T)
re=rho*V*0.7/mu
options={
    'gridFile':str(root/'meshes/ground_baseline.cgns'),
    'outputDirectory':str(out),'equationType':'RANS',
    'liftIndex':2,'turbulenceModel':'SA',
    'useANKSolver':True,'useNKSolver':True,
    'ANKSwitchTol':10.0,'NKSwitchTol':1e-6,
    'ANKCoupledSwitchTol':1e-2,'ANKSecondOrdSwitchTol':1e-2,
    'ANKADPC':True,'NKADPC':True,'NKViscPC':True,
    'NKSubspaceSize':120,'NKPCILUFill':3,
    'nCycles':5000,'L2Convergence':1e-8,
    'L2ConvergenceRel':1e-8,'mgcycle':'sg',
    'useZipperMesh':False,'nearWallDist':0.02,
    'monitorVariables':['resrho','resturb','cl','cd'],
    'volumeVariables':['cp','mach','blank'],
    'surfaceVariables':['cp','cf','yplus'],
    'writeTecplotSurfaceSolution':True,
    'writeVolumeSolution':True,'writeSurfaceSolution':True,
}
if args.restart:
    options['restartFile']=str(out/'baseline_a0_000_vol.cgns')
solver=ADFLOW(options=options,comm=MPI.COMM_WORLD)
solver.addFamilyGroup('airfoil_only',['airfoil'])
solver.addFunction('cl','airfoil_only',name='cl_airfoil')
solver.addFunction('cd','airfoil_only',name='cd_airfoil')
solver.addSlices('z', [0.035], sliceType='absolute')
ap=AeroProblem(name='baseline_a0',mach=float(mach),reynolds=float(re),
    reynoldsLength=0.7,T=T,alpha=0.0,areaRef=0.7*0.07,chordRef=0.7,
    evalFuncs=['cl_airfoil','cd_airfoil'])
solver(ap)
funcs={}
solver.evalFunctions(ap,funcs)
initial,start,final=solver.getResNorms()
if MPI.COMM_WORLD.rank == 0:
    cl=float(funcs['baseline_a0_cl_airfoil'])
    cd=float(funcs['baseline_a0_cd_airfoil'])
    report={'model':'fully turbulent SA RANS; ideal slip ground; isolated quasi-2D airfoil',
        'alpha_deg':0.0,'speed_ms':V,'chord_m':0.7,'clearance_m':0.2,
        'temperature_K':T,'pressure_Pa':P,'density_kg_m3':rho,'viscosity_Pa_s':mu,
        'reynolds':re,'mach':float(mach),'cl_upward_positive':cl,'cd':cd,
        'downforce_to_drag':-cl/cd if cd>0 else None,
        'downforce_N_per_m':-cl*0.5*rho*V*V*0.7,
        'drag_N_per_m':cd*0.5*rho*V*V*0.7,
        'residual_initial':float(initial),'residual_start':float(start),
        'residual_final':float(final),'residual_reduction':float(final/initial),
        'solve_failed':bool(ap.solveFailed),'fatal_fail':bool(ap.fatalFail),
        'mesh_refinement_checked':False}
    report['accepted_convergence']=bool(not report['solve_failed'] and not report['fatal_fail'] and final/initial<1e-7)
    (root/'reports/baseline_a0.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

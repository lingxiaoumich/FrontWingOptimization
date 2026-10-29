"""Build and solve a quasi-2D NACA 0012 ground-effect case in MACH-Aero.

Run inside the MDOLab container with --case-root pointing to a writable mount.
Negative aerodynamic incidence means nose down: the LE-to-TE chord is rotated
counterclockwise while the freestream remains horizontal.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
parser = argparse.ArgumentParser()
parser.add_argument('phase', choices=['build', 'solve'])
parser.add_argument('--case-root', type=Path, required=True)
parser.add_argument('--restart', action='store_true')
args = parser.parse_args()
root = args.case_root
cfg = json.loads((root / 'case_config.json').read_text())
ground_model = cfg.get('ground_model', 'slip')
if ground_model not in {'slip', 'moving_no_slip', 'stationary_no_slip'}:
    raise ValueError(f'Unknown ground_model: {ground_model}')
for folder in ['geometry', 'meshes', 'reports', 'runs', 'results']:
    if comm.rank == 0:
        (root / folder).mkdir(exist_ok=True)
comm.Barrier()

def write_json(path, data):
    if comm.rank == 0:
        path.write_text(json.dumps(data, indent=2) + '\n')

def build():
    from pyhyp import pyHyp
    from cgnsutilities.cgnsutilities import Grid, Block, Boco, readGrid, combineGrids
    started = time.perf_counter()
    chord, clearance, span = [cfg[k] for k in ['chord_m', 'clearance_m', 'span_extrusion_m']]
    theta = np.linspace(0, 2*np.pi, 273)
    x = (1+np.cos(theta))/2
    yt = .6*(.2969*np.sqrt(x)-.126*x-.3516*x**2+.2843*x**3-.1015*x**4)
    y = np.where(theta <= np.pi, yt, -yt)
    closing_te = np.linspace(y[-1], y[0], 17)[1:]
    x = np.r_[x, np.ones(closing_te.size)]*chord
    y = np.r_[y, closing_te]*chord
    angle = np.deg2rad(-cfg['incidence_deg'])
    pivot = .25*chord
    xr = np.cos(angle)*(x-pivot)-np.sin(angle)*y+pivot
    yr = np.sin(angle)*(x-pivot)+np.cos(angle)*y
    offset = clearance-float(yr.min())
    yr += offset
    leading_edge = [float(pivot-pivot*np.cos(angle)), float(offset-pivot*np.sin(angle))]
    trailing_edge = [float(pivot+(chord-pivot)*np.cos(angle)),
                     float(offset+(chord-pivot)*np.sin(angle))]
    xyz = np.empty((x.size, 2, 1, 3))
    xyz[...,0], xyz[...,1] = xr[:,None,None], yr[:,None,None]
    xyz[:,0,:,2], xyz[:,1,:,2] = 0, span
    if comm.rank == 0:
        with (root/'geometry/airfoil_surface.xyz').open('w') as stream:
            stream.write(f'1\n{x.size} 2 1\n')
            for dimension in range(3):
                np.savetxt(stream, xyz[...,dimension].ravel(order='F'), fmt='%.15e')
        np.savetxt(root/'geometry/airfoil_coordinates.csv', np.column_stack([xr,yr]),
                   delimiter=',',header='x_m,y_m',comments='')
    comm.Barrier()
    hyp = pyHyp(options={
        'inputFile':str(root/'geometry/airfoil_surface.xyz'),
        'unattachedEdgesAreSymmetry':False,'outerFaceBC':'overset','autoConnect':True,
        'BC':{1:{'jLow':'zSymm','jHigh':'zSymm'}},'families':'airfoil',
        'N':cfg.get('nearfield_normal_nodes',81),'s0':cfg['first_wall_spacing_m'],'marchDist':cfg['nearfield_march_distance_m'],
        'ps0':-1.,'pGridRatio':-1.,'cMax':3.,'epsE':1.,'epsI':2.,
        'theta':3.,'volCoef':.25,'volBlend':.0001,'volSmoothIter':100})
    hyp.run()
    hyp.writeCGNS(str(root/'meshes/airfoil_near.cgns'))
    if comm.rank == 0:
        def extend(start,stop,step,ratio=1.08):
            values=[start]
            direction=np.sign(stop-start)
            while abs(stop-values[-1])>step:
                values.append(values[-1]+direction*step)
                step*=ratio
            values.append(stop)
            return np.array(values)
        xx=np.r_[extend(-.15,-10*chord,.005)[::-1][:-1],
                 np.linspace(-.15,.85,201),extend(.85,15*chord,.005)[1:]]
        if ground_model == 'slip':
            yy=np.r_[np.linspace(0,.35,141),extend(.35,10*chord,.0025)[1:]]
        else:
            # Resolve viscous road shear; the old slip grid had 2.5 mm cells at the ground.
            ground_y=[0.]
            step=cfg.get('ground_first_wall_spacing_m', cfg['first_wall_spacing_m'])
            while ground_y[-1]+step < .025:
                ground_y.append(ground_y[-1]+step)
                step=min(1.12*step,.0025)
            ground_y.append(.025)
            yy=np.r_[ground_y,np.linspace(.025,.35,131)[1:],extend(.35,10*chord,.0025)[1:]]
        coordinates=np.stack(np.meshgrid(xx,yy,[0,span],indexing='ij'),axis=-1)
        dims=list(coordinates.shape[:3])
        block=Block('background',dims,coordinates)
        def boundary(name,kind,axis,index,family):
            point_range=np.array([[1,dims[0]],[1,dims[1]],[1,dims[2]]])
            point_range[axis]=index
            block.addBoco(Boco(name,kind,point_range,family))
        boundary('inlet','bcfarfield',0,1,'farfield')
        boundary('outlet','bcfarfield',0,dims[0],'farfield')
        boundary('top','bcfarfield',1,dims[1],'farfield')
        ground_bc = 'bcsymmetryplane' if ground_model == 'slip' else 'bcwallviscous'
        boundary('ground',ground_bc,1,1,'ground')
        boundary('span0','bcsymmetryplane',2,1,'symmetry')
        boundary('span1','bcsymmetryplane',2,2,'symmetry')
        background=Grid()
        background.blocks=[block]
        grid=combineGrids([readGrid(str(root/'meshes/airfoil_near.cgns')),background])
        grid.connect()
        grid.writeToCGNS(str(root/'meshes/ground.cgns'))
        checks=[]
        for b in grid.blocks:
            p=b.coords
            a=p[1:,:-1,:-1]-p[:-1,:-1,:-1]
            d=p[:-1,1:,:-1]-p[:-1,:-1,:-1]
            e=p[:-1,:-1,1:]-p[:-1,:-1,:-1]
            jac=np.einsum('...i,...i->...',np.cross(a,d),e)
            checks.append({'block':b.name,'vertices':list(p.shape[:3]),
                           'cells':int(jac.size),'min_corner_jacobian_m3':float(jac.min()),
                           'nonpositive_corners':int((jac<=0).sum()),
                           'lowest_grid_y_m':float(p[...,1].min())})
        assert all(b['nonpositive_corners']==0 for b in checks), 'Invalid corner Jacobian'
        assert float(yr.min()) >= clearance-1e-12, 'Insufficient surface clearance'
        assert float(grid.blocks[0].coords[...,1].min()) > 0, 'Near grid crosses ground'
        info={'case':cfg,'leading_edge_m':leading_edge,'trailing_edge_midpoint_m':trailing_edge,
              'LE_to_TE_chord_angle_deg':float(np.rad2deg(np.arctan2(
                  trailing_edge[1]-leading_edge[1],trailing_edge[0]-leading_edge[0]))),
              'minimum_surface_clearance_m':float(yr.min()),'blocks':checks,
              'total_cells_before_blanking':sum(c['cells'] for c in checks),
              'ground_model':ground_model,'ground_boundary_type':ground_bc,
              'ground_first_cell_height_m':float(yy[1]-yy[0]),
              'mesh_build_wall_time_s':time.perf_counter()-started}
        write_json(root/'reports/mesh_check.json',info)
        print(json.dumps(info,indent=2))

def solve():
    from baseclasses import AeroProblem
    from adflow import ADFLOW
    chord, span = cfg['chord_m'], cfg['span_extrusion_m']
    velocity=cfg['speed_kmh']/3.6
    temperature=cfg['temperature_K']
    pressure=cfg['pressure_Pa']
    density=pressure/(287.058*temperature)
    viscosity=1.716e-5*(temperature/273.15)**1.5*(273.15+110.4)/(temperature+110.4)
    mach=velocity/np.sqrt(1.4*287.058*temperature)
    reynolds=density*velocity*chord/viscosity
    options={
        'gridFile':str(root/'meshes/ground.cgns'),'outputDirectory':str(root/'runs'),
        'equationType':'RANS','equationMode':'steady','liftIndex':2,'turbulenceModel':'SA',
        'turbulenceOrder':'first order','discretization':'central plus scalar dissipation',
        'useWallFunctions':False,'useANKSolver':True,'useNKSolver':True,
        'ANKSwitchTol':10.,'NKSwitchTol':1e-6,
        'ANKCoupledSwitchTol':cfg.get('ank_coupled_switch_tol',1e-2),
        'ANKSecondOrdSwitchTol':cfg.get('ank_second_order_switch_tol',1e-2),
        'ANKADPC':True,'NKADPC':True,'NKViscPC':True,
        'NKSubspaceSize':cfg.get('nk_subspace_size',120),
        'NKPCILUFill':cfg.get('nk_pc_ilu_fill',3),
        'NKLinearSolveTol':cfg.get('nk_linear_solve_tol',.3),
        'NKUseEW':cfg.get('nk_use_ew',True),'turbResScale':10000.,
        'nCycles':cfg.get('n_cycles',8000),
        'L2Convergence':1e-8,'L2ConvergenceRel':1e-8,'MGCycle':'sg',
        'useZipperMesh':False,'nearWallDist':.01,
        'monitorVariables':['resrho','resturb','cl','cd'],
        'volumeVariables':['cp','mach','blank'],
        'surfaceVariables':['cp','cf','yplus'],
        'writeTecplotSurfaceSolution':True,'writeVolumeSolution':True,'writeSurfaceSolution':True}
    if args.restart:
        options['restartFile']=str(root/'runs'/f'{cfg["case_id"]}_000_vol.cgns')
    started=time.perf_counter()
    solver=ADFLOW(options=options,comm=comm)
    solver.addFamilyGroup('airfoil_only',['airfoil'])
    solver.addFunction('cl','airfoil_only',name='cl_airfoil')
    solver.addFunction('cd','airfoil_only',name='cd_airfoil')
    solver.addSlices('z',[span/2],sliceType='absolute',groupName='airfoil_only')
    if ground_model == 'moving_no_slip':
        if not hasattr(solver, 'setWallVelocity'):
            raise RuntimeError('Moving ground requires the aero588/adflow-moving-wall:2.12.2 image')
        solver.setWallVelocity([velocity,0.,0.],groupName='ground')
    flow_inputs=({'mach':float(mach),'reynolds':float(reynolds),'reynoldsLength':chord}
                 if ground_model=='slip' else {'V':velocity,'P':pressure})
    ap=AeroProblem(name=cfg['case_id'],**flow_inputs,T=temperature,alpha=0.,areaRef=chord*span,
                   chordRef=chord,evalFuncs=['cl_airfoil','cd_airfoil'])
    if ground_model!='slip':
        density,viscosity,mach=float(ap.rho),float(ap.mu),float(ap.mach)
        reynolds=density*velocity*chord/viscosity
    solver(ap)
    funcs={}
    solver.evalFunctions(ap,funcs)
    initial,start,final=solver.getResNorms()
    wall_stats = None
    if ground_model == 'moving_no_slip':
        wall_stats={'ground':solver.getWallVelocityStats('ground'),
                    'airfoil':solver.getWallVelocityStats('airfoil')}
    if comm.rank==0:
        cl=float(funcs[f'{cfg["case_id"]}_cl_airfoil'])
        cd=float(funcs[f'{cfg["case_id"]}_cd_airfoil'])
        failed=bool(ap.solveFailed) or bool(ap.fatalFail)
        reduction=float(final/initial)
        ground_description={'slip':'ideal slip ground','moving_no_slip':'moving no-slip ground',
                            'stationary_no_slip':'stationary no-slip ground'}[ground_model]
        result={**cfg,'model':f'fully turbulent SA RANS; {ground_description}; isolated quasi-2D airfoil',
                'ground_model':ground_model,'wall_velocity_stats':wall_stats,
                'speed_ms':velocity,'alpha_deg':cfg['incidence_deg'],'freestream_alpha_deg':0.,
                'density_kg_m3':density,'viscosity_Pa_s':viscosity,
                'reynolds':float(reynolds),'mach':float(mach),
                'cl_upward_positive':cl,'cd':cd,'downforce_to_drag':-cl/cd if cd>0 else None,
                'downforce_N_per_m':-cl*.5*density*velocity**2*chord,
                'drag_N_per_m':cd*.5*density*velocity**2*chord,
                'residual_initial':float(initial),'residual_start':float(start),
                'residual_final':float(final),'residual_reduction':reduction,
                'solve_failed':bool(ap.solveFailed),'fatal_fail':bool(ap.fatalFail),
                'accepted_convergence':bool(not failed and reduction<1e-8),
                'mesh_refinement_checked':False,'mpi_ranks':comm.size,
                'solve_including_setup_and_output_wall_time_s':time.perf_counter()-started,
                'restarted':bool(args.restart)}
        write_json(root/'results/cfd_result.json',result)
        print(json.dumps(result,indent=2))

if args.phase=='build':
    build()
else:
    solve()

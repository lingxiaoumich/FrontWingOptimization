"""NACA 0012 RANS overset grid with ideal slip ground, in metres."""
import json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from pyhyp import pyHyp
from cgnsutilities.cgnsutilities import Grid, Block, Boco, readGrid, combineGrids

root = Path('/project')
for folder in ['meshes', 'geometry', 'reports', 'runs']:
    (root / folder).mkdir(exist_ok=True)
c, h, span = 0.7, 0.2, 0.07
# Standard finite-thickness NACA trailing edge, closed by a small vertical face.
theta = np.linspace(0, 2*np.pi, 273)
x = (1 + np.cos(theta))/2
yt = 0.6*(0.2969*np.sqrt(x)-0.126*x-0.3516*x*x+0.2843*x**3-0.1015*x**4)
y = np.where(theta <= np.pi, yt, -yt)
te = np.linspace(y[-1], y[0], 17)[1:]
x = np.r_[x, np.ones(te.size)]*c
y = np.r_[y, te]*c
y += h-y.min()
xyz = np.empty((x.size, 2, 1, 3))
xyz[..., 0] = x[:, None, None]
xyz[..., 1] = y[:, None, None]
xyz[:, 0, :, 2] = 0
xyz[:, 1, :, 2] = span
if MPI.COMM_WORLD.rank == 0:
    with (root/'geometry/airfoil_surface.xyz').open('w') as f:
        f.write(f'1\n{x.size} 2 1\n')
        for dim in range(3):
            np.savetxt(f, xyz[...,dim].ravel(order='F'), fmt='%.15e')
MPI.COMM_WORLD.Barrier()
hyp = pyHyp(options={'inputFile': str(root/'geometry/airfoil_surface.xyz'),
    'unattachedEdgesAreSymmetry':False, 'outerFaceBC':'overset',
    'autoConnect':True, 'BC':{1:{'jLow':'zSymm','jHigh':'zSymm'}},
    'families':'airfoil', 'N':81, 's0':5e-6, 'marchDist':0.11,
    'ps0':-1.0,'pGridRatio':-1.0,'cMax':3.0,'epsE':1.0,'epsI':2.0,
    'theta':3.0,'volCoef':0.25,'volBlend':0.0001,'volSmoothIter':100})
hyp.run()
hyp.writeCGNS(str(root/'meshes/airfoil_near.cgns'))
if MPI.COMM_WORLD.rank == 0:
    def extend(start, stop, step, ratio):
        out=[start]
        direction=np.sign(stop-start)
        while abs(stop-out[-1])>step:
            out.append(out[-1]+direction*step)
            step*=ratio
        out.append(stop)
        return np.array(out)
    xx=np.r_[extend(-0.35,-7,0.01,1.08)[::-1][:-1], np.linspace(-0.35,1.4,176), extend(1.4,10.5,0.01,1.08)[1:]]
    yy=np.r_[np.linspace(0,0.65,131),extend(0.65,7,0.005,1.08)[1:]]
    coords=np.stack(np.meshgrid(xx,yy,[0,span],indexing='ij'),axis=-1)
    dims=list(coords.shape[:3])
    block=Block('background',dims,coords)
    def bc(name,kind,axis,idx,family):
        pr=np.array([[1,dims[0]],[1,dims[1]],[1,dims[2]]])
        pr[axis]=idx
        block.addBoco(Boco(name,kind,pr,family))
    bc('inlet','bcfarfield',0,1,'farfield')
    bc('outlet','bcfarfield',0,dims[0],'farfield')
    bc('top','bcfarfield',1,dims[1],'farfield')
    bc('ground','bcsymmetryplane',1,1,'ground')
    bc('span0','bcsymmetryplane',2,1,'symmetry')
    bc('span1','bcsymmetryplane',2,2,'symmetry')
    background=Grid()
    background.blocks=[block]
    grid=combineGrids([readGrid(str(root/'meshes/airfoil_near.cgns')),background])
    grid.connect()
    grid.writeToCGNS(str(root/'meshes/ground_baseline.cgns'))
    checks=[]
    for b in grid.blocks:
        p=b.coords
        # Signed corner Jacobian; zero/negative values reject the mesh.
        a=p[1:,:-1,:-1]-p[:-1,:-1,:-1]
        d=p[:-1,1:,:-1]-p[:-1,:-1,:-1]
        e=p[:-1,:-1,1:]-p[:-1,:-1,:-1]
        jac=np.einsum('...i,...i->...',np.cross(a,d),e)
        checks.append({'block':b.name,'cells':int(jac.size),'min_corner_jacobian_m3':float(jac.min()),'nonpositive_corners':int((jac<=0).sum())})
    report={'ground_model':'ideal slip / symmetry; no road shear', 'chord_m':c,'lowest_point_m':float(y.min()),'span_extrusion_m':span,'first_layer_m':5e-6,'trailing_edge':'finite thickness, standard -0.1015 coefficient; differs from earlier sharp-TE coordinate file','trailing_edge_thickness_m':float(2*0.6*c*(0.2969-0.126-0.3516+0.2843-0.1015)), 'blocks':checks}
    (root/'reports/mesh_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    assert all(b['nonpositive_corners']==0 for b in checks), 'Invalid mesh Jacobians'

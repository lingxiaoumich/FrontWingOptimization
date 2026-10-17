from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cgnsutilities.cgnsutilities import readGrid

root=Path('/project')
grid=readGrid(str(root/'meshes/ground_baseline.cgns'))
fig,ax=plt.subplots(figsize=(11,6),layout='constrained')
for b in grid.blocks:
    p=b.coords[:,0,:,:] if 'near' in b.name else b.coords[:,:,0,:]
    stride=3 if 'near' in b.name else 4
    for i in range(0,p.shape[0],stride):
        ax.plot(p[i,:,0],p[i,:,1],color='#91a7b6',linewidth=0.35)
    for j in range(0,p.shape[1],stride):
        ax.plot(p[:,j,0],p[:,j,1],color='#91a7b6',linewidth=0.35)
    if 'near' in b.name:
        ax.fill(p[:,0,0],p[:,0,1],color='#243e55',zorder=3)
ax.axhline(0,color='#222222',linewidth=2)
ax.set(xlim=(-0.2,0.95),ylim=(-0.04,0.62),xlabel='x (m)',ylabel='Height above ground (m)',title='NACA 0012 baseline: 0.70 m chord, 0.20 m clearance\nOverset RANS mesh; ideal slip ground')
ax.set_aspect('equal')
fig.savefig(root/'reports/mesh_preview.png',dpi=180)

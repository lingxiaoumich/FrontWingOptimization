"""XDSM of one optimization (Lambe & Martins 2012 convention, as in Martins & Ning, Engineering Design Optimization):
components on the diagonal, data in parallelograms (row = from, column = to), grey lines for data, thin black lines
for the order of steps. Writes optimizations/comparisons/report_20261122/xdsm.png.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

OUT = Path(__file__).resolve().parents[1]/'optimizations/comparisons/report_20261122/xdsm.png'
DX, DY, W, H = 3.0, 1.3, 2.6, .8
pos = lambda r, c: (c*DX, -r*DY)
ink, grey, blue = '#2b2b2b', '#c9c9c9', '#1f62bd'

fig, ax = plt.subplots(figsize=(12.5, 7.2))
ax.set_aspect('equal'); ax.axis('off')
plt.rcParams.update({'font.family': 'DejaVu Sans'})

comps = {0: ('0, 5→1: Optimizer\n(PSQP)', True), 1: ('1: Geometry\n(pyGeo FFD)', False), 2: ('2: Mesh warping\n(IDWarp)', False),
         3: ('3: Flow solve\n(ADflow, RANS SA)', False), 4: ('4: Adjoint\n(ADflow, pyGeo)', False)}
data = {(-1, 0): 'x⁰: start design', (0, -1): 'x*: best feasible\ndesign',
        (0, 1): '1: x\n(40 moves)', (1, 0): '5: geometric\nconstraints',
        (1, 2): '2: surface points', (2, 3): '3: volume mesh', (3, 0): '5: $-C_Lc$, $C_Dc$,\nloss, $C_p$ min',
        (3, 4): '4: flow state', (4, 0): '5: gradients of\nall functions'}

# data lines: each row from its component to its outputs, each column from its inputs down/up to its component
for (r, c) in data:
    x0, y0 = pos(r, c)
    if r >= 0 and c >= 0 or c == -1:
        xr, yr = pos(r, r)
        ax.plot([x0, xr], [y0, yr], c=grey, lw=9, solid_capstyle='butt', zorder=1)
    if c >= 0:
        xc, yc = pos(c, c)
        ax.plot([x0, xc], [y0, yc], c=grey, lw=9, solid_capstyle='butt', zorder=1)
# order of steps: right along a row, then down into the next component; from the adjoint back left and up to the optimizer
off = .22
def path(pts):
    xs, ys = zip(*pts)
    ax.plot(xs[:-1], ys[:-1], c=ink, lw=1.1, zorder=3)
    ax.annotate('', pts[-1], pts[-2], arrowprops=dict(arrowstyle='-|>', color=ink, lw=1.1, shrinkA=0, shrinkB=0), zorder=3)
for k in range(4):
    (xa, ya), (xb, yb) = pos(k, k), pos(k+1, k+1)
    path([(xa+W/2, ya+off), (xb-off, ya+off), (xb-off, yb+H/2)])
(xa, ya), (xb, yb) = pos(4, 4), pos(0, 0)
path([(xa-W/2, ya-off), (xb+off, ya-off), (xb+off, yb-H/2)])
for i, (txt, opt) in comps.items():
    x, y = pos(i, i)
    ax.add_patch(FancyBboxPatch((x-W/2, y-H/2), W, H, boxstyle='round,pad=0,rounding_size=.3' if opt else 'square,pad=0',
                                fc='#dbe9fb' if opt else '#eef5fd', ec=blue, lw=1.4, zorder=4))
    ax.text(x, y, txt, ha='center', va='center', fontsize=10.5, zorder=5)
for (r, c), txt in data.items():
    x, y = pos(r, c); s = .16
    ax.add_patch(Polygon([(x-W/2+s, y+H/2), (x+W/2, y+H/2), (x+W/2-s, y-H/2), (x-W/2, y-H/2)], fc='#f3f3f3', ec='#8a8a8a', lw=1, zorder=4))
    ax.text(x, y, txt, ha='center', va='center', fontsize=9.5, zorder=5)
ax.set(xlim=(-DX-W/2-.1, 4*DX+W/2+.1), ylim=(-4*DY-H/2-.9, DY+H/2+.1))
ax.text(-DX-W/2, -4*DY-H/2-.45, 'Grey: data (row = from, column = to). Black: order of steps. Not drawn: when the mesh wears out, a new mesh is '
        'grown around the latest design and the optimizer restarts from it.', fontsize=9.5, color='#6b6b6b', va='top')
fig.savefig(OUT, dpi=150, bbox_inches='tight'); print(OUT)

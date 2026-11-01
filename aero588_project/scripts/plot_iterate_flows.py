"""Plot the flow at the PSQP iterates whose full solutions were saved (--save-iterate-flows).

Each source is a chain folder (all its stages' iterates/eval_NNN) or one stage folder,
taken in the order given; consecutive iterates with the same forces (a stage restart,
or a stalled stage) are kept once. For a sequence of iterates this writes, under <out>/plots:
  cpt_iterates_<name>   total-pressure coefficient around the section at eight iterates
  cp_iterates_<name>    surface Cp at those iterates, the section where it pinches, and
                        the true surface gap against the 2 mm pair measure for every iterate
  cpt_iterates_<name>.gif  every iterate, one frame each
and <out>/iterate_flows_<name>.json with the per-iterate numbers.

Cpt = 1 - (p0,inf - p0)/q_inf (1 in the freestream, about Cp in dead air). The
recirculation zone is where u < 0, outlined and reported as an area over computed
(iblank 1) cells of both grids. Reversed wall flow is the sign of the first-cell velocity
along the surface, away from the stagnation point; the blunt trailing-edge base
(vertices 272-288) is left out. The true gap is the smallest distance
from a vertex k = 15..124 (or its mirror 272-k) to the opposite surface polyline, which
is what the chord-normal pair constraint is meant to bound. Run with the MACH-Aero env:
  python scripts/plot_iterate_flows.py --out optimizations/continuation_summary_20261031 \
      --name remesh2_xmax1.0 optimizations/continuation_tmin2mm_xmax5.0_20261031/stage_02 \
      optimizations/continuation_tmin2mm_xmax5.0_20261031/stage_03 optimizations/continuation_remesh2_xmax1.0_20261031
"""
import argparse
import ctypes as ct
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from PIL import Image

p = argparse.ArgumentParser()
p.add_argument('sources', type=Path, nargs='+')
p.add_argument('--out', type=Path, required=True)
p.add_argument('--name', required=True)
p.add_argument('--panels', type=int, default=8)
a = p.parse_args()
out = a.out.resolve(); (out/'plots').mkdir(parents=True, exist_ok=True)

lib = ct.CDLL(str(Path(sys.prefix)/'lib/libcgns.so')); cgsize = ct.c_int64
def read_zones(path):
    """Vertex x, y and cell-centred fields of each zone's first z plane."""
    fn = ct.c_int()
    if lib.cg_open(str(path).encode(), 0, ct.byref(fn)): raise RuntimeError(path)
    zones = {}
    for z in (1, 2):
        name = ct.create_string_buffer(64); size = (cgsize*9)(); lib.cg_zone_read(fn, 1, z, name, size)
        n, m = list(size[:3]), list(size[3:6])
        near = b'airfoil' in name.value
        plane = (lambda v: v[:, 0, :]) if near else (lambda v: v[:, :, 0])  # near grid is (i, z, j)
        def read(field, dims, coord):
            buf = np.empty(int(np.prod(dims))); lo = (cgsize*3)(1, 1, 1); hi = (cgsize*3)(*dims)
            args = (fn, 1, z, field, 4, lo, hi, buf.ctypes.data_as(ct.c_void_p))
            code = lib.cg_coord_read(*args) if coord else lib.cg_field_read(fn, 1, z, 1, *args[3:])
            if code: raise RuntimeError(field)
            return plane(buf.reshape(dims, order='F'))
        f = {k: read(k.encode(), n, True) for k in ['CoordinateX', 'CoordinateY']}
        f.update({k: read(k.encode(), m, False) for k in ['CoefPressure', 'Mach', 'VelocityX', 'VelocityY', 'Iblank']})
        zones['near' if near else 'background'] = f
    lib.cg_close(fn)
    return zones

def seg_dist(q, A, B):
    AB = B-A; t = np.clip(((q-A)*AB).sum(1)/(AB*AB).sum(1), 0, 1)
    return np.linalg.norm(A+t[:, None]*AB-q, axis=1)
def true_gap(xy):
    best = (np.inf, 0, 0)
    for k in range(15, 125):
        for v, side, off in [(k, xy[136:273], 136), (272-k, xy[:137], 0)]:
            d = seg_dist(xy[v], side[:-1], side[1:]); j = int(d.argmin())
            if d[j] < best[0]: best = (float(d[j]), v, off+j)
    return best

def slices(path):
    rows = []
    for line in path.read_text().splitlines():
        v = line.split()
        if len(v) == 9:
            try: rows.append([float(x) for x in v])
            except ValueError: pass
    return np.array(rows)

def iterate_dirs(src):
    stages = sorted(s for s in src.glob('stage_*') if s.is_dir()) or [src]
    return [d for s in stages for d in sorted((s/'iterates').glob('eval_*')) if (d/'flow_000_vol.cgns').exists()]

seq = []
for src in a.sources:
    for d in iterate_dirs(src.resolve()):
        r = json.loads((d/'result.json').read_text())
        if seq and abs(r['cl']-seq[-1]['r']['cl']) < 1e-4 and abs(r['cd']-seq[-1]['r']['cd']) < 1e-4: continue
        seq.append({'dir': d, 'r': r})
settings = json.loads((seq[0]['dir'].parents[1]/'settings.json').read_text())
mach = settings['flow_derived']['mach']; g = 1.4
p0_inf = (1+.2*mach**2)**3.5  # over p_inf; q_inf/p_inf = g M^2/2

def load(it):
    zones = read_zones(it['dir']/'flow_000_vol.cgns')
    for f in zones.values():
        pp = 1+f['CoefPressure']*g*mach**2/2
        f['Cpt'] = 1-(p0_inf-pp*(1+.2*f['Mach']**2)**3.5)/(g*mach**2/2)
    xy = np.loadtxt(it['dir']/'surface_geometry.csv', delimiter=',', skiprows=1)[:, :2]
    near = zones['near']; wall = np.c_[near['CoordinateX'][:, 0], near['CoordinateY'][:, 0]]
    assert np.abs(wall[:len(xy)]-xy).max() < 1e-9, 'near-grid wall is not the surface vertex order'
    t = np.diff(wall[:273], axis=0); t /= np.linalg.norm(t, axis=1)[:, None]
    ut = near['VelocityX'][:272, 0]*t[:, 0]+near['VelocityY'][:272, 0]*t[:, 1]
    stag = int(np.argmin(np.abs(ut[100:172])))+100  # flow runs away from the stagnation point near the leading edge
    reversed_ = np.r_[ut[:stag] > 0, ut[stag:] < 0]
    reversed_[stag-2:stag+3] = False
    back = 0.
    for f in zones.values():
        x, y = f['CoordinateX'], f['CoordinateY']
        area = .5*np.abs((x[1:, 1:]-x[:-1, :-1])*(y[:-1, 1:]-y[1:, :-1])-(x[:-1, 1:]-x[1:, :-1])*(y[1:, 1:]-y[:-1, :-1]))
        back += float(area[(f['Iblank'] == 1) & (f['VelocityX'] < 0)].sum())
    gap, v, w = true_gap(xy)
    s = slices(it['dir']/'flow_000_slices.dat'); s = s[np.abs(s[:, 2]-s[0, 2]) < 1e-9]
    idx = np.argmin(((s[:, None, :2]-xy[None, :273])**2).sum(-1), axis=1)
    return zones, xy, {'reversed_cells': np.flatnonzero(reversed_), 'gap': gap, 'gap_vertices': (v, w), 'reversed_area': back, 'slice': s, 'slice_index': idx}

def centers(x): return .25*(x[:-1, :-1]+x[1:, :-1]+x[1:, 1:]+x[:-1, 1:])
def tri(f, near):
    x, y = centers(f['CoordinateX']), centers(f['CoordinateY'])
    ids = np.arange(x.size).reshape(x.shape)
    quads = [[ids[:-1, :-1], ids[1:, :-1], ids[1:, 1:], ids[:-1, 1:]]]
    if near: quads.append([ids[-1, :-1], ids[0, :-1], ids[0, 1:], ids[-1, 1:]])  # O-grid seam
    t = np.concatenate([np.r_[np.c_[q[0].ravel(), q[1].ravel(), q[2].ravel()], np.c_[q[0].ravel(), q[2].ravel(), q[3].ravel()]] for q in quads])
    live = np.isin(f['Iblank'], [1, -1]).ravel()
    return mtri.Triangulation(x.ravel(), y.ravel(), t[live[t].all(axis=1)])

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'font.family': 'DejaVu Sans'})
ink, muted, grid = '#2b2b2b', '#6b6b6b', '#e4e4e4'
ramp = ['#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
order_cmap = LinearSegmentedColormap.from_list('early_to_late', ramp)
cpt_cmap = LinearSegmentedColormap.from_list('cpt', ['#0d366b', '#2a78d6', '#86b6ef', '#f4f8fd'])  # loss dark -> freestream light
reverse_c, pinch_c = '#eb6834', '#d03b3b'
levels = np.linspace(-1.2, 1.02, 75); cnorm = Normalize(-1.2, 1.02); sm = ScalarMappable(norm=cnorm, cmap=cpt_cmap)
box = dict(xmin=-.125, xmax=.625, ymin=.05, ymax=.35)
def save(fig, name):
    fig.savefig(out/'plots'/(name+'.png'), dpi=170, bbox_inches='tight')
    fig.savefig(out/'plots'/(name+'.svg'), bbox_inches='tight'); plt.close(fig)
def title(n, it, extra):
    r = it['r']
    return (f"iterate {n+1}/{len(seq)}{' · remeshed grid' if 'remesh' in str(it['dir']) else ' · warped grid'}\n"
            f"$C_L c$ {r['cl']*.5:.3f} m, $C_D c$ {r['cd']*.5:.3f} m, gap {extra['gap']*1e3:.2f} mm, u < 0 over {extra['reversed_area']*1e4:.0f} cm²")
def field(ax, zones, xy, extra, xlim=(-.2, 1.05), ylim=(0, .42)):
    for tag in ['background', 'near']:
        f = zones[tag]
        ax.tricontourf(tri(f, tag == 'near'), np.clip(f['Cpt'], -50, 2).ravel(), levels=levels, cmap=cpt_cmap,
                       norm=cnorm, extend='both', zorder=2, rasterized=True)
        ax.tricontour(tri(f, tag == 'near'), f['VelocityX'].ravel(), levels=[0], colors=reverse_c, linewidths=1.1, zorder=3)
    ax.fill(xy[:, 0], xy[:, 1], color='#9a9a9a', lw=.6, ec=ink, zorder=3)
    ax.add_patch(Rectangle((box['xmin'], box['ymin']), box['xmax']-box['xmin'], box['ymax']-box['ymin'], fill=False, ls='--', ec=ink, lw=.8, zorder=4))
    ax.axhline(0, c=ink, lw=2, zorder=4)
    if extra['gap'] < 1e-3:
        v = xy[extra['gap_vertices'][0]]
        ax.scatter(*v, s=90, facecolor='none', edgecolor=pinch_c, lw=1.4, zorder=6)
    ax.set(xlim=xlim, ylim=ylim, xlabel='x (m)', ylabel='height (m)'); ax.set_aspect('equal')

# All iterates: numbers, and the GIF frames
data = []; frames = []
for n, it in enumerate(seq):
    zones, xy, extra = load(it)
    it.update(xy=xy, extra=extra)
    rc = extra['reversed_cells']
    data.append({'iterate': str(it['dir'].relative_to(out.parent)), 'cl_c_m': it['r']['cl']*.5, 'cd_c_m': it['r']['cd']*.5,
                 'true_gap_mm': extra['gap']*1e3, 'gap_vertices': list(extra['gap_vertices']),
                 'pair_min_thickness_mm': it['r']['minimum_local_thickness_m']*1e3,
                 'reversed_area_cm2': extra['reversed_area']*1e4, 'reversed_wall_cells': int(len(rc)), 'reversed_x_m': [float(xy[rc, 0].min()), float(xy[rc, 0].max())] if len(rc) else None,
                 'min_cpt': float(min(z['Cpt'][np.isin(z['Iblank'], [1, -1])].min() for z in zones.values()))})
    fig, ax = plt.subplots(figsize=(9, 4.1))
    field(ax, zones, xy, extra); ax.set_title(title(n, it, extra), fontsize=10, loc='left')
    fig.colorbar(sm, ax=ax, shrink=.8, pad=.02, label='$C_{pt}$', extend='both')
    fig.canvas.draw(); frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')); plt.close(fig)
    print(f"{n+1:3d} {data[-1]['iterate']} CLc {data[-1]['cl_c_m']:.4f} gap {data[-1]['true_gap_mm']:.4f} mm reversed {len(rc)}", flush=True)
frames[0].save(out/'plots'/f'cpt_iterates_{a.name}.gif', save_all=True, append_images=frames[1:],
               duration=[450]*(len(frames)-1)+[2500], loop=0, optimize=True)
(out/f'iterate_flows_{a.name}.json').write_text(json.dumps(data, indent=1)+'\n')

# Panels: the first and last iterates and the ones nearest evenly spaced C_L c in between
clc = np.array([d['cl_c_m'] for d in data])
pick = sorted({int(np.argmin(np.abs(clc-v))) for v in np.linspace(clc[0], clc[-1], a.panels)})
fig, axs = plt.subplots((len(pick)+1)//2, 2, figsize=(14, 2.75*((len(pick)+1)//2)+.6), constrained_layout=True)
for ax, n in zip(axs.flat, pick):
    zones, _, _ = load(seq[n]); field(ax, zones, seq[n]['xy'], seq[n]['extra']); ax.set_title(title(n, seq[n], seq[n]['extra']), fontsize=9, loc='left')
for ax in axs.flat[len(pick):]: ax.axis('off')
fig.colorbar(sm, ax=axs, shrink=.6, pad=.01, extend='both',
             label='$C_{pt} = 1-(p_{0,\\infty}-p_0)/q_\\infty$')
fig.suptitle('Total-pressure coefficient at saved iterates — orange line: u = 0 (encloses reversed flow); red circle: surfaces within 1 mm', fontsize=11)
save(fig, f'cpt_iterates_{a.name}')

# Surface Cp, the pinch, and the gap history
colors = order_cmap(Normalize(0, len(seq)-1)(pick))
fig = plt.figure(figsize=(15, 8.6)); gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=.34, wspace=.28)
ax = fig.add_subplot(gs[0, :])
for n, c in zip(pick, colors):
    e = seq[n]['extra']; s, k = e['slice'], e['slice_index']
    for side, ls in [(k <= 136, '-'), (k >= 136, '--')]:
        o = np.argsort(k[side]); ax.plot(s[side][o, 0], s[side][o, 6], ls, c=c, lw=1.4)
    ax.plot([], [], c=c, lw=2, label=f"it {n+1}: {data[n]['cl_c_m']:.3f} m")
ax.plot([], [], c=ink, ls='-', label='upper side (vertices 0–136)'); ax.plot([], [], c=ink, ls='--', label='ground side (136–272)')
ax.invert_yaxis(); ax.grid(c=grid); ax.set(xlabel='x (m)', ylabel='$C_p$', title='Surface pressure at the panel iterates (sections move and rotate; plotted at actual x)')
ax.legend(ncol=2, fontsize=8, frameon=False, loc='upper right')
ax = fig.add_subplot(gs[1, 0])
w = seq[-1]['extra']['gap_vertices']; xy = seq[-1]['xy']; centre = xy[w[0]]
for side, ls, name in [(slice(0, 137), '-', 'upper side'), (slice(136, 273), '--', 'ground side')]:
    ax.plot(xy[side, 0], xy[side, 1], ls, c=ramp[-1], lw=1.3, marker='.', ms=4, label=name)
ax.scatter(*centre, s=110, facecolor='none', edgecolor=pinch_c, lw=1.4, zorder=5)
ax.annotate(f"vertex {w[0]} lies on the other surface\n(segment {w[1]}–{w[1]+1}); its pair {272-w[0]}\nreads {data[-1]['pair_min_thickness_mm']:.2f} mm", centre,
            xytext=(.03, .97), textcoords='axes fraction', ha='left', va='top', fontsize=8.5, color=ink, arrowprops=dict(arrowstyle='->', color=muted, lw=.8))
ax.set(xlim=(centre[0]-.04, centre[0]+.04), ylim=(centre[1]-.006, centre[1]+.006), xlabel='x (m)', ylabel='height (m)',
       title=f'Last iterate near vertex {w[0]} (height stretched)'); ax.grid(c=grid); ax.legend(fontsize=8, frameon=False, loc='lower right')
ax = fig.add_subplot(gs[1, 1])
ax.plot(clc, [d['pair_min_thickness_mm'] for d in data], c=muted, lw=1.4, marker='s', ms=3, label='constrained: chord-normal pair k / 272−k')
ax.plot(clc, [d['true_gap_mm'] for d in data], c=ramp[5], lw=1.8, marker='o', ms=3.5, label='true gap: nearest point on the other surface')
ax.axhline(2, c=ink, ls='--', lw=.9); ax.text(clc[0], 2.04, '2 mm floor', fontsize=8, color=ink, va='bottom')
ax.invert_xaxis(); ax.grid(c=grid); ax.set_ylim(bottom=0)
ax.set(xlabel='$C_L c$ (m), each point a saved iterate', ylabel='minimum thickness (mm)', title='What the constraint measured vs the real gap')
ax.legend(fontsize=8, frameon=False, loc='lower left')
ax = fig.add_subplot(gs[1, 2])
ax.plot(clc, [d['reversed_area_cm2'] for d in data], c=reverse_c, lw=1.8, marker='o', ms=3.5)
ax.invert_xaxis(); ax.grid(c=grid); ax.set_ylim(bottom=0)
ax.set(xlabel='$C_L c$ (m)', ylabel='area with u < 0 (cm²)', title='Recirculation zone behind the section')
save(fig, f'cp_iterates_{a.name}')

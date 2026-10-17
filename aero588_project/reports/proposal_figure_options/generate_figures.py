"""Two proposal figure options: saved CFD grid and an illustrative design schematic."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/aer588-proposal-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CASE = ROOT / 'cases/naca0012_c050_h005_u200_am3'
FIELDS = np.load(CASE / 'reports/plots/baseline_fields.npz')
SURFACE = np.genfromtxt(CASE / 'results/surface.csv', delimiter=',', skip_header=1)
NAVY, BLUE, GREY, ORANGE = '#243b4a', '#256b91', '#b2c0c9', '#bd5b27'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'svg.fonttype': 'none', 'savefig.facecolor': 'white'})

def save(fig, name):
    for ext in ['png', 'svg']:
        fig.savefig(OUT / f'{name}.{ext}', dpi=320, bbox_inches='tight', pad_inches=.055)
    plt.close(fig)

def mesh():
    fig, ax = plt.subplots(figsize=(6.5, 3.8), layout='constrained')
    for tag, stride, color, lw in [('background', 4, GREY, .38), ('near', 3, BLUE, .45)]:
        xx, yy = [FIELDS[f'{tag}_{name}'] for name in ['CoordinateX', 'CoordinateY']]
        ii = sorted(set(range(0, xx.shape[0], stride)) | {xx.shape[0] - 1})
        jj = sorted(set(range(0, xx.shape[1], stride)) | {xx.shape[1] - 1})
        segments = [np.column_stack([xx[i, :], yy[i, :]]) for i in ii]
        segments += [np.column_stack([xx[:, j], yy[:, j]]) for j in jj]
        ax.add_collection(LineCollection(segments, colors=color, linewidths=lw))
    ax.fill(SURFACE[:, 0], SURFACE[:, 1], color=NAVY, ec='#142632', lw=.9, zorder=5)
    ax.axhline(0, color=NAVY, lw=2, zorder=6)
    lowest = SURFACE[np.argmin(SURFACE[:, 1])]
    ax.annotate('', (lowest[0], lowest[1]), (lowest[0], 0),
                arrowprops={'arrowstyle': '<->', 'color': ORANGE, 'lw': 1.5}, zorder=7)
    ax.text(lowest[0]+.014, .023, '5 cm', color=ORANGE, fontsize=11,
            rotation=90, ha='left', va='center', zorder=7)
    ax.set(xlim=(-.10, .65), ylim=(-.008, .30), xlabel='$X$ (m)', ylabel='$Y$ (m)')
    ax.set_aspect('equal')
    ax.set_title('NACA 0012: existing overset mesh', loc='left', fontsize=13, pad=9)
    ax.text(.02, .97, '$c=0.50$ m, $h_{\\min}=0.05$ m', transform=ax.transAxes,
            va='top', fontsize=11, color=NAVY,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': 2})
    ax.set_xticks([-.1, 0, .1, .2, .3, .4, .5, .6])
    ax.set_yticks([0, .1, .2, .3])
    save(fig, 'option_a_actual_mesh')

def airfoil(camber, angle):
    x = (1-np.cos(np.linspace(0, np.pi, 260)))/2
    yt = .6*(.2969*np.sqrt(x)-.126*x-.3516*x*x+.2843*x**3-.1015*x**4)
    center = -camber*np.sin(np.pi*x)
    a = np.deg2rad(angle)
    def rotate(y):
        return np.cos(a)*(x-.25)-np.sin(a)*y+.25, np.sin(a)*(x-.25)+np.cos(a)*y
    xu, yu = rotate(center+yt)
    xl, yl = rotate(center-yt)
    shift = .1-min(yu.min(), yl.min())
    return xu, yu+shift, xl, yl+shift

def schematic():
    fig, ax = plt.subplots(figsize=(7.6, 3.0), layout='constrained')
    for offset, camber, angle, label in [(0, 0, 3, 'Initial section'),
                                         (1.72, .075, 10, 'Candidate shape')]:
        xu, yu, xl, yl = airfoil(camber, angle)
        ax.add_patch(Rectangle((offset-.035, .1), 1.08, .31, fill=False,
                               ec='#a1a9af', lw=1, ls=(0, (4, 3))))
        ax.text(offset+.5, .435, 'Design box', ha='center', va='bottom', fontsize=11)
        ax.fill(np.r_[xu, xl[::-1]]+offset, np.r_[yu, yl[::-1]],
                facecolor='#e4ebf7', edgecolor=NAVY, lw=1.5, zorder=4)
        # A fitted illustrative FFD cage; these are not an existing optimized design.
        sample = np.linspace(0, 1, 9)
        for xx, yy, side in [(xu, yu, 1), (xl, yl, -1)]:
            t = np.linspace(0, 1, len(xx))
            px = np.interp(sample, t, xx)+offset
            py = np.interp(sample, t, yy)+side*.028
            ax.plot(px, py, color='#a9b7cf', lw=.7, ls='--', zorder=3)
            ax.plot(px, py, 'o', mfc='white', mec='#4b70c8', ms=4.7, mew=1.05, zorder=5)
        ax.annotate('', (offset+.38, .1), (offset+.38, 0),
                    arrowprops={'arrowstyle': '<->', 'lw': 1.05, 'color': NAVY})
        ax.text(offset+.42, .04, '$h_{\\min}$', fontsize=11, color=NAVY)
        ax.text(offset+.5, .54, label, ha='center', fontsize=12, color=NAVY)
    ax.annotate('', (1.58, .26), (1.16, .26),
                arrowprops={'arrowstyle': '-|>', 'lw': 1.1, 'color': NAVY})
    ax.text(1.37, .30, 'Shape search', ha='center', va='bottom', fontsize=10)
    ax.annotate('', (-.07, .25), (-.39, .25),
                arrowprops={'arrowstyle': '-|>', 'lw': 1.3, 'color': NAVY})
    ax.text(-.24, .29, '$U_\\infty$', ha='center', fontsize=12)
    ax.plot([-.4, 2.89], [0, 0], color=NAVY, lw=1.4)
    for xx in np.arange(-.4, 2.9, .055):
        ax.plot([xx, xx-.026], [0, -.025], color='#a6afb8', lw=.65)
    ax.annotate('', (2.55, -.105), (2.21, -.105),
                arrowprops={'arrowstyle': '-|>', 'lw': 1.3, 'color': NAVY})
    ax.text(1.29, -.11, 'Moving no-slip ground: $U_{\\rm wall}=U_\\infty$',
            ha='center', va='center', fontsize=11, color=NAVY)
    ax.plot([.82], [-.205], 'o', mfc='white', mec='#4b70c8', ms=5, mew=1.1)
    ax.text(.87, -.205, 'Illustrative FFD control points', ha='left', va='center', fontsize=10)
    ax.set(xlim=(-.43, 2.95), ylim=(-.25, .62))
    ax.axis('off')
    save(fig, 'option_b_shape_schematic')

def compact_schematic():
    """Stack the two shapes so that the figure is legible in a half-page column."""
    fig, ax = plt.subplots(figsize=(3.5, 3.3), layout='constrained')
    for shift, camber, angle, label in [(0, 0, 3, 'Initial section'),
                                       (-.77, .075, 10, 'Candidate shape')]:
        xu, yu, xl, yl = airfoil(camber, angle)
        ax.add_patch(Rectangle((-.035, .1+shift), 1.08, .31, fill=False,
                               ec='#a1a9af', lw=.9, ls=(0, (4, 3))))
        ax.text(.5, .445+shift, label, ha='center', va='bottom', fontsize=10, color=NAVY)
        ax.fill(np.r_[xu, xl[::-1]], np.r_[yu, yl[::-1]]+shift,
                facecolor='#e4ebf7', edgecolor=NAVY, lw=1.2, zorder=4)
        sample = np.linspace(0, 1, 9)
        for xx, yy, side in [(xu, yu, 1), (xl, yl, -1)]:
            t = np.linspace(0, 1, len(xx))
            px = np.interp(sample, t, xx)
            py = np.interp(sample, t, yy)+side*.028+shift
            ax.plot(px, py, color='#a9b7cf', lw=.6, ls='--', zorder=3)
            ax.plot(px, py, 'o', mfc='white', mec='#4b70c8', ms=3.0, mew=.8, zorder=5)
        ax.annotate('', (1.10, .1+shift), (1.10, shift),
                    arrowprops={'arrowstyle': '<->', 'lw': .9, 'color': NAVY})
        ax.text(1.03, .04+shift, '$h_{\\min}$', ha='right', fontsize=8, color=NAVY)
        ax.plot([-.36, 1.12], [shift, shift], color=NAVY, lw=1.15)
        for xx in np.arange(-.35, 1.13, .065):
            ax.plot([xx, xx-.024], [shift, shift-.025], color='#a6afb8', lw=.5)
    ax.annotate('', (.66, -.18), (.66, -.065),
                arrowprops={'arrowstyle': '-|>', 'lw': 1.0, 'color': NAVY})
    ax.text(.60, -.12, 'Shape search', ha='right', va='center', fontsize=9, color=NAVY)
    ax.annotate('', (-.07, .25), (-.36, .25),
                arrowprops={'arrowstyle': '-|>', 'lw': 1.1, 'color': NAVY})
    ax.text(-.24, .29, '$U_\\infty$', ha='center', fontsize=9)
    ax.text(.48, .575, 'Design box + FFD', ha='center', fontsize=10, color=NAVY)
    ax.text(.30, -.90, 'Moving ground: $U_{\\rm wall}=U_\\infty$',
            ha='center', fontsize=8, color=NAVY)
    ax.annotate('', (1.10, -.875), (.87, -.875),
                arrowprops={'arrowstyle': '-|>', 'lw': 1.1, 'color': NAVY})
    ax.plot([-.16], [-1.01], 'o', mfc='white', mec='#4b70c8', ms=3.5, mew=.8)
    ax.text(-.09, -1.01, 'Illustrative FFD points', va='center', fontsize=8, color=NAVY)
    ax.set(xlim=(-.40, 1.16), ylim=(-1.07, .67))
    ax.axis('off')
    save(fig, 'option_b_compact')

if __name__ == '__main__':
    mesh()
    schematic()
    compact_schematic()
    print('Created both proposal options from the saved grid and an illustrative shape cage.')

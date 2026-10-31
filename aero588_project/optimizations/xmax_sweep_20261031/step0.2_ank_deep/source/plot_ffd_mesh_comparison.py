"""Plot the saved CFD meshes and verified pyGeo FFD coefficients; no CFD solve."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    root = parser.parse_args().root
    out = root / 'plots'
    out.mkdir(exist_ok=True)
    lattices = np.load(root / 'ffd_geometry/control_points.npz')
    labels = ['baseline', 'final']
    nodes = [lattices[label][:, :, 0, :2] for label in labels]
    nstation = nodes[0].shape[0]
    surfaces = [np.loadtxt(root / label / 'surface_geometry.csv', delimiter=',', skiprows=1)
                for label in labels]
    results = [json.loads((root / label / 'result.json').read_text()) for label in labels]
    colors = ['#175676', '#ce5927']
    titles = ['Before: initial NACA 0012', 'After: design at iteration 50']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'font.family': 'DejaVu Sans',
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white'})

    def save(fig, name):
        for suffix in ['png', 'svg']:
            fig.savefig(out / (name + '.' + suffix), dpi=220, bbox_inches='tight')
        plt.close(fig)

    def lattice(ax, points, color, style='-', alpha=1.):
        for j in range(2):
            ax.plot(points[:, j, 0], points[:, j, 1], style, color=color,
                    lw=1.1, alpha=alpha, zorder=10)
        ax.add_collection(LineCollection(points, colors=color, linewidths=.8,
                                         linestyles=style, alpha=alpha, zorder=10))
        ax.scatter(points[..., 0], points[..., 1], s=40, facecolors='white',
                   edgecolors=color, linewidths=1.5, zorder=12)

    def axes_setup(ax):
        ax.axhline(0, color='#253944', lw=1.6)
        ax.text(.545, .003, 'Slip ground', ha='right', va='bottom', fontsize=8,
                color='#253944', bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1})
        ax.set(xlim=(-.045, .55), ylim=(-.009, .18), xlabel='x (m)',
               ylabel='Height above ground (m)')
        ax.set_aspect('equal')

    handles = [Line2D([], [], color='#40778d', lw=.6, label='Near-body mesh'),
               Line2D([], [], color='#b8c6cd', lw=.6, label='Background mesh'),
               Line2D([], [], color='#7040a0', marker='o', markerfacecolor='white',
                      lw=1.1, label='FFD lattice and control points')]
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.2))
    for k, (ax, label) in enumerate(zip(axes, labels)):
        fields = np.load(root / label / 'fields.npz')
        for tag, stride, color, width in [('background', 3, '#b8c6cd', .35),
                                           ('near', 3, '#40778d', .4)]:
            x, y = [fields[tag + '_' + dim] for dim in ['CoordinateX', 'CoordinateY']]
            segments = [np.c_[x[i, :], y[i, :]] for i in range(0, x.shape[0], stride)]
            segments += [np.c_[x[:, j], y[:, j]] for j in range(0, x.shape[1], stride)]
            ax.add_collection(LineCollection(segments, colors=color, linewidths=width,
                                              zorder=1, rasterized=False))
        surface = surfaces[k]
        ax.fill(surface[:, 0], surface[:, 1], color='#263947', zorder=5)
        lattice(ax, nodes[k], '#7040a0')
        axes_setup(ax)
        ax.set_title(titles[k] + '\n' +
                     f"$C_L$ = {results[k]['cl']:.4f}   |   $C_D$ = {results[k]['cd']:.5f}",
                     fontsize=11, pad=10)
    fig.suptitle('FFD control points over the actual computational mesh', fontsize=15, y=.99)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .91),
               ncol=3, frameon=False)
    fig.subplots_adjust(left=.055, right=.99, bottom=.19, top=.76, wspace=.16)
    fig.text(.5, .055, f'{4*nstation} control points ({nstation} × 2 × 2); the two tied span planes overlap into {2*nstation} visible pairs.\n'
             'Selected mesh lines shown. The final design reached the iteration limit before optimizer convergence.',
             ha='center', fontsize=9, color='#4b5358')
    save(fig, 'ffd_mesh_comparison')

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.9),
                             gridspec_kw={'height_ratios': [2., 1.]})
    ax = axes[0]
    for k in range(2):
        surface = surfaces[k]
        ax.fill(surface[:, 0], surface[:, 1], color=colors[k], alpha=.09, zorder=2)
        ax.plot(surface[:, 0], surface[:, 1], color=colors[k], lw=1.8, zorder=6)
        lattice(ax, nodes[k], colors[k], '--' if k == 0 else '-', alpha=.85)
    start, end = nodes[0].reshape(-1, 2), nodes[1].reshape(-1, 2)
    ax.quiver(start[:, 0], start[:, 1], *(end - start).T, angles='xy',
              scale_units='xy', scale=1, color='#383e43', width=.0019,
              headwidth=3.5, headlength=4.5, zorder=13)
    axes_setup(ax)
    ax.set_title('Before / after: airfoil and FFD control-point movement', fontsize=14, pad=37)
    legend = [Line2D([], [], color=colors[k], lw=1.6, marker='o',
                     markerfacecolor='white', ls='--' if k == 0 else '-', label=titles[k])
              for k in range(2)]
    ax.legend(handles=legend, loc='lower center', bbox_to_anchor=(.5, 1.005),
              ncol=2, frameon=False)
    displacement = (nodes[1] - nodes[0]) * 1000
    for j, marker, color, label in [(0, 'o', '#175676', 'Lower FFD row'),
                                     (1, 's', '#ce5927', 'Upper FFD row')]:
        axes[1].plot(np.arange(1, nstation+1), displacement[:, j, 1], marker=marker,
                      color=color, lw=1.6, ms=5, label=label)
    axes[1].axhline(0, color='#90989e', lw=.8)
    axes[1].set(xlabel='FFD chordwise station (leading edge → trailing edge)',
                ylabel='Control-point Δy (mm)', xticks=np.arange(1, nstation+1), xlim=(.75, nstation+.25))
    axes[1].grid(alpha=.2)
    axes[1].legend(loc='lower left', frameon=False, ncol=2)
    fig.subplots_adjust(left=.10, right=.98, top=.87, bottom=.15, hspace=.45)
    fig.text(.5, .035, 'Geometry and arrows shown at true scale; Δy includes shape and incidence changes.\n'
             'Final design: 50 PSQP iterations; optimizer convergence was not reached.',
             ha='center', fontsize=9, color='#4b5358')
    save(fig, 'ffd_before_after_comparison')
    print('Saved FFD mesh and before/after comparisons (PNG and SVG).')


if __name__ == '__main__':
    main()

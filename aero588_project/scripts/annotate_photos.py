"""Arrows on the two photos of the report (optimizations/comparisons/report_20261122/photos_src, supplied by the user):
  photo_front_wing.jpg   F1 car, arrow at the front wing, mid section
  photo_gurney.jpg       front wing flap, arrow at the Gurney flap on its trailing edge (above the top of the "h")
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from PIL import Image

C = Path(__file__).resolve().parents[1]/'optimizations/comparisons/report_20261122'
ORANGE = '#eb6834'
halo = [pe.withStroke(linewidth=4, foreground='black')]

def annotate(src, out, head, tail, label, ha):
    im = Image.open(C/'photos_src'/src).convert('RGB'); w, h = im.size
    fig = plt.figure(figsize=(w/100, h/100), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(im); ax.axis('off')
    ax.annotate('', head, tail, arrowprops=dict(arrowstyle='-|>,head_length=1.1,head_width=.55', color=ORANGE, lw=4.5,
                                                 shrinkA=0, shrinkB=2, path_effects=[pe.withStroke(linewidth=7, foreground='black')]))
    ax.text(*tail, label, color='white', fontsize=20, fontweight='bold', ha=ha, va='center', path_effects=halo)
    fig.savefig(C/out, dpi=100, pil_kwargs={'quality': 90}); plt.close(fig); print(C/out)

annotate('f1_car.png', 'photo_front_wing.jpg', head=(800, 626), tail=(1010, 770), label='front wing, mid section  ', ha='right')
annotate('front_wing_gurney.png', 'photo_gurney.jpg', head=(400, 211), tail=(215, 95), label='Gurney flap  ', ha='right')

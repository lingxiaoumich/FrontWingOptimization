"""Build reports/report_20261127.html: one self-contained page (figures and GIFs embedded as base64), sections behind
buttons. Numbers come from runs.json (report_gather.py). Figures: plot_report.py for sections 1-7 (unchanged since
report_20261122, in comparisons/report_20261122/), plot_report_20261127.py for sections 8-9 (comparisons/report_20261127/),
and the runs' chain GIFs (plot_run_chain.py). From the project folder, in the MACH-Aero env:
  python scripts/report_gather.py optimizations/comparisons/report_20261127/runs.json
  python scripts/plot_report_20261127.py optimizations/comparisons/report_20261127/runs.json
  python scripts/make_report.py optimizations/comparisons/report_20261127/runs.json --no-photos   # reports/report_20261127.html (in git)
  python scripts/make_report.py optimizations/comparisons/report_20261127/runs.json --out reports/report_20261127_photos.html
The last copy is local only: the two photos (comparisons/report_20261122/photo_*.jpg) are supplied, third-party.
The previous version, reports/report_20261122.html, was built by this script as of commit a9171d15.
"""
import argparse
import base64
import json
import sys
from pathlib import Path

P = Path(__file__).resolve().parents[1]
O, CMP = P/'optimizations', P/'optimizations/comparisons/report_20261122'
CMP2 = P/'optimizations/comparisons/report_20261127'
ap = argparse.ArgumentParser(); ap.add_argument('runs', type=Path); ap.add_argument('--out', type=Path, default=P/'reports/report_20261127.html')
ap.add_argument('--no-photos', action='store_true'); A = ap.parse_args()
OUT = A.out.resolve()
R = json.loads(A.runs.read_text())
A0 = 0.021549464275835817
Q = 0.5*1.2249912464200752*55.55555555555556**2  # Pa

def b64(path):
    return base64.b64encode(Path(path).read_bytes()).decode()

def fig(path, caption, width=100):
    kind = {'.gif': 'gif', '.jpg': 'jpeg', '.jpeg': 'jpeg'}.get(Path(path).suffix.lower(), 'png')
    return (f'<figure><img style="width:{width}%" src="data:image/{kind};base64,{b64(path)}" alt="">'
            f'<figcaption>{caption}</figcaption></figure>')

def table(head, rows, cls=''):
    h = ''.join(f'<th>{c}</th>' for c in head)
    body = ''.join('<tr>'+''.join(f'<td>{c}</td>' for c in r)+'</tr>' for r in rows)
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{body}</tbody></table>'

def cpmin_surface(run):
    b = R[run]['best']; d = O/run/b['round']/'best_feasible'/'flow_000_slices.dat'
    lines = d.read_text().split('\n'); n = int(lines[3].split()[2])
    return min(float(l.split()[6]) for l in lines[5:5+n])

def B(run): return R[run]['best']
def dn(run): return f"{-B(run)['clc']:.2f}"
def dr(run): return f"{B(run)['cdc']:.3f}"
def ls(run): return '–' if B(run)['loss'] is None else f"{100*B(run)['loss']:.2f}"
def low(run): return f"{B(run)['clearance_mm']:.0f}"
def area(run): return f"{B(run)['area_pct']:.0f}"
def pct(a, b): return f"{100*(a/b-1):+.0f} %"
def meshes(run): return str(len(R[run]['rounds']))
def hours(run): return f"{R[run]['wall_h']:.1f}"

R16 = 'run_016_area050'
s16 = R[R16]['start']
# sections 3-7 keep the GIFs of report_20261122 (four panels, as their captions say), copied out of git (commit a9171d15)
# into comparisons/report_20261127/gifs/; the chain GIFs were redrawn on 11-25 with a constraint column, twice the size
gifs = {k: CMP2/'gifs'/f'{k}.gif' for k in
        ['run_016_area050', 'run_029_area050DragCap002', 'run_021_area050CptLoss025', 'run_024_area050Floor12p5mm',
         'run_035_area050FreeAir', 'run_027_area050CpMin4']}
gifs['run_038_area050MinLossClc14CpMin4CapStart'] = O/'run_038_area050MinLossClc14CpMin4CapStart/summary/plots/chain_iterates.gif'

KKTR = [('run_037_area050MinLossClc08CpMin4CapStart', '0.8'), ('run_036_area050MinLossClc10CpMin4CapStart', '1.0'),
        ('run_039_area050MinLossClc12CpMin4', '1.2'), ('run_038_area050MinLossClc14CpMin4CapStart', '1.4')]
R38, R39 = KKTR[3][0], KKTR[2][0]
R40, R41, R42, R27 = 'run_040_area050RotStart', 'run_041_area050DragCap002RotStart', 'run_042_area050MinSuctionClc12CapStart', 'run_027_area050CpMin4'
R29 = 'run_029_area050DragCap002'
PUSH = dict(evals=sum(R[r]['kkt_evals'] for r, _ in KKTR), h=sum(R[r]['kkt_wall_h'] for r, _ in KKTR), core_h=sum(R[r]['kkt_core_h'] for r, _ in KKTR))
ALL = dict(evals=sum(d.get('evals', 0) for d in R.values()), h=sum(d.get('wall_h', 0) for d in R.values()),
           core_h=sum(d.get('core_h', d.get('wall_h', 0)*(d.get('ranks') or 24)) for d in R.values()))
NEW = dict(evals=sum(R[r]['evals'] for r in (R40, R41, R42)), h=sum(R[r]['wall_h'] for r in (R40, R41, R42)))
LABEL = {'lift': 'downforce', 'section_area': 'area', 'cp_min': 'C<sub>p</sub>', 'local_thickness': 'local thickness', 'te_base': 'TE base',
         'surface_y': 'box top', 'surface_x': 'box ends', 'te_taper': 'TE taper', 'max_thickness': 'max thickness'}
def active(b):
    n = {}
    for name, lam in b['multipliers']:
        if lam > 1e-3: g = name.split('[')[0]; n[g] = n.get(g, 0)+1
    return ', '.join(LABEL.get(g, g)+('' if k == 1 else f' &times;{k}') for g, k in n.items())
def e(v):  # 1.2e-04 -> 1.2&times;10<sup>&minus;4</sup>
    m, x = f'{v:.1e}'.split('e'); return f'{m}&times;10<sup>&minus;{int(x[1:])}</sup>' if x.startswith('-') else f'{m}e{x}'

S = {}
S['update'] = ('What is new', f'''
<p>This version (27 November 2026) adds to the report of 22 November:</p>
<ul>
<li><b>Runs 036&ndash;039 converged</b> (section 8). They minimize the wake loss with a downforce limit and
C<sub>p</sub> &ge; &minus;4. On 22 November run_038 had no feasible design and run_039 had stopped just short of its
downforce limit. All four were run on until the optimizer&rsquo;s own first-order optimality (KKT) test was met.</li>
<li><b>The loss measurement was corrected</b> (section 8). Getting there meant fixing how the loss is measured and how
the far-field boundary moves with the mesh. Loss values up to run_035 (sections 3&ndash;7) use the old measurement,
which reads 9&ndash;12 % low; values marked &ldquo;corrected&rdquo; use the new one. The two should not be compared.</li>
<li><b>Runs 040&ndash;042</b> (section 9): a start section already tilted toward the top corner of the box (040, 041),
and a run that minimizes the suction peak (042).</li>
<li>Methods, the C<sub>p</sub>-limit section and the conclusion are updated to match. Sections 3&ndash;6 are unchanged
apart from a note on the loss values.</li>
</ul>
''')

S['background'] = ('1 Background', f'''
<p>A race car front wing works close to the ground. Its main job is downforce: it pushes the front tyres into the road
for grip. It also sets up the flow for everything behind it (wheels, floor, radiators, rear wing), so a good front wing
makes downforce without leaving a big loss of total pressure in its wake [1, 2].</p>
{'' if A.no_photos else fig(CMP/'photo_front_wing.jpg', 'Cadillac F1 car. This study optimizes a 2-D section of the front wing, like the one at the arrow.', 75)}
<p>The rules say where the wing may be. Here that is a box in a 2-D section: x from &minus;12.3 to 512 mm and y from
50 to 188.1 mm above the ground. The section must stay inside it.</p>
<p>Near the ground the wing works in ground effect. The flow between the suction (lower) side and the ground speeds up,
so the pressure there drops much further than in free air. In our baseline the surface pressure reaches
C<sub>p</sub> &asymp; &minus;10 near the leading edge, against &minus;1.2 for the starting section. Closer to the
ground gives more downforce, until the suction side starts to separate and the downforce falls [3].</p>
{fig(CMP/'cp_ground_vs_free.png', 'Pressure around the best designs of the baseline run (with ground) and of the same setup in free air (no ground). Different designs, same colour scale.')}
<p>So the usual goal is the most downforce with little drag or wake loss, inside the box. This study uses RANS CFD with
adjoint gradients to optimize the section shape, then sweeps the drag limit, the loss limit, the ground height and a
limit on the lowest C<sub>p</sub>.</p>
''')

# ------------------------------------------------------------------ 2 methods
xdsm = fig(CMP/'xdsm.png', 'One optimizer iteration, in the XDSM convention [6, 7].', 90)

cells = {'run_022_area050Floor37p5mm': 90138, 'run_023_area050Floor25mm': 92364, 'run_024_area050Floor12p5mm': 95226}
S['methods'] = ('2 Methods', f'''
<h3>Geometry and design variables</h3>
<p>The start is a NACA 0012 section, level, stretched to the box length (524.3 mm), with its lowest point on the box
floor (50 mm). Free-form deformation (FFD) moves it [12, 13]: 20 pairs of control points, one above and one below the
section, in a flat box around it. Each of the 40 points moves only up or down; the tilt of the section is locked, so
any tilt comes from these moves. A surface point moves with a blend of the control points around it.</p>
{fig(O/R16/'summary/plots/setup_geometry.png', 'Start section, FFD control points and the box (run_016).', 85)}
<h3>Flow solver</h3>
<p>ADflow [8], steady RANS with the Spalart&ndash;Allmaras turbulence model [15] (first-order turbulence advection),
central differencing with scalar dissipation. Each solve starts with the approximate Newton&ndash;Krylov (ANK) solver
[9] and finishes with Newton&ndash;Krylov, to a residual 10<sup>&minus;10</sup> times that of free stream (3&times;10<sup>&minus;9</sup>
in the free-air run, where solves stall near 10<sup>&minus;9</sup>; 10<sup>&minus;12</sup> in runs 036&ndash;039 from 25
November, and 10<sup>&minus;13</sup> for the last steps of run_038, section 8). Flow: 55.6 m/s (200 km/h), sea-level air
(288.15 K, 101,325 Pa), Mach 0.163, Reynolds number 1.9 million on 0.5 m. The 2-D section is a one-cell-thick slice
(0.05 m) with symmetry planes on both sides.</p>
<h3>Mesh</h3>
<p>Overset grid: an O-grid grown from the section with pyHyp [11] (289 &times; 81 nodes, 25 mm thick, first cell 2
&micro;m, so y<sup>+</sup> &lt; 0.6 on the wall), over a Cartesian background (319 &times; 207 nodes, x from &minus;5 to
7.5 m, y from 0 to 5 m, 2.5 mm &times; 5 mm cells near the wing). 88,548 cells in total. The floor runs refine the background near
the ground ({cells['run_022_area050Floor37p5mm']:,} to {cells['run_024_area050Floor12p5mm']:,} cells); the free-air
run carries the domain down to y = &minus;5 m (128,616 cells). IDWarp [11] moves the mesh with the design. When the
mesh wears out (bad cells near the trailing edge), a new O-grid is grown around the latest design and the optimizer
restarts from it: mesh 0, 1, 2 and so on. In the runs up to 035 a new mesh moves C<sub>L</sub>c by about 2 %, so
smaller differences between those runs are noise. Much of that shift probably came from the far-field boundary moving
with the mesh, which the last meshes of runs 036&ndash;039 no longer do (section 8); there a new mesh at the same design
moves C<sub>L</sub>c by 0.04&ndash;0.6 %.</p>
{fig(O/R16/'summary/plots/mesh.png', 'Mesh of the start section (run_016, mesh 0).')}
<h3>Boundary conditions</h3>
{table(['boundary', 'condition'], [['wing', 'no-slip wall'], ['ground (y = 0)', 'slip wall (symmetry); stands in for a moving ground, no road boundary layer'],
      ['outer edges', 'far field at free-stream conditions: 10 chords ahead and above, 14 chords behind (x = 7.5 m). In the last meshes of 036&ndash;039 the boundary is held in place when the mesh deforms, and its top is tilted down 0.2 % toward the outlet (section 8)'],
      ['span', 'symmetry planes'], ['free-air run (035)', 'ground replaced by far field 10 chords below']])}
<h3>Optimizer</h3>
<p>PSQP, a sequential quadratic programming method with quasi-Newton updates, through pyOptSparse [14]. Gradients:
ADflow&rsquo;s discrete adjoint for the flow functions [10], pyGeo&rsquo;s analytic derivatives for the geometric ones.
A derivative check at the start of every mesh compares the adjoint against finite differences. The diagram follows
the XDSM convention [7] as used in Martins &amp; Ning [6].</p>
<p>A chain driver ends a mesh when PSQP exits or when the design stops improving for 15 minutes (the stall rule), then
remeshes, and ends the run when a new mesh no longer helps. So the best designs of runs up to 035 and of 040&ndash;042
are close to an optimum, not proven optimal. Runs 036&ndash;039 were run without the stall rule until PSQP&rsquo;s
optimality test was met (section 8).</p>
{xdsm}
<h3>Constraints</h3>
{table(['constraint', 'value', 'runs'], [
    ['every surface point inside the box', 'x &minus;12.3&hellip;512 mm, y 50&hellip;188.1 mm', 'all (floor 37.5 / 25 / 12.5 mm in 022&ndash;024)'],
    ['local thickness', '&ge; 2 mm (110 point pairs)', 'all'],
    ['maximum thickness', '&ge; 30 mm', 'all'],
    ['trailing-edge base thickness', '&ge; 1.2 mm', 'all'],
    ['trailing-edge taper', '&ge; 8&deg; over the last 14 point pairs', 'all'],
    ['tilt', 'locked', 'all'],
    ['section area', '&ge; 50 % of the start (107.75 cm&sup2;)', '016&ndash;039'],
    ['drag C<sub>D</sub>c', '&le; limit', '017, 018, 028, 029'],
    ['outlet C<sub>pt</sub> loss', '&le; limit', '019&ndash;021, 030'],
    ['lowest C<sub>p</sub> on the wing (KS aggregate)', '&ge; limit', '025&ndash;027, 036&ndash;039'],
    ['downforce &minus;C<sub>L</sub>c (when loss or suction is the objective)', '&ge; limit', '031&ndash;039, 042']])}
<p>About 830 constraints in total (289 box points in each of x and y, 220 local thickness, 28 taper, and the rest).
A design is feasible when it meets all of them to 10<sup>&minus;6</sup>; &ldquo;best feasible&rdquo; below means the best
such design of a run.</p>
<h3>Computation effort</h3>
{table(['item', 'time'], [
    ['flow solve, warm start (median)', '16&ndash;17 s on 24 cores'],
    ['first flow solve on a new mesh', '0.5&ndash;10 min'],
    ['gradient: adjoints of 2&ndash;3 functions + geometry (median)', '59 s'],
    ['PSQP&rsquo;s own work per iteration', '&asymp; 0.05 s (same-size test problem)'],
    ['one optimizer iteration', '&asymp; 75 s, plus line-search solves (1.9 flow solves per gradient)'],
    ['one run (2&ndash;6 meshes)', '1&ndash;7 h'],
    ['runs 016&ndash;039, up to 22 November', '3,068 flow evaluations, 1,653 gradients, 78.5 h, &asymp; 1,900 core-hours'],
    ['KKT push of 036&ndash;039, 25&ndash;27 November', f"{PUSH['evals']:,} flow evaluations, {PUSH['h']:.0f} h of run time (four runs at once, 12&ndash;24 cores each), &asymp; {round(PUSH['core_h'], -2):,.0f} core-hours"],
    ['runs 040&ndash;042', f"{NEW['evals']} flow evaluations, {NEW['h']:.0f} h"],
    ['all runs 001&ndash;042', f"{ALL['evals']:,} flow evaluations, {ALL['h']:.0f} h, &asymp; {round(ALL['core_h'], -2):,.0f} core-hours"]])}
<p>Flow solves take two thirds of the time, adjoints most of the rest; the optimizer itself is negligible. Cores: AMD
EPYC 7601, two runs at once.</p>
<h3>Metrics</h3>
{table(['quantity', 'definition'], [
    ['downforce &minus;C<sub>L</sub>c', 'L&prime;/q<sub>&infin;</sub> per unit span (m), sign flipped so downforce is positive. Using C<sub>L</sub>c avoids a reference chord, since the chord changes. 1.0 is 1.89 kN per metre of span.'],
    ['drag C<sub>D</sub>c', 'D&prime;/q<sub>&infin;</sub> per unit span (m)'],
    ['pressure coefficient C<sub>p</sub>', '(p &minus; p<sub>&infin;</sub>)/q<sub>&infin;</sub>'],
    ['total-pressure coefficient C<sub>pt</sub>', '(p<sub>0</sub> &minus; p<sub>&infin;</sub>)/q<sub>&infin;</sub>: 1 in clean free stream, lower in wakes'],
    ['outlet C<sub>pt</sub> loss', 'mass-averaged (p<sub>0&infin;</sub> &minus; p<sub>0</sub>)/q<sub>&infin;</sub> on a vertical plane at x = 7.3 m (14 chords behind), full domain height; shown in % of q<sub>&infin;</sub>. <b>Corrected</b> (runs 036&ndash;039 from 26 November): measured against the mass-averaged p<sub>0</sub> on a plane at x = &minus;4 m of the same solution, with the plane points at cell midpoints. The old values read 9&ndash;12 % low (section 8)'],
    ['section area', '% of the start section (215.5 cm&sup2;)'],
    ['q<sub>&infin;</sub>', f'&frac12;&rho;V&sup2; = {Q:,.0f} Pa']])}
''')

# ------------------------------------------------------------------ 3 baseline
b = B(R16)
S['baseline'] = ('3 Baseline (016)', f'''
<p>run_016: maximize downforce with all geometric limits and the 50 % area limit, no drag or loss limit.</p>
{table(['', '&minus;C<sub>L</sub>c', 'C<sub>D</sub>c', 'outlet loss (%)', 'area (%)', 'lowest C<sub>p</sub>', 'lowest point (mm)'], [
    ['start: level NACA 0012', f"{-s16['clc']:.2f}", f"{s16['cdc']:.4f}", f"{100*s16['loss']:.2f}", '100', '&minus;1.2', '50'],
    ['best feasible', dn(R16), dr(R16), ls(R16), area(R16), f"{cpmin_surface(R16):.1f}", low(R16)]])}
<p>Downforce goes up almost ten times. The best design came after {R[R16]['evals']} flow evaluations on 4 meshes,
{hours(R16)} h on 24 cores. Most of the gain comes in the first two meshes; the last two add about 1 %.</p>
{fig(O/R16/'summary/plots/result_overview.png', '(a) Start and best sections. (b) Surface C<sub>p</sub>: the suction peak reaches &minus;9.9 near the leading edge. (c) Downforce of each optimizer iteration.')}
<p>What the optimizer did:</p>
<ul>
<li>The nose drops onto the box floor and the suction side runs close to the ground over the front half, where the
pressure is lowest.</li>
<li>The section thins to the area limit.</li>
<li>The trailing edge turns up sharply to the top corner of the box: the last 6 % of the chord is steeper than 45&deg;
and rises 64 mm (13 % of the chord).</li>
</ul>
{fig(O/R16/'summary/plots/trailing_edge.png', 'Trailing edge of the best design. The steep part of the pressure side (orange) works like a Gurney flap.', 45)}
<p><b>Gurney flap.</b> A Gurney flap is a short tab at right angles to the chord, at the trailing edge on the pressure
side, usually a few percent of the chord high [5]. It turns the flow at the trailing edge and adds lift, for some drag.
Zerihan and Zhang measured the same gain in downforce on a wing in ground effect [4]. The optimizer found this shape on
its own, from a plain NACA section: a steep turn on the pressure side, right at the trailing edge. Here it is taller
and blends into the section, closer to a flap, and the top of the box stops it.</p>
{'' if A.no_photos else fig(CMP/'photo_gurney.jpg', 'Ferrari front wing flap: a Gurney flap, a small lip along the trailing edge (arrow).', 60)}
<p><b>Top surface.</b> The top (pressure) side changes a lot from run to run (bumps, thick or thin) with little effect
on downforce; from experience its sensitivity is low. The suction side and the trailing-edge turn set the downforce.
None of the runs limits surface curvature, and some surfaces come out wavy.</p>
{fig(gifs[R16], 'run_016, every optimizer iteration. (a) C<sub>pt</sub>, (b) C<sub>p</sub>, (c) section and geometric limits, (d) downforce against drag. In the tables, bold: at a limit (within 1 %); underlined: past it (by more than 0.1 %).')}
''')

# ------------------------------------------------------------------ 4 drag sweep
DRAG = [(R16, 'none'), ('run_018_area050DragCap008', '0.08'), ('run_028_area050DragCap005', '0.05'),
        ('run_017_area050DragCap003', '0.03'), ('run_029_area050DragCap002', '0.02')]
base = -B(R16)['clc']
S['drag'] = ('4 Drag limit', f'''
<p>run_016 with a limit on C<sub>D</sub>c. Everything else is the same.</p>
{table(['run', 'C<sub>D</sub>c limit', '&minus;C<sub>L</sub>c', 'vs 016', 'C<sub>D</sub>c', 'outlet loss (%)', 'lowest point (mm)', 'area (%)', 'meshes'],
       [[r[4:7], lim, dn(r), '&ndash;' if r == R16 else pct(-B(r)['clc'], base), dr(r), ls(r), low(r), area(r), meshes(r)] for r, lim in DRAG])}
{fig(CMP/'drag_sweep.png', 'Best feasible design of each run. The drag limit is active in every run.', 62)}
<ul>
<li>Halving drag (0.16 &rarr; 0.08) costs 11 % of the downforce; a limit of 0.02 (one eighth) costs 36 %. The first
cuts in drag are cheap.</li>
<li>The best designs do not always sit on the box floor. Under a drag limit the wing moves up, to 63&ndash;106 mm, and
gives up some ground effect to cut drag.</li>
<li>The section stays at the 50 % area limit in every run.</li>
</ul>
{fig(CMP/'contact_drag.png', 'Best feasible sections, on the same grid in every panel.')}
{fig(gifs['run_029_area050DragCap002'], 'run_029 (C<sub>D</sub>c &le; 0.02): the wing lifts off the floor as the optimizer trades downforce for drag.')}
''')

# ------------------------------------------------------------------ 5 loss sweep
LOSS = [(R16, 'none'), ('run_019_area050CptLoss075', '75 %'), ('run_020_area050CptLoss050', '50 %'),
        ('run_021_area050CptLoss025', '25 %'), ('run_030_area050CptLoss015', '15 %')]
ML = [('run_031_area050MinLossClc18CapStart', '1.8', '016-type, area &ge; 50 %'), ('run_032_area050MinLossClc16CapStart', '1.6', ''),
      ('run_033_area050MinLossClc14CapStart', '1.4', ''), ('run_034_area050MinLossClc12', '1.2', 'level start')]
def mlrow(r, t):
    d = R[r]; b_ = d.get('best'); run_ = ' (still running)' if d.get('running') else ''
    start = 'level' if 'CapStart' not in r else 'loss-cap design'
    if b_: return [r[4:7], t, start, dn(r), ls(r), dr(r), 'yes'+run_]
    return [r[4:7], t, start, '&ndash;', '&ndash;', '&ndash;', 'no (stopped 0.1 % short of the downforce limit)'+run_]
S['loss'] = ('5 Loss limit', f'''
<p>run_016 with a limit on the outlet C<sub>pt</sub> loss, set as a fraction of run_016&rsquo;s best loss
({100*B(R16)['loss']:.2f} %).</p>
<p><i>Loss values in this section use the old measurement, which reads 9&ndash;12 % low and jumps between nearby
designs (section 8). Comparisons within the section hold roughly; the values should not be compared with the corrected
ones of section 8.</i></p>
{table(['run', 'loss limit', '&minus;C<sub>L</sub>c', 'vs 016', 'C<sub>D</sub>c', 'outlet loss (%)', 'lowest point (mm)', 'area (%)', 'meshes'],
       [[r[4:7], lim if r == R16 else f"{lim} ({100*R[r]['loss_cap']:.2f} %)", dn(r), '&ndash;' if r == R16 else pct(-B(r)['clc'], base), dr(r), ls(r), low(r), area(r), meshes(r)] for r, lim in LOSS])}
{fig(CMP/'loss_sweep.png', '(a) Downforce against drag; grey: the drag sweep. (b) Downforce against outlet loss; squares: runs that minimize loss with a downforce limit (036 and 037 as of 22 November, old measurement; see section 8).')}
<ul>
<li>The loss limit is active in every run. Cutting the loss by a quarter costs 5 % of the downforce; by 85 %,
41 %.</li>
<li>Capping loss gives the same downforce-for-drag trade as capping drag (plot a): the two fronts lie within the 2 %
mesh noise.</li>
<li>Tighter loss limits lower the trailing edge (188 &rarr; 105 mm) and flatten the section.</li>
</ul>
<p><b>Loss as the objective.</b> Runs 031&ndash;034 swap the roles: minimize loss with a downforce limit. Started from
the level section, they ended worse than the loss-limit runs: run_034 reached 0.40 % loss at downforce 1.20, while
run_030 passed a feasible design with 0.30 % at 1.20 on its way up. Started from a design of the loss-limit runs
(&ldquo;loss-cap start&rdquo;), they end on the same front and improve their start by 5&ndash;6 %.</p>
{table(['run', 'downforce limit', 'start', '&minus;C<sub>L</sub>c', 'outlet loss (%)', 'C<sub>D</sub>c', 'feasible design'],
       [mlrow(r, t) for r, t, _ in ML])}
{fig(CMP/'contact_loss.png', 'Best feasible sections, on the same grid in every panel.')}
{fig(gifs['run_021_area050CptLoss025'], 'run_021 (loss &le; 25 % of run_016&rsquo;s).')}
''')

# ------------------------------------------------------------------ 6 ground height
FL = [(R16, '50'), ('run_022_area050Floor37p5mm', '37.5'), ('run_023_area050Floor25mm', '25'), ('run_024_area050Floor12p5mm', '12.5')]
FR = 'run_035_area050FreeAir'
S['floor'] = ('6 Ground height', f'''
<p>run_016 with the box floor lowered to 37.5, 25 and 12.5 mm (the top stays at 188.1 mm), and once with no ground at
all (free air, box unchanged).</p>
{table(['run', 'box floor (mm)', '&minus;C<sub>L</sub>c', 'vs 016', 'C<sub>D</sub>c', 'vs 016', 'lowest point (mm)', 'meshes', 'hours'],
       [[r[4:7], f, dn(r), '&ndash;' if r == R16 else pct(-B(r)['clc'], base), dr(r), '&ndash;' if r == R16 else pct(B(r)['cdc'], B(R16)['cdc']), low(r), meshes(r), hours(r)] for r, f in FL]
       + [['035', 'free air', dn(FR), pct(-B(FR)['clc'], base), dr(FR), pct(B(FR)['cdc'], B(R16)['cdc']), low(FR)+' (no ground)', meshes(FR), hours(FR)]])}
{fig(CMP/'floor_sweep.png', '(a) Downforce against drag. (b) Downforce against the height of the lowest point of the best design (h/c with c = 0.5 m); dashed: free air.')}
<ul>
<li>Lower is better only down to about 25 mm: +5 % downforce at 25 mm, for +53 % drag.</li>
<li>With 12.5 mm allowed, the best design sits at 23 mm and makes less downforce than at 25 mm. Zerihan and Zhang saw
the same in a wind tunnel [3]: downforce rises as the wing nears the ground, peaks, then drops at small heights as
the flow on the suction side separates.</li>
<li>Without the ground the same setup reaches only {dn(FR)}, {pct(-B(FR)['clc'], base)} against run_016, at a quarter
of the drag. The outlet loss is not meaningful in the free-air domain and is not used there.</li>
</ul>
{fig(CMP/'contact_floor.png', 'Best feasible sections, on the same grid in every panel.')}
{fig(gifs['run_024_area050Floor12p5mm'], 'run_024 (floor 12.5 mm).')}
{fig(gifs[FR], 'run_035 (free air).')}
''')

# ------------------------------------------------------------------ 7 Cp min
CP = [(R16, 'none'), ('run_025_area050CpMin8', '&minus;8'), ('run_026_area050CpMin6', '&minus;6'), ('run_027_area050CpMin4', '&minus;4')]
CPML = [('run_036_area050MinLossClc10CpMin4CapStart', '1.0'), ('run_037_area050MinLossClc08CpMin4CapStart', '0.8'),
        ('run_038_area050MinLossClc14CpMin4CapStart', '1.4'), ('run_039_area050MinLossClc12CpMin4', '1.2')]
S['cpmin'] = ('7 C<sub>p</sub> limit', f'''
<p>run_016 with a lower limit on C<sub>p</sub> anywhere on the wing. A very low suction peak is followed by a steep
pressure rise, which risks separation; the limit keeps the peak mild. In the saved flows the lowest C<sub>p</sub> of
the whole field sits on the wing, so a limit on the wing covers the field.</p>
{table(['run', 'C<sub>p</sub> limit', '&minus;C<sub>L</sub>c', 'vs 016', 'C<sub>D</sub>c', 'outlet loss (%)', 'area (%)', 'lowest C<sub>p</sub> on the surface'],
       [[r[4:7], lim, dn(r), '&ndash;' if r == R16 else pct(-B(r)['clc'], base), dr(r), ls(r), area(r), f"{cpmin_surface(r):.1f}"] for r, lim in CP])}
{fig(CMP/'cpmin_sweep.png', '(a) Downforce against drag. (b) Suction-side C<sub>p</sub> of the best designs: the limit cuts the peak into a flat top. (c) Downforce against section area.')}
<ul>
<li>The limit flattens the suction peak near the leading edge into a plateau at the limit.</li>
<li>Downforce drops fast: &minus;5 % at &minus;8, &minus;18 % at &minus;6, &minus;37 % at &minus;4. Drag drops with it.</li>
<li>The area limit stops being active (63&ndash;83 % of the start): without a strong suction peak, a thin section no
longer helps.</li>
<li>With C<sub>p</sub> &ge; &minus;4, run_027 reaches only {dn(R27)}. Run_038 (least loss with downforce &ge; 1.4 and
the same C<sub>p</sub> limit) reaches 1.40, the same within the mesh noise, with its lowest point
{B(R38)['clearance_mm']:.0f} mm above the ground instead of {B(R27)['clearance_mm']:.0f} mm.</li>
</ul>
<p>Loss minimized with C<sub>p</sub> &ge; &minus;4 and a downforce limit (runs 036&ndash;039, at their KKT points,
section 8). On 22 November run_038 had no feasible design and run_039 had stopped 0.1 % short of its downforce limit.</p>
{table(['run', 'downforce limit', '&minus;C<sub>L</sub>c', 'outlet loss (%, corrected)', 'C<sub>D</sub>c', 'lowest C<sub>p</sub> (KS)'],
       [[r[4:7], t, f"{-B(r)['clc']:.2f}", f"{100*B(r)['loss']:.3f}", f"{B(r)['cdc']:.4f}", f"{B(r)['cp_min_ks']:.2f}"] for r, t in sorted(CPML, key=lambda q: q[1])])}
{fig(CMP/'contact_cpmin.png', 'Best feasible sections, on the same grid in every panel.')}
{fig(gifs['run_027_area050CpMin4'], 'run_027 (C<sub>p</sub> &ge; &minus;4).')}
''')

# ------------------------------------------------------------------ 8 KKT
def krow(lab, lim, b, meshes_evals):
    return [lab, lim, f"{-b['clc']:.3f}", f"{b['cdc']:.4f}", f"{100*b['loss']:.3f}", f"{b['cp_min_ks']:.2f}", f"{b['clearance_mm']:.0f}",
            active(b), e(b['psqp_g']), e(b['kkt_residual']), meshes_evals]
KROWS = [krow(r[4:7], t, B(r), f"{R[r]['kkt_rounds']} / {R[r]['kkt_evals']}") for r, t in KKTR]
alt38 = R[R38]['kkt_alternatives'][0]
KROWS.append(krow('038, second point', '1.4', alt38, '(5 mm branch)'))
loss_lo, loss_hi = 100*B(KKTR[0][0])['loss'], 100*B(R38)['loss']
clr = sorted(B(r)['clearance_mm'] for r, _ in KKTR)
S['kkt'] = ('8 Converging to KKT (036&ndash;039)', f'''
<p>Runs 036&ndash;039 minimize the outlet loss with a downforce limit (0.8, 1.0, 1.2, 1.4) and C<sub>p</sub> &ge; &minus;4.
By 22 November none had ended on its own: the stall rule ended every mesh, and run_038 had no feasible design. From 25
to 27 November they were run on, without the stall rule, until PSQP&rsquo;s own optimality test was met: the largest
component of the gradient of the Lagrangian (in the scaled design variables) at most 10<sup>&minus;4</sup>, with every
constraint met to 10<sup>&minus;6</sup>. These are the first-order (KKT) conditions: no small change of the shape that
keeps the active limits can lower the loss.</p>
<h3>Results</h3>
{table(['run', 'downforce limit', '&minus;C<sub>L</sub>c', 'C<sub>D</sub>c', 'loss (%, corrected)', 'lowest C<sub>p</sub> (KS)',
        'lowest point (mm)', 'at a limit (multiplier &gt; 0.001)', 'PSQP G', 'independent check', 'meshes / evaluations in the push'], KROWS)}
<p>&ldquo;Independent check&rdquo;: the KKT residual recomputed from the stored gradients, with sign-constrained
multipliers fitted by least squares (scripts/kkt_residual.py). It agrees with PSQP within the tolerance; for run_038 it
gives {e(B(R38)['kkt_residual'])}, at the tolerance rather than inside it.</p>
{fig(CMP2/'kkt_sections.png', 'The KKT designs.')}
{fig(CMP2/'kkt_results.png', '(a) Loss at the KKT point, corrected measurement. (b) The KKT designs against the drag sweep of section 4 and run_027. (c) Suction-side C<sub>p</sub>.')}
<ul>
<li>All four meet the downforce limit exactly, sit on the 50 % area limit, and put the trailing edge in the top corner
of the box. The lowest point is {clr[0]:.0f}&ndash;{clr[-1]:.0f} mm above the ground, far from the 50 mm floor: with
downforce limited rather than maximized, the least-loss wing stays away from the ground.</li>
<li>The loss rises with the downforce limit, from {loss_lo:.2f} % at 0.8 to {loss_hi:.2f} % at 1.4.</li>
<li>The C<sub>p</sub> limit is active only at 1.2 and 1.4 (039, 038), where the suction side shows the plateau at
&minus;4 seen in section 7. At 0.8 and 1.0 the peak stays at {B(KKTR[0][0])['cp_min_ks']:.1f} and {B(KKTR[1][0])['cp_min_ks']:.1f}.</li>
<li>In the downforce&ndash;drag plane the KKT designs lie on or just below the drag sweep: run_038 makes 1.40 at
C<sub>D</sub>c {B(R38)['cdc']:.4f}, against {dn(R29)} at {B(R29)['cdc']:.4f} for run_029 (drag limit 0.02). Least
loss and least drag give nearly the same design here, as the loss and drag limits did in section 5.</li>
<li>Run_038 has two KKT points, {100*B(R38)['loss']:.3f} % and {100*alt38['loss']:.3f} %, with different thickness
stations at their limits; the lower is reported. With about 830 constraints and up to 75 near a bound, more than one
local optimum is to be expected.</li>
</ul>
<h3>What stopped them, and the fixes</h3>
<p>Near an optimum the loss differs between the trial points of a line search by as little as 10<sup>&minus;7</sup> of
its value. Any error in the loss that jumps by that much between nearby designs stops the optimizer: PSQP rejects the step,
shortens it, and finally gives up (exit &minus;6). Four such errors were found and removed, in this order.</p>
{table(['problem', 'what it did', 'fix'], [
 ['pressure level', 'the loss jumped by up to 3 % between designs 10<sup>&minus;8</sup> m apart. Part of it was the overall pressure level of the field, which the low-Mach solve leaves slightly unconverged and which moves with the design; most of it turned out to be the loss planes (next row)',
  'loss measured against the mass-averaged total pressure on an upstream plane (x = &minus;4 m) of the same solution; solves to 10<sup>&minus;12</sup> instead of 10<sup>&minus;10</sup>. At the KKT designs that plane is within 0.4 Pa of the free-stream value, so the reference itself changes the loss by under 0.0002'],
 ['loss-plane points', 'every other point of the outlet and inlet planes (1.25 mm apart, over 2.5 mm cells) sat exactly on a row of background cell centres. Whether ADflow found those points depended on round-off and on the number of MPI ranks. The loss read 9&ndash;12 % low, by an amount that changed with tiny design changes',
  'points moved to the cell midpoints. The same saved flow then reads the same on 4 and 12 ranks (0.004764) and matches an independent integration of the field (0.00477)'],
 ['far field moving with the wing', 'the mesh deformation (IDWarp) moves every node of the background grid, so the domain boundary moved with the wing, by up to 56 mm within a mesh. The far-field condition treats each boundary face as inflow or outflow by the sign of the local normal velocity, so faces on the top boundary switched as the design changed, and the converged flow jumped (figure below, a). PSQP sat in a &ldquo;pocket&rdquo;: every trial point read worse than the current design',
  'the far-field boundary is held in place in the deformation, and its top is tilted down by 0.2 % toward the outlet so every top face is an outflow. The replay is then smooth (a). Holding it changes the result little near a mesh&rsquo;s start and more further on: 1.3 % less loss and 0.25 % more downforce for run_037 at 30 mm from its mesh&rsquo;s start'],
 ['solver tolerance (038)', 'with solves to 10<sup>&minus;12</sup> the loss still scattered by 3&times;10<sup>&minus;7</sup> of its value when nearby designs were solved in sequence (b): as large as the gain left near the optimum. Both 038 branches ended their line searches there',
  'solves to 10<sup>&minus;13</sup> (scatter 1&ndash;3&times;10<sup>&minus;8</sup>, solves 10&ndash;25 % longer), continued on the same mesh. Both branches then met the test within 2 iterations']])}
{fig(CMP2/'replays.png', 'Replays. (a) The failed line search of run_037 mesh 7, evaluation 172, solved again with the old far field (orange) and with the far field held (blue). (b) Near 038&rsquo;s KKT point, a 0.36 mm segment solved point after point, each solve starting from the previous one, to 10<sup>&minus;12</sup> (orange) and 10<sup>&minus;13</sup> (blue).')}
<p><b>Step length (038).</b> At downforce 1.4 the C<sub>p</sub> limit is strongly curved: the slope of the KS aggregate
changes by 20 % over 0.2 mm. With steps of up to 10 mm (PSQP&rsquo;s XMAX 1.0 in the scaled variables, used by all four
runs) each step&rsquo;s linear prediction put C<sub>p</sub> on its limit while the solved C<sub>p</sub> missed by
0.06&ndash;0.16, and the next step spent itself undoing the miss: 038 cycled, never feasible, for 40 iterations over two
meshes. Steps of at most 0.5 mm cured it (misses below 10<sup>&minus;3</sup>) at the cost of slow travel. Once the path
was straight, branches with 2 mm and 5 mm steps covered the rest faster; both 038 KKT points come from them.</p>
{fig(CMP2/'kkt_history.png', 'The push, 25&ndash;27 November: (a) loss and (b) PSQP&rsquo;s G of every iterate; one line per mesh. Circled numbers: the changes, in order (listed under the plot). The loss definition changes at 1 and 3, so the levels on either side are not comparable.')}
<p>Cost: {PUSH['evals']:,} flow evaluations and {PUSH['h']:.0f} h of run time over the four runs (run at the same time on
12&ndash;24 cores each), about {round(PUSH['core_h'], -2):,.0f} core-hours; run_038 alone {R[R38]['kkt_evals']:,} evaluations.</p>
<h3>What this means for the other results</h3>
<ul>
<li><b>Loss values</b> up to run_035 (sections 3&ndash;7) use the old plane points and read 9&ndash;12 % low, by an
amount that varies with the design. Comparisons within those sections hold roughly. Re-evaluating the saved best designs
with the corrected planes would put them on the same footing as this section.</li>
<li><b>Runs that optimized the loss</b> (019&ndash;021 and 030 with a loss limit, 031&ndash;034 minimizing it)
optimized the flawed measurement, and their line searches stopped on its jumps.</li>
<li><b>Downforce, drag and C<sub>p</sub></b> do not depend on the loss planes. The far-field drift affects them a
little: it grows within a mesh and resets at each new one, which probably explains much of the ~2 % shift of
C<sub>L</sub>c at a new mesh (section 2). With the far field held, three remeshes of run_038 at fixed designs moved
C<sub>L</sub>c by 0.04&ndash;0.6 % and the loss by 0.1&ndash;0.6 %.</li>
</ul>
{fig(gifs[R38], 'run_038, every optimizer iteration from its first mesh to the KKT point (its main chain, then the 2 mm branch), in the newer layout: (a) C<sub>pt</sub>, (b) C<sub>p</sub> with the surface C<sub>p</sub> inset, (c) section with its FFD control points and the surface sensitivity of the loss (maize teeth: the way a surface move lowers it), (d) downforce against drag, (e) every constraint as the share of its limit used (blue: on the limit, maize: past it).')}
''')

# ------------------------------------------------------------------ 9 rotated start, least suction
b42 = B(R42); s42 = json.loads((O/R42/'r00/baseline/result.json').read_text()); l41 = R[R41]['last']
S['rotsuction'] = ('9 Rotated start, least suction (040&ndash;042)', f'''
<h3>Rotated start (040, 041)</h3>
<p>Most remeshes in runs 016&ndash;029 come from the trailing edge turning up into the top corner of the box: the mesh
deformation folds cells as the tail rotates. Runs 040 and 041 start instead from the level NACA 0012 rotated 8.3&deg;
nose down (and scaled by 1.01) so that it still spans the box, its nose touches the floor and its trailing edge sits at
2/3 of the box height. Run_040 has run_016&rsquo;s limits, run_041 run_029&rsquo;s (C<sub>D</sub>c &le; 0.02).</p>
{table(['run', 'start', 'limits', '&minus;C<sub>L</sub>c', 'C<sub>D</sub>c', 'meshes', 'hours'], [
    ['016', 'level', 'no drag limit', dn(R16), dr(R16), meshes(R16), hours(R16)],
    ['040', 'rotated', 'as 016', dn(R40), dr(R40), meshes(R40), hours(R40)],
    ['029', 'level', 'C<sub>D</sub>c &le; 0.02', dn(R29), dr(R29), meshes(R29), hours(R29)],
    ['041', 'rotated', 'as 029', f"{-l41['clc']:.2f} (last iterate, outside the drag limit)", f"{l41['cdc']:.4f}", meshes(R41), hours(R41)]])}
{fig(CMP2/'rotsuction.png', '(a) Level and rotated start. (b) Best designs of 016 and 040, same limits. (c, d) Run_042: start and least-suction design, sections and surface C<sub>p</sub>; dotted in (d): 039, the least-loss design at the same downforce with C<sub>p</sub> &ge; &minus;4.')}
<ul>
<li>Run_040 ends at the same design as run_016 ({-B(R40)['clc']:.3f} against {-B(R16)['clc']:.3f}), with one remesh fewer, and neither of its
remeshes came from a failed mesh. The optimum curls the tail up about 72&deg;, far beyond the start&rsquo;s 8&deg;, so
a rotated start delays mesh failure but does not remove remeshing.</li>
<li>Run_041 was ended by the stall rule after 6 optimizer iterations, still outside its drag limit, at
&minus;C<sub>L</sub>c {-l41['clc']:.2f} against {dn(R29)} for run_029. It is inconclusive and was not continued.</li>
</ul>
<h3>Least suction peak (042)</h3>
<p>Run_042 minimizes the suction peak (raises the lowest C<sub>p</sub>, KS aggregate) with downforce &ge; 1.2, starting
from run_027&rsquo;s best design on a new mesh.</p>
{table(['', '&minus;C<sub>L</sub>c', 'C<sub>D</sub>c', 'lowest C<sub>p</sub> (KS)', 'area (%)', 'lowest point (mm)'], [
    ['start (027 best, new mesh)', f"{-0.5*s42['cl']:.2f}", f"{0.5*s42['cd']:.3f}", f"{s42['cp_min_ks']:.2f}", f"{100*s42['section_area_m2']/A0:.0f}", f"{1e3*s42['minimum_clearance_m']:.0f}"],
    ['042 best', dn(R42), dr(R42), f"{b42['cp_min_ks']:.2f}", area(R42), low(R42)],
    ['039 KKT (least loss, C<sub>p</sub> &ge; &minus;4)', dn(R39), f"{B(R39)['cdc']:.3f}", f"{B(R39)['cp_min_ks']:.2f}", area(R39), low(R39)]])}
<ul>
<li>The peak rises from {s42['cp_min_ks']:.2f} to {b42['cp_min_ks']:.2f} at the downforce limit: the suction side
becomes a flat plateau near &minus;2.2 over most of the chord, ending in a steep pressure rise near the trailing edge.</li>
<li>It costs drag: C<sub>D</sub>c {b42['cdc']:.3f} against {0.5*s42['cd']:.3f} at the start and {B(R39)['cdc']:.3f}
for run_039, which reaches the same downforce with C<sub>p</sub> &ge; &minus;4 at about a seventh of the drag. A mild
peak alone is a poor goal; a C<sub>p</sub> limit with a loss or drag objective does better.</li>
</ul>
''')

# ------------------------------------------------------------------ 10 misc
MISC = [('run_001_baseline', 'large box (x &minus;125&hellip;625, y 50&hellip;350 mm), 10 pairs, x + y moves, tilt free, NACA 0012 at 3&deg;'),
        ('run_002_reducedBox', 'as 001 in the final box'),
        ('run_003_ffd20Pairs', 'as 002 with 20 pairs'),
        ('run_011_naca0010LevelFront', 'as 003, start NACA 0010 at 0&deg;'),
        ('run_007_floor30mm', 'as 003, floor 30 mm'),
        ('run_006_levelNoTilt', 'y moves only, tilt locked, level NACA 0012 start (016 without the area limit)'),
        (R16, 'as 006, area &ge; 50 %')]
def mrow(r, desc):
    b_ = B(r); ar = '&ndash;' if b_['area_pct'] is None else f"{b_['area_pct']:.0f}"
    return [r[4:7], desc, dn(r), dr(r), ar, meshes(r), hours(r)]
a006 = 'about 14'
S['misc'] = ('10 Other runs', f'''
<p>Runs from before run_016 that changed the setup, and what they showed.</p>
{table(['run', 'setup', '&minus;C<sub>L</sub>c', 'C<sub>D</sub>c', 'area (%)', 'meshes', 'hours'], [mrow(r, d) for r, d in MISC])}
<ul>
<li><b>Control points moving in x too</b> (003, tilt also free): 6 % more downforce than y moves only with the tilt
locked (006), same drag. We kept y moves only: simpler, and the y moves still tilt the section (12&deg; in 006). The
leading and trailing edges then stay at the ends of the box.</li>
<li><b>No area limit</b> (006): the section thins to {a006} % of the start. Keeping at least 50 % (016) costs only
3.5 % of the downforce.</li>
<li><b>Other starting points</b>: NACA 0010 at 0&deg; (011) and NACA 0012 at 3&deg; (003) end at the same downforce and
a similar shape. For the loss-objective runs the start mattered (section 5).</li>
<li><b>Box size</b>: the large box (001) gave &minus;C<sub>L</sub>c 3.96 but took 14 h; the final box cut a run to 2&ndash;3 h.</li>
<li><b>More pairs</b>: 20 pairs (003) against 10 (002) gave 1 % more downforce, inside the mesh noise, and wore out
the mesh twice as fast.</li>
</ul>
{fig(CMP/'contact_misc.png', 'Best feasible sections. The x + y runs (003, 011, 007) end with a thin section and a near-vertical flap at the trailing edge.')}
''')

# ------------------------------------------------------------------ 11 conclusion
S['conclusion'] = ('11 Conclusion', f'''
<ul>
<li>From a level NACA 0012, the optimizer finds a thin section with its nose on the box floor and a steep,
Gurney-like turn at the trailing edge: downforce {-s16['clc']:.2f} &rarr; {dn(R16)}.</li>
<li>Drag and loss limits trade downforce smoothly, cheap at first and steep later. Capping loss and capping drag give
the same trade. Under a drag limit the wing rises off the floor.</li>
<li>Lowering the floor helps down to about 25 mm, then hurts, as Zerihan and Zhang found [3]. In free air the same setup
makes {100*B(FR)['clc']/B(R16)['clc']:.0f} % of the downforce.</li>
<li>A C<sub>p</sub> limit flattens the suction peak and costs downforce fast (&minus;37 % at &minus;4).</li>
<li>With the loss as the objective, a downforce limit and C<sub>p</sub> &ge; &minus;4, all four runs (036&ndash;039) now
end at verified first-order optima: the loss rises from {loss_lo:.2f} % at downforce 0.8 to {loss_hi:.2f} % at 1.4
(corrected measurement), the C<sub>p</sub> limit binds from 1.2 up, and the designs lie on or near the drag front.</li>
<li>Getting there took fixing the measurement, not the optimizer: loss-plane points on cell-centre rows (9&ndash;12 % low),
a far-field boundary that moved with the mesh, and the solver tolerance. Loss values up to run_035 should be
re-evaluated before they are compared with the corrected ones.</li>
<li>A start tilted toward the top corner reaches the same optimum as the level start with one remesh fewer; it does not
remove remeshing. Minimizing the suction peak alone flattens it to &minus;2.2 at about seven times the drag of the
least-loss design.</li>
<li>The suction side and the trailing edge set the downforce; the top surface matters little and varies between runs.</li>
<li>Some surfaces come out wavy. A curvature limit would be the next step toward designs that can be built.</li>
<li>A new mesh moved the results of runs up to 035 by about 2 %; with the far field held it moves them by under 1 %.</li>
<li>One optimizer iteration takes about 75 s on 24 cores, two thirds in flow solves; the optimizer&rsquo;s own work is
negligible.</li>
</ul>
<h3>References</h3>
<ol class="refs">
<li>Katz, J. Aerodynamics of race cars. <i>Annual Review of Fluid Mechanics</i> 38, 27&ndash;63, 2006.</li>
<li>Zhang, X., Toet, W., Zerihan, J. Ground effect aerodynamics of race cars. <i>Applied Mechanics Reviews</i> 59(1), 33&ndash;49, 2006.</li>
<li>Zerihan, J., Zhang, X. Aerodynamics of a single element wing in ground effect. <i>Journal of Aircraft</i> 37(6), 1058&ndash;1064, 2000.</li>
<li>Zerihan, J., Zhang, X. Aerodynamics of Gurney flaps on a wing in ground effect. <i>AIAA Journal</i> 39(5), 772&ndash;780, 2001.</li>
<li>Liebeck, R. H. Design of subsonic airfoils for high lift. <i>Journal of Aircraft</i> 15(9), 547&ndash;561, 1978.</li>
<li>Martins, J. R. R. A., Ning, A. <i>Engineering Design Optimization</i>. Cambridge University Press, 2021.</li>
<li>Lambe, A. B., Martins, J. R. R. A. Extensions to the design structure matrix for the description of multidisciplinary design, analysis, and optimization processes. <i>Structural and Multidisciplinary Optimization</i> 46, 273&ndash;284, 2012.</li>
<li>Mader, C. A., Kenway, G. K. W., Yildirim, A., Martins, J. R. R. A. ADflow: an open-source computational fluid dynamics solver for aerodynamic and multidisciplinary optimization. <i>Journal of Aerospace Information Systems</i> 17(9), 508&ndash;527, 2020.</li>
<li>Yildirim, A., Kenway, G. K. W., Mader, C. A., Martins, J. R. R. A. A Jacobian-free approximate Newton&ndash;Krylov startup strategy for RANS simulations. <i>Journal of Computational Physics</i> 397, 108741, 2019.</li>
<li>Kenway, G. K. W., Mader, C. A., He, P., Martins, J. R. R. A. Effective adjoint approaches for computational fluid dynamics. <i>Progress in Aerospace Sciences</i> 110, 100542, 2019.</li>
<li>Secco, N. R., Kenway, G. K. W., He, P., Mader, C. A., Martins, J. R. R. A. Efficient mesh generation and deformation for aerodynamic shape optimization. <i>AIAA Journal</i> 59(4), 1151&ndash;1168, 2021.</li>
<li>Hajdik, H. M., et al. pyGeo: a geometry package for multidisciplinary design optimization. <i>Journal of Open Source Software</i> 8(87), 5319, 2023.</li>
<li>Kenway, G. K. W., Kennedy, G. J., Martins, J. R. R. A. A CAD-free approach to high-fidelity aerostructural optimization. AIAA 2010-9231, 2010.</li>
<li>Wu, N., Kenway, G., Mader, C. A., Jasa, J., Martins, J. R. R. A. pyOptSparse: a Python framework for large-scale constrained nonlinear optimization of sparse systems. <i>Journal of Open Source Software</i> 5(54), 2564, 2020.</li>
<li>Spalart, P. R., Allmaras, S. R. A one-equation turbulence model for aerodynamic flows. AIAA Paper 92-0439, 1992.</li>
</ol>
''')

# ------------------------------------------------------------------ page
nav = ''.join(f'<button data-s="{k}">{v[0]}</button>' for k, v in S.items())
secs = ''.join(f'<section id="{k}"><h2>{v[0]}</h2>{v[1]}</section>' for k, v in S.items())
html = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Front wing in ground effect: shape optimization</title>
<style>
:root {{ --ink:#2b2b2b; --muted:#6b6b6b; --line:#e4e4e4; --blue:#1f62bd; --bg:#ffffff; --tint:#f4f8fd; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font:16px/1.55 -apple-system, "Segoe UI", "DejaVu Sans", Arial, sans-serif; }}
header {{ padding:22px 16px 10px; max-width:1100px; margin:0 auto; }}
header h1 {{ font-size:24px; margin:0 0 4px; }} header p {{ margin:0; color:var(--muted); }}
nav {{ position:sticky; top:0; z-index:5; background:var(--bg); border-bottom:1px solid var(--line); }}
nav div {{ max-width:1100px; margin:0 auto; padding:8px 16px; display:flex; flex-wrap:wrap; gap:6px; }}
nav button {{ font:inherit; font-size:14px; padding:6px 11px; border:1px solid var(--line); background:#fff; border-radius:6px; cursor:pointer; color:var(--ink); }}
nav button.on {{ background:var(--blue); border-color:var(--blue); color:#fff; }}
main {{ max-width:1100px; margin:0 auto; padding:6px 16px 60px; }}
section {{ display:none; }} section.on {{ display:block; }}
h2 {{ font-size:21px; margin:18px 0 8px; }} h3 {{ font-size:17px; margin:22px 0 6px; }}
p, li {{ max-width:860px; }}
figure {{ margin:16px 0 22px; }} figure img {{ display:block; max-width:100%; height:auto; border:1px solid var(--line); border-radius:4px; }}
figcaption {{ color:var(--muted); font-size:14px; margin-top:6px; max-width:860px; }}
table {{ border-collapse:collapse; margin:10px 0 16px; font-size:14.5px; display:block; overflow-x:auto; max-width:100%; }}
th, td {{ border-bottom:1px solid var(--line); padding:5px 10px; text-align:left; vertical-align:top; }}
th {{ background:var(--tint); font-weight:600; }}
svg.xdsm {{ width:100%; max-width:900px; height:auto; display:block; margin:10px 0 18px; }}
.refs li {{ font-size:14.5px; margin-bottom:4px; }}
.pager {{ margin-top:28px; display:flex; justify-content:space-between; }}
.pager button {{ font:inherit; font-size:14px; padding:6px 12px; border:1px solid var(--line); background:#fff; border-radius:6px; cursor:pointer; }}
@media print {{ nav, .pager {{ display:none; }} section {{ display:block; page-break-before:always; }} }}
</style></head>
<body>
<header><h1>Front wing in ground effect: shape optimization</h1>
<p>2-D section in a regulation box, RANS CFD with adjoint gradients (MACH-Aero). Runs 001&ndash;042, November 2026;
updated 27 November.</p></header>
<nav><div>{nav}</div></nav>
<main>{secs}<div class="pager"><button id="prev">&larr; previous</button><button id="next">next &rarr;</button></div></main>
<script>
const keys = [...document.querySelectorAll('nav button')].map(b => b.dataset.s);
function show(k) {{
  document.querySelectorAll('section').forEach(s => s.classList.toggle('on', s.id === k));
  document.querySelectorAll('nav button').forEach(b => b.classList.toggle('on', b.dataset.s === k));
  try {{ history.replaceState(null, '', '#' + k); }} catch (e) {{}}
  window.scrollTo(0, 0);
}}
document.querySelectorAll('nav button').forEach(b => b.addEventListener('click', () => show(b.dataset.s)));
const cur = () => keys.findIndex(k => document.getElementById(k).classList.contains('on'));
document.getElementById('prev').onclick = () => show(keys[Math.max(0, cur() - 1)]);
document.getElementById('next').onclick = () => show(keys[Math.min(keys.length - 1, cur() + 1)]);
show(keys.includes(location.hash.slice(1)) ? location.hash.slice(1) : keys[0]);
</script>
</body></html>'''
OUT.write_text(html)
print(OUT, f'{OUT.stat().st_size/1e6:.1f} MB')

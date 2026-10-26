# How the section shape is controlled, and how it is optimized

Scope: the quasi-2-D NACA 0012 front-wing section of this repo (chord 0.70 m, 100 km/h,
minimum clearance 0.20 m, objective = downforce / drag). The numbers that verify the set-up
(solver settings, mesh level, derivative accuracy) are in [`../study/README.md`](../study/README.md).
Book references are to Martins & Ning, *Engineering Design Optimization* (2022).

## 1. Conventions (from the repo, kept everywhere)

| Item | Convention |
| --- | --- |
| Axes | x downstream, y up, z spanwise. Ground is the plane y = 0. |
| Freestream | Always horizontal (ADflow `alpha = 0`). Incidence is put into the geometry, never into the inflow, because tilting the inflow would tilt it into the road. |
| Incidence α | **Negative α = nose down = leading edge below trailing edge.** Rotation about the quarter chord, c/4. A rotation of the geometry by −α (counter-clockwise) gives incidence α. |
| Clearance h | Distance from the lowest surface point to the ground, measured after rotation. |
| Forces | Airfoil only (`airfoil` family), lift positive up. Downforce coefficient is −C_L. Per metre of span, with reference area chord × extruded span. |
| Ground BC | Ideal slip, i.e. a symmetry plane on the bottom face of the background block. Inlet, outlet and top are `bcfarfield`; the two span faces are symmetry planes. |

## 2. Parameterization: three layers of design variables

The airfoil surface is embedded in a **free-form deformation (FFD) box**, a B-spline volume
handled by pyGeo `DVGeometry`. The surface nodes are parametric points inside that volume, so
moving the box's control points moves the surface smoothly. The box has

* **10 chordwise stations** (index i, LE → TE, slightly wider than the chord),
* **2 rows** (index j, lower and upper), fitted to the section's thickness plus a margin of 0.02 c,
* **2 spanwise planes** (index k, just outside z = 0 and z = span).

Chordwise the B-spline is cubic, so neighbouring stations blend smoothly. Across the
thickness and the span it is linear.

To get a feel for it, open [`ffd_explorer.html`](ffd_explorer.html) in a browser. It
reproduces this box, the embedding, the global/local order and the constraints below, and
lets you drag the control points.

| Layer | Design variables | What moves | How it is set up in pyGeo |
| --- | --- | --- | --- |
| **Incidence** | α (1 global DV, deg) | The whole box rotates rigidly about a spanwise reference axis through c/4 | `addRefAxis("c4", xFraction=0.25, alignIndex="k")` and `addGlobalDV("alpha", ...)`, which sets `rot_z = -α` (repo sign) |
| **Ride height** | Δh (1 global DV, m) | The whole box translates vertically | Same reference axis: `coef[:, 1] += Δh` |
| **Section shape** | 20 local DVs, Δy of each control point (10 stations × 2 rows) | Local camber and thickness | `addSpanwiseLocalDV("shape", "k", axis="y")`: the two spanwise planes are linked, so the section stays 2-D |

Two linear constraints keep the layers independent:

* **LE/TE constraints** (`DVCon.addLeTeConstraints`): at the first and last station, the upper
  and lower control points must move equal and opposite. Without them the local variables
  could also shear (rotate) or shift the section, duplicating α and Δh. That would make the
  optimum non-unique and the Hessian singular.
* pyGeo applies the global variables first (rotation, then translation, of the control
  points) and then adds the local Δy along **global y** (`DVGeometry.update`). So a local
  variable always moves its control point vertically, whatever the current α. At |α| ≤ 4°
  that differs from "normal to the chord" by less than 0.25 %.

From design variables to flow:

```
x = (α, Δh, Δy_1..20) --pyGeo FFD--> surface nodes --IDWarp--> volume nodes of both overset blocks
   --ADflow--> hole cutting / interpolation --> RANS-SA solution --> C_L, C_D on the airfoil
```

IDWarp moves every volume node, weighted by distance from the moved surface. Nodes on the
ground and span symmetry planes stay in those planes. The near-field O-mesh is 0.11 m thick
and the minimum clearance is 0.20 m, so its outer boundary never reaches the ground.

**Overset caveat.** Once the shape moves, the near-field moves relative to the
background. ADflow's `oversetUpdateMode` decides what happens:

* `frozen` keeps the interpolation weights. It is consistent with the adjoint, but stale
  after large moves.
* `full` re-cuts the holes. It is correct, but this re-cut is *not* differentiated.

The plan is `frozen` inside a design iteration and a `full` re-cut whenever the geometry has
moved more than a fraction of a background cell. The derivative study measures how much this
matters.

## 3. Optimization problem

In the book's notation (x design variables, f objective, g ≤ 0 inequality constraints).
Maximizing downforce/drag = (−C_L)/C_D is written as minimizing its negative:

```
minimize    f(x) = C_L(x) / C_D(x)                       (C_L < 0 for downforce)
by varying  α                    −4° ≤ α ≤ +4°              (repo bounds)
            Δh                   −0.05 m ≤ Δh ≤ 0.10 m       (g_1 keeps the gap ≥ 0.20 m)
            Δy_1..20             −0.05 c ≤ Δy_i ≤ 0.05 c    (to be confirmed)
subject to  g_1 = 0.20 m − KS_min(y_surface(x)) ≤ 0       minimum clearance
            g_2 = C_L(x) − C_L,max ≤ 0                    downforce at least the baseline value (to be confirmed)
            g_3 = τ − t_k(x)/t_k,0 ≤ 0, k = 1..33         thickness at 2 % … 99 % c, τ to be confirmed
            g_4 = 1 − A(x)/A_0 ≤ 0                        section area at least the baseline
            g_5 = g_min − g_i(x) ≤ 0, i = 1..10           control gap per FFD station (fold-free box)
            h   = LE/TE linear constraints (above)        linear
```

* **Clearance** (g_1). The lowest surface point is a min over points, which is not smooth
  as the lowest point jumps from one node to another. It is replaced by a
  Kreisselmeier–Steinhauser aggregate of the lower-surface node heights (book §5.7,
  constraint aggregation). Its gradient comes from `DVGeo.totalSensitivity` on the surface
  point set: cheap, no flow solve. Downforce in ground effect grows as the gap closes, so this
  constraint is expected to be **active**. The optimizer then holds the section at exactly
  0.20 m, which settles the repo's open "exact-height repositioning" question without a
  separate rule.
* **Thickness and area** (g_3, g_4) come from `DVCon.addThicknessConstraints2D` and
  `addVolumeConstraint`. They keep the section buildable and keep IDWarp from collapsing
  cells at the TE. The stations reach 99 % c because the finite TE (1.8 mm) is where the
  margin is smallest: the ±0.05 c variable bounds alone allow a crossing anywhere the section
  is thinner than about 70 mm. pyGeo's thickness is the *distance* between embedded point
  pairs, so it is unsigned. A step that pushes a station through zero thickness would read
  as positive again, so g_3 alone does not guard against an inverted section.
* **Control gaps** (g_5) guard against inversion. g_i = y_upper,i − y_lower,i is the vertical
  gap between the two control points of FFD station i. At a given x, the upper and lower
  surface points share the box coordinate u, so the local thickness is
  (v_upper − v_lower) · Σ_i N_i(u) g_i. The first factor is fixed and positive and the
  N_i are non-negative, so g_i ≥ g_min > 0 at every station keeps the FFD map fold-free
  and the thickness positive **everywhere**, not only at the sampled stations. This is
  exact for vertical moves at zero incidence. The incidence rotation is applied rigidly
  before the local moves, so within ±4° it only tilts the argument.

  The constraints are linear in the variables:
  `DVCon.addLinearConstraintsShape(upper, lower, factorA=1, factorB=-1, lower=g_min − g_i,0)`,
  using `DVGeo.getLocalIndex(0)[i, 1, 0]` and `[i, 0, 0]`. SLSQP satisfies linear constraints at
  every iterate from a feasible start, so the flow solver never sees a folded mesh. A failed
  IDWarp/ADflow evaluation is only a back-stop: SLSQP handles failures poorly.
* [`ffd_explorer.html`](ffd_explorer.html) enforces g_1 (as the regulation-box floor), g_3
  and g_5 interactively and shows the signed area.
* **Downforce floor** (g_2). The ratio alone can be raised by shedding both downforce and
  drag. The floor stops the optimizer from trading the wing's purpose for efficiency.

### Algorithm

* **SQP** (book §5.5) through pyOptSparse **SLSQP**, the gradient-based constrained optimizer
  available in the image (no SNOPT; IPOPT, an interior-point method in book §5.6, is the
  alternative). The problem has 22 variables and roughly 15 constraints, all smooth when the
  flow is attached: a good fit for a quasi-Newton SQP.
* **Gradients by the discrete adjoint** (book §6.7). One flow solve plus **two adjoint solves
  (C_L and C_D)** per design iteration, whatever the number of variables. The objective
  gradient follows from the quotient rule, ∇f = (C_D ∇C_L − C_L ∇C_D) / C_D². The geometric
  constraints have analytic Jacobians from pyGeo.
* **Scaling**: α in degrees and Δy in units of 0.01 c, so all variables are O(1). The objective
  is scaled so that |f| ≈ 1 at the start.
* **Warm starts**: each flow solve restarts from the previous state, and the overset re-cut is
  scheduled as above.
* **Verification before optimizing** (this study):
  1. The solver settings converge to the tolerance the derivatives need.
  2. The mesh level sits in the asymptotic range, with a known discretization error.
  3. Adjoint gradients match finite differences (book §6.4) for α, Δh and the local variables.

  The optimum is then re-evaluated on the finer mesh with a `full` overset update.

### Starting point and expected behaviour

The start is NACA 0012 at α = 0 and h = 0.20 m (the repo baseline, downforce/drag ≈ 7.9). A
symmetric section makes downforce at zero incidence only through the Venturi effect under its
thickness. The optimizer is expected to:

* drive α nose-down towards the −4° bound,
* add camber with the suction side towards the ground,
* hold the clearance constraint active.

The limits of this model: the slip ground has no road boundary layer, the flow is fully
turbulent, there are no finite-span effects, and lower-surface separation at small gaps is
only as good as SA predicts it.

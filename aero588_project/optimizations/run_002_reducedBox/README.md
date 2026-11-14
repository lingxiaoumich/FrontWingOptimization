# run_002_reducedBox: baseline, Nov 14 2026

Downforce optimization of the NACA 0012 section over a slip ground at 200 km/h, as
run_001_baseline but in a smaller regulation box, so that a full optimization takes about
2 h instead of about 6 h.

## Setup (r00/command.txt)

- Start: NACA 0012, chord 0.5 m, 3° nose down, lowest point 50 mm above the ground;
  original grid; initial flow from setups/slip_ground_xy5_cd005.
- Box for every surface vertex: x from -0.0123 to 0.512 m (1.05 × the starting length)
  and height from 0.05 to 0.1881 m (75 mm above the starting top). run_001 used
  -0.125 to 0.625 m by 0.05 to 0.35 m.
- Design: 10 FFD upper/lower pairs moving in x and y, plus the incidence (41 variables).
- Constraints: surface gap ≥ 2 mm (vertices 15..124 and mirrors), 8° trailing-edge taper,
  trailing-edge base ≥ 1.2 mm, maximum thickness ≥ 6 % chord, C_D c ≤ 10 (never active).
- PSQP with XMAX 1.0; ADflow with the built-in ank_deep settings and NKSwitchTol 1e-6;
  48 MPI ranks (6 per NUMA node, leaving cores free for the desktop).
- run_chain.py remeshes when the grid wears out and stops when the PSQP iterates gain
  less than 0.005 in C_L c in 15 min.

## Result

| run | grid | solver time | C_L c |
|---|---|---|---|
| r00 | original | 45.7 min | -0.576 → -1.938 |
| r01 | remesh of r00 evaluation 20 | 64.8 min | -1.861 → -2.437 |
| r02 | remesh of r01 evaluation 94 | none | died at start, see below |

Final design: C_L c -2.437, C_D c 0.171 (r01 evaluation 94). The best design inside
every limit to the letter is -2.435 (evaluation 84): most iterates sit a few µm under the
1.2 mm trailing-edge base. The section touches every box wall, the 2 mm gap floor and the
trailing-edge base, with its trailing edge turned up into the top-rear corner.

Share of the final C_L c reached after 60 min: 85 %; 75 min: 95 %; 90 min: 98 %.

## Notes

- r01 stopped after 6 of its last 8 evaluations failed (the warped volume grid could not
  follow the trial shapes). It was still gaining about 0.03 per 15 min, so the design is
  close to, but not fully, converged.
- r02 died before its first solve: its deformed FFD lattice lost 22 leading-edge
  vertices (re-embedding error 0.46 m), which gave a self-crossing start shape.
  remesh_design.py now refuses such a lattice; going further needs a fresh lattice built
  around the new shape.
- Flow solutions (*.cgns) and optimizer databases (*.hst) are not in git.

## Plots

- summary/plots/chain_history, chain_sections: progress, constraints, drag polar, sections.
- summary/plots/chain_iterates.gif: every iterate with Cpt, Cp, the FFD control points and
  C_L c against C_D c (scripts/plot_chain_gif.py).
- r00/plots, r01/plots: the flow at the iterates (scripts/plot_iterate_flows.py).

"""Generate the unrotated NACA 0012 baseline; no CFD is performed."""
import json
import math
from pathlib import Path

root = Path(__file__).resolve().parents[1]
chord = 0.7
clearance = 0.2
n = 200
xs = [(1 - math.cos(math.pi * i / n)) / 2 for i in range(n + 1)]
# Closed trailing edge using coefficient -0.1036.
yt = [5 * 0.12 * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4) for x in xs]
height = clearance + chord * max(yt)
upper = [(chord * x, height + chord * y) for x, y in zip(xs, yt)]
lower = [(chord * x, height - chord * y) for x, y in zip(xs, yt)]
points = list(reversed(upper)) + lower[1:]
assert abs(min(y for _, y in points) - clearance) < 1e-12
assert abs(max(x for x, _ in points) - chord) < 1e-12
(root / 'geometry' / 'naca0012_baseline.dat').write_text('NACA0012 chord=0.7m clearance=0.2m incidence=0deg closed_TE\n' + ''.join(f'{x:.12f} {y:.12f}\n' for x, y in points))
config = {
    'airfoil': 'NACA0012', 'dimension': 'quasi-2D', 'chord_m': chord,
    'speed_kmh': 100, 'speed_ms': 100 / 3.6,
    'clearance_m': clearance, 'clearance_reference': 'lowest surface point',
    'ground_y_m': 0, 'clearance_over_chord': clearance / chord,
    'objective': 'maximize downforce / drag',
    'coefficient_objective_if_upward_lift_positive': '-CL/CD',
    'force_scope': 'airfoil only; exclude ground forces',
    'baseline_geometry_incidence_deg': 0,
    'baseline_incidence_note': 'Geometry datum only; operating incidence is not selected.',
    'trailing_edge': 'closed',
    'atmospheric_conditions': None,
    'incidence_bounds_deg': None, 'shape_bounds': None,
    'minimum_downforce': None,
    'ground_boundary': 'moving road intended; ADflow implementation unverified',
    'clearance_during_optimization': 'minimum 0.20 m; exact-height repositioning policy pending',
    'solver_status': 'unverified; WSL launch failed',
    'cfd_results_available': False
}
config_path = root / 'case_config.json'
if config_path.exists():
    config.update(json.loads(config_path.read_text(encoding='utf-8-sig')))
config_path.write_text(json.dumps(config, indent=2) + '\n')
print(f'Generated {len(points)} baseline points; lowest sampled point = {min(y for _, y in points):.6f} m. No CFD run.')


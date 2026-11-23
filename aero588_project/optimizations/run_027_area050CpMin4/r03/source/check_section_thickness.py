"""Analytical examples and derivatives for the actual-surface thickness routine."""
import json
import numpy as np
from airfoil_geometry import maximum_chord_normal_thickness, boundary_self_intersects


def run_checks():
    rectangle = np.array([[.5,.03],[0,.03],[0,0],[0,-.03],[.5,-.03],[.5,.03]])
    t = maximum_chord_normal_thickness(rectangle, 2, (0,4))
    assert abs(t-.06) < 1e-12
    # Upper/lower material vertices are horizontally offset. Pairing them would
    # return 0.3, although the actual maximum chord-normal thickness is 0.246667.
    diamond = np.array([[1.,0],[.65,.2],[0,0],[.25,-.1],[1.,0]])
    expected = .2 + .1*(1-.65)/(1-.25)
    value, gradient = maximum_chord_normal_thickness(diamond, 2, (0,0), True)
    assert abs(value-expected) < 1e-12
    assert not boundary_self_intersects(diamond)
    theta = .47
    rotation = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    rotated = diamond @ rotation.T + [.9,-.7]
    assert abs(maximum_chord_normal_thickness(rotated,2,(0,0))-expected) < 1e-12
    # A concave section can have multiple disjoint interior intervals along a
    # normal. Do not count the exterior gap between them as airfoil thickness.
    concave = np.array([[3,0],[0,0],[0,1],[1,2],[3,4],[3,3],[1,1],[3,1],[3,0]],dtype=float)
    assert not boundary_self_intersects(concave)
    assert abs(maximum_chord_normal_thickness(concave,1,(0,0))-2) < 1e-12
    assert boundary_self_intersects(np.array([[0,0],[1,1],[0,1],[1,0],[0,0]]))
    assert boundary_self_intersects(np.array([[0,0],[1,0],[1,0],[0,1],[0,0]]))
    rng = np.random.default_rng(638)
    errors = []
    for _ in range(8):
        direction = rng.normal(size=diamond.shape)
        direction[-1] = direction[0]  # preserve exact polygon closure
        direction /= np.linalg.norm(direction)
        step = 1e-7
        fd = (maximum_chord_normal_thickness(diamond+step*direction,2,(0,0))-
              maximum_chord_normal_thickness(diamond-step*direction,2,(0,0)))/(2*step)
        errors.append(abs(fd-np.sum(gradient*direction)))
    assert max(errors) < 1e-7, errors
    return {'rectangle_thickness_m':float(t),
            'offset_upper_lower_vertices_thickness_m':float(value),
            'offset_vertices_incorrect_paired_estimate_m':.3,
            'rigid_rotation_and_translation_invariant':True,
            'disjoint_intervals_exclude_exterior_gaps':True,
            'crossed_and_collapsed_sections_rejected':True,
            'max_surface_gradient_error':float(max(errors))}


if __name__ == '__main__':
    print(json.dumps(run_checks(),indent=2))

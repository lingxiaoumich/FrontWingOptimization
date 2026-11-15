"""Thickness of a closed, piecewise-linear CFD section, including x deformation.

Thickness is the longest continuous interval inside the polygon on a line normal
to the current chord. It is evaluated at projected surface vertices, where the
maximum of a piecewise-linear thickness function must occur. Multiple interior
intervals are treated separately, so an overhang does not include exterior gaps.
The active vertex/edge feature gives a piecewise-smooth geometry derivative.
"""
import numpy as np


def _section(points, leading_edge, trailing_edge):
    xy = np.asarray(points)[:, :2]
    if not np.isfinite(xy).all():
        raise ValueError('Nonfinite section coordinates')
    chord = (xy[trailing_edge[0]] + xy[trailing_edge[1]]) / 2 - xy[leading_edge]
    length = np.sqrt(np.dot(chord, chord))
    if float(np.real(length)) < 1e-10:
        raise ValueError('Degenerate section chord')
    tangent = chord / length
    normal = np.array([-tangent[1], tangent[0]])
    # A duplicate closing vertex is not a separate polygon vertex.
    closed = np.linalg.norm(np.real(xy[-1] - xy[0])) < 1e-12
    polygon = xy[:-1] if closed else xy
    relative = polygon - xy[leading_edge]
    return relative @ tangent, relative @ normal


def boundary_self_intersects(points):
    """Detect nonadjacent edge crossings/touches and collapsed boundary edges."""
    xy = np.real(np.asarray(points)[:, :2])
    if not np.isfinite(xy).all():
        return True
    if np.linalg.norm(xy[-1] - xy[0]) < 1e-12:
        xy = xy[:-1]
    end = np.roll(xy, -1, axis=0)
    edge = end - xy
    if np.any(np.linalg.norm(edge, axis=1) < 1e-12):
        return True

    def cross(a, b):
        return a[..., 0]*b[..., 1] - a[..., 1]*b[..., 0]

    delta = xy[None, :, :] - xy[:, None, :]
    denominator = cross(edge[:, None, :], edge[None, :, :])
    nonparallel = np.abs(denominator) > 1e-20
    safe_denominator = np.where(nonparallel, denominator, 1.)
    t = cross(delta, edge[None, :, :]) / safe_denominator
    u = cross(delta, edge[:, None, :]) / safe_denominator
    pairs = np.triu(np.ones(denominator.shape, dtype=bool), k=2)
    pairs[0, -1] = False  # first and last edges are adjacent
    touching = nonparallel & (t >= -1e-12) & (t <= 1+1e-12) & (u >= -1e-12) & (u <= 1+1e-12)
    # Handle collinear overlap separately, including nonadjacent vertex contact.
    collinear = (~nonparallel) & (np.abs(cross(delta, edge[:, None, :])) < 1e-20)
    overlap_x = np.minimum(np.maximum(xy[:, None, 0], end[:, None, 0]),
                           np.maximum(xy[None, :, 0], end[None, :, 0])) >= \
                np.maximum(np.minimum(xy[:, None, 0], end[:, None, 0]),
                           np.minimum(xy[None, :, 0], end[None, :, 0])) - 1e-13
    overlap_y = np.minimum(np.maximum(xy[:, None, 1], end[:, None, 1]),
                           np.maximum(xy[None, :, 1], end[None, :, 1])) >= \
                np.maximum(np.minimum(xy[:, None, 1], end[:, None, 1]),
                           np.minimum(xy[None, :, 1], end[None, :, 1])) - 1e-13
    return bool(np.any(pairs & (touching | (collinear & overlap_x & overlap_y))))


def _feature_value(points, feature, leading_edge, trailing_edge):
    s, h = _section(points, leading_edge, trailing_edge)
    vertex, lower_edge, upper_edge = feature

    def height(edge):
        other = (edge + 1) % len(s)
        fraction = (s[vertex] - s[edge]) / (s[other] - s[edge])
        return h[edge] + fraction*(h[other] - h[edge])

    return height(upper_edge) - height(lower_edge)


def maximum_chord_normal_thickness(points, leading_edge=136,
                                   trailing_edge=(0, 272), return_gradient=False):
    """Return maximum physical thickness, optionally its surface-coordinate Jacobian.

    The caller rejects self-intersecting sections before CFD. Branch selection
    uses real coordinates; complex arithmetic is retained for the active value.
    Like a maximum over stations, the function can have a kink when its active
    location changes. At a tie, a derivative of the selected branch is returned.
    """
    s, h = _section(points, leading_edge, trailing_edge)
    sr, hr = np.real(s), np.real(h)
    order = np.argsort(sr, kind='stable')
    tolerance = 1e-12 * max(1., float(np.ptp(sr)))
    events = [int(order[0])]
    for index in order[1:]:
        if sr[index] - sr[events[-1]] > tolerance:
            events.append(int(index))
    other = np.roll(np.arange(len(s)), -1)
    low, high = np.minimum(sr, sr[other]), np.maximum(sr, sr[other])
    best, feature = -np.inf, None
    for left, right in zip(events[:-1], events[1:]):
        mid = (sr[left] + sr[right]) / 2
        edges = np.flatnonzero((low < mid) & (high > mid))
        if len(edges) % 2:
            raise ValueError('Invalid section: odd number of line/surface intersections')
        heights = hr[edges] + (mid - sr[edges]) / (sr[other[edges]] - sr[edges]) * (hr[other[edges]] - hr[edges])
        edges = edges[np.argsort(heights, kind='stable')]
        # Even/odd filling gives each continuous interval inside the section.
        for bottom, top in edges.reshape(-1, 2):
            for vertex in [left, right]:
                candidate = (vertex, int(bottom), int(top))
                value = float(np.real(_feature_value(points, candidate, leading_edge, trailing_edge)))
                if value > best:
                    best, feature = value, candidate
    if feature is None or best <= 0:
        raise ValueError('Invalid section: no positive interior thickness')
    value = _feature_value(points, feature, leading_edge, trailing_edge)
    if not return_gradient:
        return value
    gradient = np.zeros(np.asarray(points).shape)
    # Only the active vertex, intersected edges, and chord endpoints contribute.
    vertex, bottom, top = feature
    active = {vertex, bottom, (bottom+1) % len(s), top, (top+1) % len(s),
              leading_edge, *trailing_edge}
    for i in active:
        for axis in range(2):
            perturbed = np.asarray(points, dtype=complex).copy()
            perturbed[i, axis] += 1e-30j
            gradient[i, axis] = np.imag(_feature_value(
                perturbed, feature, leading_edge, trailing_edge)) / 1e-30
    return float(np.real(value)), gradient

"""Pure geometry helpers shared by the 2D and 3D lattice previews."""

from __future__ import annotations

import math


def lattice_graph(lattice_type, lattice_constant_nm, half_extent_nm):
    """Return clipped helix centres and nearest-neighbour lattice bonds.

    ``lattice_constant_nm`` is the Square nearest-neighbour spacing and the
    Kagome triangular Bravais constant used by the prediction model.  Kagome
    nearest neighbours are therefore separated by half that value.
    """
    lattice_type = str(lattice_type or "square")
    spacing = max(float(lattice_constant_nm), 1e-9)
    half_extent = max(float(half_extent_nm), spacing)
    raw_points = []
    if lattice_type == "kagome":
        root3 = math.sqrt(3.0)
        basis = ((0.0, 0.0), (.5, 0.0), (.25, root3/4.0))
        count = int(math.ceil(2.0*half_extent/spacing))+4
        for row in range(-count, count+1):
            for column in range(-count, count+1):
                origin_x = spacing*(column+.5*row)
                origin_y = spacing*root3*.5*row
                for basis_x, basis_y in basis:
                    x = origin_x+basis_x*spacing
                    y = origin_y+basis_y*spacing
                    if abs(x) <= half_extent+1e-7 and \
                            abs(y) <= half_extent+1e-7:
                        raw_points.append((x, y))
        neighbour_distance = spacing*.5
    else:
        count = int(math.ceil(half_extent/spacing))+1
        for row in range(-count, count+1):
            for column in range(-count, count+1):
                x, y = column*spacing, row*spacing
                if abs(x) <= half_extent+1e-7 and \
                        abs(y) <= half_extent+1e-7:
                    raw_points.append((x, y))
        neighbour_distance = spacing

    # Adjacent primitive cells can emit the same boundary node.  Quantizing
    # at sub-picometre precision removes only exact construction duplicates.
    point_lookup = {}
    for x, y in raw_points:
        point_lookup.setdefault((round(x, 9), round(y, 9)), (x, y))
    points = sorted(point_lookup.values(), key=lambda value: (value[1], value[0]))

    # A spatial hash avoids an O(n^2) search for the 100-cell 2D preview.
    bucket_size = max(neighbour_distance*1.05, 1e-9)
    buckets = {}
    for index, (x, y) in enumerate(points):
        key = (math.floor(x/bucket_size), math.floor(y/bucket_size))
        buckets.setdefault(key, []).append(index)
    tolerance = max(1e-7, neighbour_distance*.035)
    edges = []
    for first, (x, y) in enumerate(points):
        key_x = math.floor(x/bucket_size)
        key_y = math.floor(y/bucket_size)
        for bucket_x in range(key_x-1, key_x+2):
            for bucket_y in range(key_y-1, key_y+2):
                for second in buckets.get((bucket_x, bucket_y), ()):
                    if second <= first:
                        continue
                    other_x, other_y = points[second]
                    distance = math.hypot(other_x-x, other_y-y)
                    if abs(distance-neighbour_distance) <= tolerance:
                        edges.append((first, second))
    return points, edges


def rotated_graph_outside_square(points, edges, rotation_deg,
                                 square_half_extent_nm):
    """Rotate a lattice graph and its central Square-Seed footprint together.

    Each SST layer shares the in-plane coordinate frame of the Seed support
    attached to that layer.  The exclusion test must therefore happen before
    rotation.  Testing the already-rotated points against an axis-aligned
    square leaves the hole fixed in the lower-layer frame and visibly cuts
    into the upper Seed/SST interface at large twist angles.

    Edges touching a removed helix are removed as well, so no SST bond enters
    the reserved Seed area.
    """
    angle = math.radians(float(rotation_deg))
    cosine, sine = math.cos(angle), math.sin(angle)
    half_extent = max(0.0, float(square_half_extent_nm))
    rotated = [
        (cosine*x-sine*y, sine*x+cosine*y)
        for x, y in points
    ]
    visible_indices = [
        index for index, (x, y) in enumerate(points)
        if not (abs(x) <= half_extent and abs(y) <= half_extent)
    ]
    remap = {old: new for new, old in enumerate(visible_indices)}
    visible_points = [rotated[index] for index in visible_indices]
    visible_edges = [
        (remap[first], remap[second])
        for first, second in edges
        if first in remap and second in remap
    ]
    return visible_points, visible_edges


def rotated_outer_boundary_points(points, rotation_deg, field_half_extent_nm,
                                  boundary_width_nm):
    """Return only the rotated outer ring of a square-clipped lattice field."""
    half_extent = max(0.0, float(field_half_extent_nm))
    cutoff = max(0.0, half_extent-max(0.0, float(boundary_width_nm)))
    angle = math.radians(float(rotation_deg))
    cosine, sine = math.cos(angle), math.sin(angle)
    return [
        (cosine*x-sine*y, sine*x+cosine*y)
        for x, y in points
        if max(abs(x), abs(y)) >= cutoff
    ]

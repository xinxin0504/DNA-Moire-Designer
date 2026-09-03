#!/usr/bin/env python3
"""NumPy-only TEM/FFT analyzer executed outside the Qt process.

The cadnano virtualenv intentionally stays small and does not ship NumPy.
The GUI converts any Qt-readable image to PGM and invokes this worker with the
Homebrew Python that is already available on the target Mac.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from moire_runtime import tool_executable


def read_pgm(filename):
    data = Path(filename).read_bytes()
    if not data.startswith(b"P5"):
        raise ValueError("Only GUI-generated P5 PGM input is supported.")
    position = 2
    tokens = []
    while len(tokens) < 3:
        while position < len(data) and data[position] in b" \t\r\n":
            position += 1
        if position < len(data) and data[position] == 35:
            while position < len(data) and data[position] not in b"\r\n":
                position += 1
            continue
        start = position
        while position < len(data) and data[position] not in b" \t\r\n":
            position += 1
        tokens.append(int(data[start:position]))
    width, height, maximum = tokens
    while position < len(data) and data[position] in b" \t\r\n":
        position += 1
    dtype = np.uint8 if maximum < 256 else ">u2"
    image = np.frombuffer(data[position:], dtype=dtype, count=width * height)
    return image.reshape((height, width)).astype(np.float64)


def write_pgm(filename, image):
    """Write an 8-bit grayscale image without requiring Pillow."""
    array = np.asarray(image, dtype=np.float64)
    array = np.clip(array, 0.0, 255.0).astype(np.uint8)
    height, width = array.shape
    Path(filename).write_bytes(
        ("P5\n%d %d\n255\n" % (width, height)).encode("ascii") +
        array.tobytes())


def _gaussian_filter_fft(array, sigma):
    fy = np.fft.fftfreq(array.shape[0])[:, None]
    fx = np.fft.fftfreq(array.shape[1])[None, :]
    transfer = np.exp(
        -2.0 * np.pi * np.pi * sigma * sigma * (fx * fx + fy * fy))
    return np.real(np.fft.ifft2(
        np.fft.fft2(array) * transfer)).astype(np.float64)


def _write_selected_spot_ifft(image, spots, filename):
    """Phase-preserving IFFT for an explicitly validated set of FFT spots."""
    gray = np.asarray(image, dtype=np.float64) / 255.0
    height, width = gray.shape
    yy, xx = np.mgrid[:height, :width]
    full_spectrum = np.fft.fftshift(np.fft.fft2(
        gray - float(np.mean(gray))))
    aperture_mask = np.zeros((height, width), dtype=np.float64)
    for peak in spots:
        theta = math.radians(float(peak.get("angle", 0.0)))
        cosine, sine = math.cos(theta), math.sin(theta)
        dx, dy = xx - float(peak["x"]), yy - float(peak["y"])
        along = cosine * dx + sine * dy
        across = -sine * dx + cosine * dy
        sigma_x = max(1.35, float(peak.get("rx", 4.5)) / 1.9)
        sigma_y = max(1.05, float(peak.get("ry", 3.2)) / 1.9)
        distance = (along * along / (sigma_x * sigma_x) +
                    across * across / (sigma_y * sigma_y))
        aperture_mask = np.maximum(aperture_mask, np.exp(-.5 * distance))
    reconstruction = np.real(np.fft.ifft2(np.fft.ifftshift(
        full_spectrum * aperture_mask)))
    limit = float(np.percentile(np.abs(reconstruction), 99.5))
    reconstruction = np.clip(
        .5 + .5 * reconstruction / max(limit, 1e-12), 0, 1)
    reconstruction = np.power(reconstruction, .92)
    write_pgm(filename, reconstruction * 255.0)
    return str(filename)


def _reflection_selection_summary(spots):
    """Count accepted Friedel pairs by crystallographic reflection class."""
    classes, families = {}, {}
    for index in range(0, len(spots), 2):
        point = spots[index]
        reflection_class = str(point.get("reflection_class", "unclassified"))
        family = str(point.get("basis_family", point.get(
            "lattice_role", "unclassified")))
        classes[reflection_class] = classes.get(reflection_class, 0) + 1
        families[family] = families.get(family, 0) + 1
    return {"pair_count_by_class": classes,
            "pair_count_by_basis": families}


def _refine_fft_spot_estimates(image, estimates):
    """Snap coarse secondary-lattice peaks to clear full-resolution maxima."""
    gray = np.asarray(image, dtype=np.float64) / 255.0
    height, width = gray.shape
    window = np.outer(np.hanning(height), np.hanning(width))
    magnitude = np.log1p(np.abs(np.fft.fftshift(np.fft.fft2(
        (gray - float(np.mean(gray))) * window))))
    search = max(6, int(round(min(width, height) * 0.009)))
    refined = []
    for estimate in estimates:
        cx, cy = int(round(float(estimate["x"]))), int(round(float(estimate["y"])))
        x0, x1 = max(0, cx-search), min(width, cx+search+1)
        y0, y1 = max(0, cy-search), min(height, cy+search+1)
        patch = magnitude[y0:y1, x0:x1]
        if patch.size == 0:
            continue
        peak_y, peak_x = np.unravel_index(int(np.argmax(patch)), patch.shape)
        px, py = x0 + int(peak_x), y0 + int(peak_y)
        fit = max(4, int(round(search * .55)))
        fx0, fx1 = max(0, px-fit), min(width, px+fit+1)
        fy0, fy1 = max(0, py-fit), min(height, py+fit+1)
        local = magnitude[fy0:fy1, fx0:fx1]
        background = float(np.percentile(local, 45.0))
        weights = np.maximum(local-background, 0.0) ** 1.5
        total = float(np.sum(weights))
        if total > 0:
            gy, gx = np.mgrid[fy0:fy1, fx0:fx1]
            precise_x = float(np.sum(weights*gx)/total)
            precise_y = float(np.sum(weights*gy)/total)
        else:
            precise_x, precise_y = float(px), float(py)
        refined.append({
            "x": round(precise_x, 3), "y": round(precise_y, 3),
            "rx": max(3.5, fit*.62), "ry": max(2.5, fit*.44),
            "angle": 0.0, "score": float(estimate.get("score", 0.0)),
            "lattice_role": "secondary_a",
        })
    return refined


def _indexed_square_spots(assets, basis_families, source_shape,
                          residual_tolerance=.20, score_fraction=.42):
    """Return clear peaks belonging to validated reciprocal bases.

    Every accepted Friedel pair must be expressible as

        G(h,k) = h*g1 + k*g2,  h,k integers.

    For Square this includes the familiar 45-degree (h,h)/(h,-h) diagonal
    shell.  Honeycomb/Kagome/hexagonal families use their measured oblique
    reciprocal basis, so their legal higher-order angles emerge from integer
    indices rather than from a hard-coded 45-degree rule.  Missing or weak
    theoretical reflections are never synthesized.
    """
    visible = assets.get("all_visible_spots") or []
    if not basis_families or len(visible) < 2:
        return []
    height, width = source_shape
    center = assets.get("center") or [width/2.0, height/2.0]
    cx, cy = float(center[0]), float(center[1])
    pair_scores = [float(visible[index].get("score", 0.0))
                   for index in range(0, len(visible), 2)]
    threshold = max(pair_scores) * float(score_fraction) if pair_scores else 0.0
    bases = []
    for family_index, family in enumerate(basis_families):
        angles = family.get("axis_angles_deg") or []
        periods = family.get("axis_periods_px") or []
        if len(periods) < 2 and family.get("lattice_constant_px"):
            periods = [family["lattice_constant_px"]] * 2
        if len(angles) < 2 or len(periods) < 2:
            continue
        columns = []
        for angle, period in zip(angles[:2], periods[:2]):
            radians = math.radians(float(angle))
            period = float(period)
            if period <= 0:
                columns = []
                break
            columns.append((math.cos(radians)/period,
                            math.sin(radians)/period))
        if len(columns) != 2:
            continue
        matrix = np.asarray([[columns[0][0], columns[1][0]],
                             [columns[0][1], columns[1][1]]], dtype=float)
        if abs(float(np.linalg.det(matrix))) < 1e-9:
            continue
        bases.append({
            "matrix": matrix,
            "symmetry": str(family.get("symmetry") or "Square"),
            "family_index": int(family_index),
            "family_id": str(family.get("family_id", "square_%d" % (
                family_index + 1))),
            "lattice_role": str(family.get("lattice_role", "helix_lattice")),
            "residual_tolerance": float(family.get(
                "residual_tolerance", residual_tolerance)),
        })
    if not bases:
        return []
    indexed = {}
    for index in range(0, len(visible), 2):
        pair = visible[index:index+2]
        if len(pair) < 2:
            continue
        score = float(pair[0].get("score", 0.0))
        if score < threshold:
            continue
        fx = (float(pair[0]["x"])-cx) / max(float(width), 1.0)
        fy = -(float(pair[0]["y"])-cy) / max(float(height), 1.0)
        best = None
        for basis in bases:
            coordinates = np.linalg.solve(
                basis["matrix"], np.asarray([fx, fy], dtype=float))
            h, k = int(round(coordinates[0])), int(round(coordinates[1]))
            residual = math.hypot(coordinates[0]-h, coordinates[1]-k)
            if (h == 0 and k == 0 or
                    residual > basis["residual_tolerance"]):
                continue
            candidate = (residual, basis, h, k, coordinates)
            if best is None or candidate[0] < best[0]:
                best = candidate
        if best is None:
            continue
        residual, basis, h, k, coordinates = best
        # Friedel pairs (h,k) and (-h,-k) are one indexed reflection pair.
        if h < 0 or (h == 0 and k < 0):
            h, k = -h, -k
        key = (basis["family_id"], h, k)
        old = indexed.get(key)
        if old is None or score > old[0]:
            indexed[key] = (score, residual, pair, basis)
    selected = []
    for (family_id, h, k), (unused_score, residual, pair, basis) in sorted(
            indexed.items(), key=lambda item: (
                item[0][0], item[0][1]**2+item[0][2]**2,
                item[0][1], item[0][2])):
        symmetry = basis.get("symmetry", "Square")
        if h == 0 or k == 0 or (symmetry != "Square" and h == -k):
            reflection_class = "axis"
        elif symmetry == "Square" and abs(h) == abs(k):
            reflection_class = "diagonal_45"
        elif symmetry != "Square":
            reflection_class = "oblique_integer_order"
        else:
            reflection_class = "mixed_integer"
        for point in pair:
            marked = dict(point)
            marked.update({"lattice_role": basis["lattice_role"],
                           "basis_family": family_id,
                           "miller_h": h, "miller_k": k,
                           "index_residual": residual,
                           "reflection_class": reflection_class,
                           "radial_order": math.hypot(h, k)})
            selected.append(marked)
    return selected


def _indexed_supercell_spots(assets, pore_lattice, source_shape):
    """Return clear Square reflections indexed in a supercell basis."""
    single = assets.get("single_lattice") or {}
    angles = (pore_lattice.get("axis_angles_deg") or
              single.get("axis_angles_deg") or [])
    periods = pore_lattice.get("axis_periods_px") or []
    return _indexed_square_spots(
        assets, [{
            "family_id": "square_supercell",
            "lattice_role": "square_supercell_reflection",
            "axis_angles_deg": angles,
            "axis_periods_px": periods,
            "residual_tolerance": .12,
        }], source_shape, residual_tolerance=.12, score_fraction=.42)


def _local_fft_disk_mean(array, x, y, radius=2.2):
    height, width = array.shape
    x0, x1 = max(0, int(x-radius-1)), min(width, int(x+radius+2))
    y0, y1 = max(0, int(y-radius-1)), min(height, int(y+radius+2))
    yy, xx = np.mgrid[y0:y1, x0:x1]
    mask = (xx-x)**2 + (yy-y)**2 <= radius**2
    values = array[y0:y1, x0:x1][mask]
    return float(np.mean(values)) if values.size else 0.0


def _smooth_boundary_segments(segments, iterations=2):
    """Join grid-edge segments and return antialiased-looking smooth chains."""
    if not segments:
        return [], None
    adjacency = {}
    edges = set()
    for x1, y1, x2, y2 in segments:
        first = (round(float(x1), 4), round(float(y1), 4))
        second = (round(float(x2), 4), round(float(y2), 4))
        if first == second:
            continue
        edge = tuple(sorted((first, second)))
        edges.add(edge)
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)
    chains = []
    unused = set(edges)
    while unused:
        available_degree = {}
        for edge in unused:
            for point in edge:
                available_degree[point] = available_degree.get(point, 0) + 1
        starts = [point for point, degree in available_degree.items()
                  if degree == 1]
        start = starts[0] if starts else next(iter(unused))[0]
        chain, current, previous = [start], start, None
        while True:
            options = []
            for neighbor in adjacency.get(current, ()):
                edge = tuple(sorted((current, neighbor)))
                if edge in unused:
                    options.append(neighbor)
            if not options:
                break
            if previous is not None and len(options) > 1:
                # Continue as straight as possible at a rare grid junction.
                old_angle = math.atan2(current[1]-previous[1],
                                       current[0]-previous[0])
                options.sort(key=lambda point: abs(math.atan2(
                    math.sin(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle),
                    math.cos(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle))))
            neighbor = options[0]
            unused.remove(tuple(sorted((current, neighbor))))
            chain.append(neighbor)
            previous, current = current, neighbor
            if current == start:
                break
        if len(chain) >= 2:
            chains.append(chain)
    smoothed_chains = []
    for chain in chains:
        closed = len(chain) > 2 and chain[0] == chain[-1]
        points = chain[:-1] if closed else chain
        for unused_iteration in range(iterations):
            if len(points) < 3:
                break
            refined = [] if closed else [points[0]]
            count = len(points) if closed else len(points)-1
            for index in range(count):
                first = points[index]
                second = points[(index+1) % len(points)]
                refined.extend((
                    (.75*first[0]+.25*second[0],
                     .75*first[1]+.25*second[1]),
                    (.25*first[0]+.75*second[0],
                     .25*first[1]+.75*second[1]),
                ))
            if not closed:
                refined.append(points[-1])
            points = refined
        if closed and points:
            points = points + [points[0]]
        smoothed_chains.append(points)
    output = []
    for chain in smoothed_chains:
        output.extend([[float(first[0]), float(first[1]),
                        float(second[0]), float(second[1])]
                       for first, second in zip(chain[:-1], chain[1:])])
    longest = max(smoothed_chains, key=lambda chain: sum(
        math.hypot(second[0]-first[0], second[1]-first[1])
        for first, second in zip(chain[:-1], chain[1:])), default=[])
    label = None
    if longest:
        label = longest[len(longest)//2]
    return output, label


def _extend_open_boundary_to_frame(segments, width, height, maximum_gap):
    """Extend only near-frame open endpoints to the image boundary.

    Local FFT windows cannot be centred closer than half a patch to an image
    edge.  Consequently, a real domain seam that exits the image is otherwise
    drawn with two short missing end pieces.  This function adds those pieces
    before TEM evidence refinement.  Closed/internal contours are untouched.
    """
    if not segments:
        return segments
    adjacency = {}
    for x1, y1, x2, y2 in segments:
        first = (round(float(x1), 4), round(float(y1), 4))
        second = (round(float(x2), 4), round(float(y2), 4))
        if first == second:
            continue
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)
    endpoints = [point for point, neighbours in adjacency.items()
                 if len(neighbours) == 1]
    if not endpoints:
        return segments
    extensions = []
    frame_width = max(0.0, float(width)-1.0)
    frame_height = max(0.0, float(height)-1.0)
    for endpoint in endpoints:
        x, y = float(endpoint[0]), float(endpoint[1])
        distances = [(x, (0.0, y)),
                     (frame_width-x, (frame_width, y)),
                     (y, (x, 0.0)),
                     (frame_height-y, (x, frame_height))]
        gap, target = min(distances, key=lambda item: item[0])
        if gap < -1e-6 or gap > float(maximum_gap):
            continue
        # Do not continue the final local-FFT cell orthogonally.  The last
        # cell centre is commonly 1--2 strides away from the image frame and
        # its final edge is horizontal/vertical even when the physical seam
        # is oblique.  Estimate the recent global trend from the connected
        # component and extrapolate that trend to the selected frame edge.
        component = set()
        pending = [endpoint]
        while pending:
            point = pending.pop()
            if point in component:
                continue
            component.add(point)
            pending.extend(adjacency.get(point, ()))
        target_x, target_y = float(target[0]), float(target[1])
        trend_span = max(48.0, float(maximum_gap)*2.5)
        if target_y in (0.0, frame_height):
            slopes = []
            for point_x, point_y in component:
                delta = y-float(point_y)
                if (abs(delta) < max(24.0, float(maximum_gap)*.35) or
                        abs(delta) > trend_span):
                    continue
                slopes.append((x-float(point_x))/delta)
            if slopes:
                slope = float(np.median(np.clip(slopes, -2.0, 2.0)))
                target_x = float(np.clip(
                    x+slope*(target_y-y), 0.0, frame_width))
        else:
            slopes = []
            for point_x, point_y in component:
                delta = x-float(point_x)
                if (abs(delta) < max(24.0, float(maximum_gap)*.35) or
                        abs(delta) > trend_span):
                    continue
                slopes.append((y-float(point_y))/delta)
            if slopes:
                slope = float(np.median(np.clip(slopes, -2.0, 2.0)))
                target_y = float(np.clip(
                    y+slope*(target_x-x), 0.0, frame_height))
        target = (target_x, target_y)
        length = math.hypot(target[0]-x, target[1]-y)
        if length <= 1e-6 or length > float(maximum_gap)*2.75:
            continue
        count = max(1, int(math.ceil(length/max(8.0,
                                                float(maximum_gap)/8.0))))
        points = [(x+(target[0]-x)*index/count,
                   y+(target[1]-y)*index/count)
                  for index in range(count+1)]
        extensions.extend([[float(first[0]), float(first[1]),
                            float(second[0]), float(second[1])]
                           for first, second in zip(points[:-1], points[1:])])
    return list(segments)+extensions


def _ordered_boundary_chains(segments):
    """Return connected boundary segments as ordered point chains."""
    adjacency, unused = {}, set()
    for x1, y1, x2, y2 in segments or []:
        first = (round(float(x1), 4), round(float(y1), 4))
        second = (round(float(x2), 4), round(float(y2), 4))
        if first == second:
            continue
        edge = tuple(sorted((first, second)))
        unused.add(edge)
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)
    chains = []
    while unused:
        degree = {}
        for edge in unused:
            for point in edge:
                degree[point] = degree.get(point, 0)+1
        endpoints = [point for point, value in degree.items() if value == 1]
        start = endpoints[0] if endpoints else next(iter(unused))[0]
        chain, current, previous = [start], start, None
        while True:
            options = [neighbor for neighbor in adjacency.get(current, ())
                       if tuple(sorted((current, neighbor))) in unused]
            if not options:
                break
            if previous is not None and len(options) > 1:
                old_angle = math.atan2(current[1]-previous[1],
                                       current[0]-previous[0])
                options.sort(key=lambda point: abs(math.atan2(
                    math.sin(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle),
                    math.cos(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle))))
            following = options[0]
            unused.remove(tuple(sorted((current, following))))
            chain.append(following)
            previous, current = current, following
            if current == start:
                break
        if len(chain) >= 2:
            chains.append([(float(x), float(y)) for x, y in chain])
    return chains


def _point_in_polygon(point, polygon):
    """Even/odd point-in-polygon test without optional geometry packages."""
    x, y = float(point[0]), float(point[1])
    inside = False
    for first, second in zip(polygon, polygon[1:]+polygon[:1]):
        x1, y1 = first
        x2, y2 = second
        if ((y1 > y) != (y2 > y)):
            crossing = (x2-x1)*(y-y1)/(y2-y1)+x1
            if x < crossing:
                inside = not inside
    return inside


def _polygon_area(polygon):
    return abs(sum(first[0]*second[1]-second[0]*first[1]
                   for first, second in zip(
                       polygon, polygon[1:]+polygon[:1])))/2.0


def _precise_two_domain_area_polygons(domains, boundaries, width, height):
    """Close one measured open seam along both possible frame routes.

    This makes the numerical area and translucent overlay use the same
    precise TEM-refined boundary, rather than the earlier local-FFT cell
    staircase.  More complex multi-domain topologies retain the conservative
    cell-based fallback.
    """
    if len(domains) != 2 or len(boundaries) != 1:
        return False
    chains = _ordered_boundary_chains(boundaries[0].get("segments") or [])
    if not chains:
        return False
    chain = max(chains, key=lambda points: sum(math.hypot(
        second[0]-first[0], second[1]-first[1])
        for first, second in zip(points[:-1], points[1:])))
    if len(chain) < 3 or chain[0] == chain[-1]:
        return False
    frame_width = max(1.0, float(width)-1.0)
    frame_height = max(1.0, float(height)-1.0)
    perimeter = 2.0*(frame_width+frame_height)

    def snap_to_frame(point):
        x, y = float(point[0]), float(point[1])
        candidates = [
            (y, (float(np.clip(x, 0, frame_width)), 0.0),
             float(np.clip(x, 0, frame_width))),
            (frame_width-x, (frame_width,
             float(np.clip(y, 0, frame_height))),
             frame_width+float(np.clip(y, 0, frame_height))),
            (frame_height-y, (float(np.clip(x, 0, frame_width)),
             frame_height), frame_width+frame_height+
             (frame_width-float(np.clip(x, 0, frame_width)))),
            (x, (0.0, float(np.clip(y, 0, frame_height))),
             2.0*frame_width+frame_height+
             (frame_height-float(np.clip(y, 0, frame_height)))),
        ]
        unused_distance, snapped, coordinate = min(
            candidates, key=lambda item: abs(item[0]))
        return snapped, coordinate % perimeter

    def clockwise_route(first, first_t, second, second_t):
        target_t = second_t
        if target_t <= first_t+1e-9:
            target_t += perimeter
        corners = [
            (0.0, (0.0, 0.0)),
            (frame_width, (frame_width, 0.0)),
            (frame_width+frame_height, (frame_width, frame_height)),
            (2.0*frame_width+frame_height, (0.0, frame_height)),
            (perimeter, (0.0, 0.0)),
        ]
        points = [first]
        for offset in (0.0, perimeter):
            for corner_t, corner in corners:
                value = corner_t+offset
                if first_t+1e-9 < value < target_t-1e-9:
                    points.append(corner)
        points.append(second)
        return points

    start, start_t = snap_to_frame(chain[0])
    end, end_t = snap_to_frame(chain[-1])
    chain = [start]+chain[1:-1]+[end]
    clockwise = clockwise_route(end, end_t, start, start_t)
    reverse_clockwise = list(reversed(
        clockwise_route(start, start_t, end, end_t)))
    first_polygon = chain+clockwise[1:]
    second_polygon = chain+reverse_clockwise[1:]
    polygons = [first_polygon, second_polygon]
    assignments = []
    for domain in domains:
        marker = (float(domain.get("marker_x", width/2.0)),
                  float(domain.get("marker_y", height/2.0)))
        containing = [index for index, polygon in enumerate(polygons)
                      if _point_in_polygon(marker, polygon)]
        assignments.append(containing[0] if len(containing) == 1 else -1)
    if sorted(assignments) != [0, 1]:
        return False
    areas = [_polygon_area(polygons[index]) for index in assignments]
    area_sum = sum(areas)
    if area_sum <= 1e-9:
        return False
    for domain, polygon_index, area in zip(domains, assignments, areas):
        domain["area_polygons"] = [[list(point)
                                    for point in polygons[polygon_index]]]
        domain["area_fraction"] = float(area/area_sum)
        domain["area_overlay_basis"] = "precise_tem_domain_boundary"
    return True


def _polygon_mask_at_analysis_resolution(polygon, width, height,
                                         maximum_side=760):
    """Rasterize a vector polygon only for stable area integration.

    Rendering remains fully vectorial.  This bounded-resolution mask is used
    solely to normalize domain percentages inside an arbitrary curved sample
    boundary without requiring an optional geometry package.
    """
    if len(polygon) < 3:
        return None
    step = max(1.0, max(float(width), float(height))/float(maximum_side))
    xs = np.arange(step*.5, float(width), step, dtype=float)
    ys = np.arange(step*.5, float(height), step, dtype=float)
    mask = np.zeros((len(ys), len(xs)), dtype=bool)
    edges = list(zip(polygon, polygon[1:]+polygon[:1]))
    for row, y in enumerate(ys):
        intersections = []
        for first, second in edges:
            x1, y1 = float(first[0]), float(first[1])
            x2, y2 = float(second[0]), float(second[1])
            if (y1 > y) == (y2 > y):
                continue
            intersections.append(x1+(y-y1)*(x2-x1)/(y2-y1))
        intersections.sort()
        for start, end in zip(intersections[0::2], intersections[1::2]):
            mask[row] ^= ((xs >= min(start, end)) &
                          (xs < max(start, end)))
    return mask


def _normalize_domain_areas_to_sample_boundary(domains, sample_boundary,
                                               width, height):
    """Normalize spatial-domain areas to the exact analyzed sample area."""
    if not domains or not (sample_boundary or {}).get("valid"):
        return False
    sample_chains = _ordered_boundary_chains(
        sample_boundary.get("segments") or [])
    if not sample_chains:
        return False
    sample_polygon = max(sample_chains, key=lambda points: _polygon_area(
        points) if len(points) >= 3 else 0.0)
    sample_mask = _polygon_mask_at_analysis_resolution(
        sample_polygon, width, height)
    if sample_mask is None or not np.any(sample_mask):
        return False
    if len(domains) == 1:
        domains[0]["area_fraction"] = 1.0
        domains[0]["area_normalization_basis"] = (
            "precise_sample_boundary_100_percent")
        return True
    domain_masks = []
    for domain in domains:
        polygons = domain.get("area_polygons") or []
        if not polygons:
            return False
        combined = np.zeros(sample_mask.shape, dtype=bool)
        for polygon in polygons:
            current = _polygon_mask_at_analysis_resolution(
                [(float(point[0]), float(point[1])) for point in polygon],
                width, height)
            if current is not None:
                combined |= current
        domain_masks.append(combined & sample_mask)
    areas = [int(np.count_nonzero(mask)) for mask in domain_masks]
    total = sum(areas)
    if total <= 0:
        return False
    for domain, area in zip(domains, areas):
        domain["area_fraction"] = float(area)/float(total)
        domain["area_normalization_basis"] = (
            "precise_domain_polygons_within_sample_boundary")
    return True


def _reduce_boundary_staircase(segments, tolerance):
    """Remove local-FFT grid staircases while retaining measured bends.

    TEM domain contours arrive as connected vector segments.  Douglas-Peucker
    reduction is applied independently to each open chain, followed by two
    conservative corner-cut passes.  The tolerance is tied to the local-FFT
    stride/search scale, so a visible large bend survives while dozens of
    horizontal/vertical micro-steps become one oblique or gently curved run.
    """
    if not segments:
        return segments
    adjacency, edges = {}, set()
    for x1, y1, x2, y2 in segments:
        first = (round(float(x1), 4), round(float(y1), 4))
        second = (round(float(x2), 4), round(float(y2), 4))
        if first == second:
            continue
        edge = tuple(sorted((first, second)))
        edges.add(edge)
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)

    chains, unused = [], set(edges)
    while unused:
        degree = {}
        for edge in unused:
            for point in edge:
                degree[point] = degree.get(point, 0)+1
        endpoints = [point for point, value in degree.items() if value == 1]
        current = endpoints[0] if endpoints else next(iter(unused))[0]
        chain, previous = [current], None
        while True:
            options = [neighbor for neighbor in adjacency.get(current, ())
                       if tuple(sorted((current, neighbor))) in unused]
            if not options:
                break
            if previous is not None and len(options) > 1:
                old_angle = math.atan2(current[1]-previous[1],
                                       current[0]-previous[0])
                options.sort(key=lambda point: abs(math.atan2(
                    math.sin(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle),
                    math.cos(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle))))
            following = options[0]
            unused.remove(tuple(sorted((current, following))))
            chain.append(following)
            previous, current = current, following
            if current == chain[0]:
                break
        if len(chain) >= 2:
            chains.append(chain)

    def reduce_points(points):
        if len(points) < 3:
            return points
        first = np.asarray(points[0], dtype=float)
        last = np.asarray(points[-1], dtype=float)
        vector = last-first
        norm = float(np.hypot(vector[0], vector[1]))
        values = np.asarray(points, dtype=float)
        if norm <= 1e-9:
            distances = np.hypot(values[:, 0]-first[0],
                                 values[:, 1]-first[1])
        else:
            distances = np.abs(
                vector[0]*(first[1]-values[:, 1])-
                (first[0]-values[:, 0])*vector[1])/norm
        index = int(np.argmax(distances))
        if float(distances[index]) <= float(tolerance):
            return [points[0], points[-1]]
        return (reduce_points(points[:index+1])[:-1]+
                reduce_points(points[index:]))

    reduced_segments = []
    for chain in chains:
        closed = len(chain) > 2 and chain[0] == chain[-1]
        if closed:
            # A closed contour has coincident RDP endpoints.  Preserve it here
            # rather than accidentally collapsing it; the preceding evidence
            # refinement and corner-cut smoothing have already regularized it.
            points = chain
        else:
            points = reduce_points(chain)
        reduced_segments.extend([
            [float(first[0]), float(first[1]),
             float(second[0]), float(second[1])]
            for first, second in zip(points[:-1], points[1:])])
    if not reduced_segments:
        return segments
    rounded, unused_label = _smooth_boundary_segments(
        reduced_segments, iterations=2)
    return rounded or reduced_segments


def _tem_texture_boundary_evidence(image):
    """Return multi-scale TEM evidence for a visible grain/domain boundary.

    Local FFT determines which lattice family lies on each side.  The visible
    seam itself is positioned from the TEM: slow intensity edges, changes in
    local texture orientation/anisotropy, variance transitions and dark grain
    boundary ridges.  The calculation is deliberately performed at a bounded
    resolution and never invents a spline-smoothed boundary.
    """
    gray = np.asarray(image, dtype=np.float64)
    if float(np.max(gray)) > 1.5:
        gray = gray/255.0
    step = max(1, int(math.ceil(max(gray.shape)/1400.0)))
    work = gray[::step, ::step]

    def unit_feature(array):
        array = np.asarray(array, dtype=float)
        low = float(np.percentile(array, 35.0))
        high = float(np.percentile(array, 99.0))
        return np.clip((array-low)/max(high-low, 1e-9), 0.0, 1.0)

    smooth = _gaussian_filter_fft(work, 1.5)
    grad_y, grad_x = np.gradient(smooth)
    intensity_edge = np.hypot(grad_x, grad_y)

    high_pass = work-_gaussian_filter_fft(work, 2.4)
    texture_y, texture_x = np.gradient(high_pass)
    tensor_sigma = 5.0
    jxx = _gaussian_filter_fft(texture_x*texture_x, tensor_sigma)
    jyy = _gaussian_filter_fft(texture_y*texture_y, tensor_sigma)
    jxy = _gaussian_filter_fft(texture_x*texture_y, tensor_sigma)
    trace = jxx+jyy+1e-12
    orientation_x = (jxx-jyy)/trace
    orientation_y = 2.0*jxy/trace
    ox_y, ox_x = np.gradient(orientation_x)
    oy_y, oy_x = np.gradient(orientation_y)
    orientation_edge = np.sqrt(
        ox_x*ox_x+ox_y*ox_y+oy_x*oy_x+oy_y*oy_y)

    local_variance = _gaussian_filter_fft(high_pass*high_pass, 5.0)
    variance_y, variance_x = np.gradient(local_variance)
    variance_edge = np.hypot(variance_x, variance_y)
    dark_ridge = np.maximum(_gaussian_filter_fft(work, 6.0)-smooth, 0.0)
    evidence = (
        .42*unit_feature(orientation_edge)+
        .24*unit_feature(intensity_edge)+
        .20*unit_feature(variance_edge)+
        .14*unit_feature(dark_ridge))
    return evidence, float(step)


def _tem_dark_seam_boundary_evidence(image):
    """Return TEM evidence centred on a visibly dark twin/grain seam.

    An orientation-gradient map has two maxima, one on either shoulder of a
    dark physical seam.  Snapping to that map alone can therefore draw a line
    parallel to, but several pixels away from, the boundary.  For same-a,
    same-symmetry twins we combine the texture evidence with a dark
    difference-of-Gaussians ridge so the reported vector follows the seam
    centre rather than a local-FFT cell centre or one of its shoulders.
    """
    texture, step = _tem_texture_boundary_evidence(image)
    gray = np.asarray(image, dtype=np.float64)
    if float(np.max(gray)) > 1.5:
        gray = gray/255.0
    stride = max(1, int(round(float(step))))
    work = gray[::stride, ::stride]
    fine = _gaussian_filter_fft(work, 1.15)
    broad = _gaussian_filter_fft(work, 7.0)
    dark_ridge = np.maximum(broad-fine, 0.0)
    low = float(np.percentile(dark_ridge, 38.0))
    high = float(np.percentile(dark_ridge, 99.4))
    dark_ridge = np.clip(
        (dark_ridge-low)/max(high-low, 1e-9), 0.0, 1.0)
    return .18*np.asarray(texture, dtype=float)+.82*dark_ridge, float(step)


def _trace_monotone_tem_seam(segments, evidence, evidence_scale,
                             search_pixels=72):
    """Trace a mostly vertical/horizontal seam at near-pixel resolution.

    A coarse local-FFT boundary is a staircase of cell edges.  Offsetting
    nodes only along each edge normal cannot turn that staircase into the
    oblique physical seam.  When the observed boundary is predominantly
    monotone, use every TEM row (or column) as a free path state and solve a
    continuity-regularized dynamic program through the dark-seam evidence.
    """
    if not segments or evidence is None or not np.size(evidence):
        return None
    scale = max(float(evidence_scale), 1e-9)
    points = []
    for x1, y1, x2, y2 in segments:
        points.extend(((float(x1)/scale, float(y1)/scale),
                       (float(x2)/scale, float(y2)/scale)))
    points = np.asarray(points, dtype=float)
    if len(points) < 4:
        return None
    span_x = float(np.ptp(points[:, 0]))
    span_y = float(np.ptp(points[:, 1]))
    if max(span_x, span_y) < 12.0 or max(span_x, span_y) < 1.18*min(
            span_x, span_y):
        return None
    primary_axis = 1 if span_y >= span_x else 0
    secondary_axis = 1-primary_axis
    primary = points[:, primary_axis]
    secondary = points[:, secondary_axis]
    rounded_primary = np.rint(primary).astype(int)
    unique_primary = np.unique(rounded_primary)
    if len(unique_primary) < 4:
        return None
    secondary_centres = np.asarray([
        np.median(secondary[rounded_primary == value])
        for value in unique_primary], dtype=float)
    height, width = evidence.shape
    primary_limit = height if primary_axis == 1 else width
    secondary_limit = width if primary_axis == 1 else height
    primary_values = np.arange(
        max(1, int(unique_primary.min())),
        min(primary_limit-2, int(unique_primary.max()))+1, dtype=int)
    if len(primary_values) < 8:
        return None
    centres = np.interp(primary_values, unique_primary, secondary_centres)
    search = max(8, min(int(round(float(search_pixels)/scale)), 88))
    offsets = np.arange(-search, search+1, dtype=int)
    local_scores = np.full((len(primary_values), len(offsets)),
                           -1e9, dtype=float)
    for index, (primary_value, centre) in enumerate(zip(primary_values,
                                                         centres)):
        secondary_values = np.rint(centre+offsets).astype(int)
        valid = ((secondary_values >= 1) &
                 (secondary_values < secondary_limit-1))
        if primary_axis == 1:
            local_scores[index, valid] = evidence[
                primary_value, secondary_values[valid]]
        else:
            local_scores[index, valid] = evidence[
                secondary_values[valid], primary_value]
        local_scores[index, valid] -= .0008*np.abs(offsets[valid])
    objective = local_scores[0].copy()
    parents = np.zeros(local_scores.shape, dtype=np.int16)
    state_indices = np.arange(len(offsets))
    delta = np.abs(state_indices[:, None]-state_indices[None, :])
    transition_cost = .10*delta
    for index in range(1, len(primary_values)):
        transitions = objective[:, None]-transition_cost
        parents[index] = np.argmax(transitions, axis=0)
        objective = local_scores[index]+np.max(transitions, axis=0)
    choice = int(np.argmax(objective))
    choices = [choice]
    for index in range(len(primary_values)-1, 0, -1):
        choices.append(int(parents[index, choices[-1]]))
    choices.reverse()
    secondary_values = centres+offsets[np.asarray(choices, dtype=int)]
    # Suppress single-row lattice noise without moving the broad physical
    # path away from the dark seam.
    for unused in range(5):
        old = secondary_values.copy()
        secondary_values[1:-1] = (
            .20*old[:-2]+.60*old[1:-1]+.20*old[2:])
    if primary_axis == 1:
        traced = np.column_stack((secondary_values*scale,
                                  primary_values*scale))
    else:
        traced = np.column_stack((primary_values*scale,
                                  secondary_values*scale))
    # Three working pixels per SVG segment retain the measured bends while
    # avoiding thousands of redundant vector nodes on high-resolution TEMs.
    traced = traced[::3]
    return [[float(first[0]), float(first[1]),
             float(second[0]), float(second[1])]
            for first, second in zip(traced[:-1], traced[1:])]


def _tem_sample_boundary_evidence(image):
    """Return real-space evidence specialized for specimen/background edges.

    Domain seams are dominated by lattice-orientation changes, whereas a
    specimen outline is usually expressed by a broad mass-thickness/tone edge
    plus a change in periodic texture.  Keeping these evidence maps separate
    prevents the sample contour from snapping to a strong internal lattice
    transition, as happened in high-resolution cryo/TEM fields.
    """
    gray = np.asarray(image, dtype=np.float64)
    if float(np.max(gray)) > 1.5:
        gray = gray/255.0
    step = max(1, int(math.ceil(max(gray.shape)/1400.0)))
    work = gray[::step, ::step]

    def unit_feature(array):
        array = np.asarray(array, dtype=float)
        low = float(np.percentile(array, 30.0))
        high = float(np.percentile(array, 99.2))
        return np.clip((array-low)/max(high-low, 1e-9), 0.0, 1.0)

    fine = _gaussian_filter_fft(work, 2.2)
    broad = _gaussian_filter_fft(work, 8.0)
    fine_y, fine_x = np.gradient(fine)
    broad_y, broad_x = np.gradient(broad)
    fine_edge = np.hypot(fine_x, fine_y)
    broad_edge = np.hypot(broad_x, broad_y)
    high_pass = work-_gaussian_filter_fft(work, 2.0)
    local_texture = _gaussian_filter_fft(high_pass*high_pass, 7.0)
    texture_y, texture_x = np.gradient(local_texture)
    texture_edge = np.hypot(texture_x, texture_y)
    evidence = (
        .48*unit_feature(broad_edge)+
        .30*unit_feature(fine_edge)+
        .22*unit_feature(texture_edge))
    return evidence, float(step)


def _tem_precise_sample_boundary_evidence(image):
    """Combine specimen contrast with the proven domain-seam evidence.

    Sample/background classification and domain classification remain
    independent.  Only the final sub-pixel/vector placement is shared: broad
    mass-thickness contrast keeps the path on the specimen rim, while the
    texture-orientation and dark-ridge terms make it hug the visible TEM edge
    as accurately as a domain/domain boundary.
    """
    sample, step = _tem_sample_boundary_evidence(image)
    texture, texture_step = _tem_texture_boundary_evidence(image)
    dark, dark_step = _tem_dark_seam_boundary_evidence(image)
    if (sample.shape != texture.shape or sample.shape != dark.shape or
            abs(float(step)-float(texture_step)) > 1e-6 or
            abs(float(step)-float(dark_step)) > 1e-6):
        return sample, float(step)
    evidence = (.50*np.asarray(sample, dtype=float)+
                .35*np.asarray(texture, dtype=float)+
                .15*np.asarray(dark, dtype=float))
    return evidence, float(step)


def _radial_sample_boundary(image, scale_bar=None):
    """Trace a closed specimen/background envelope from the raw TEM.

    A specimen outline and an internal crystallographic boundary are not the
    same segmentation problem.  The former is normally a single, large,
    centre-containing envelope.  Searching a continuous radius around the
    full 360 degrees prevents a coarse local-FFT mask from folding back onto
    internal Moire cells or producing self-intersecting loops.  A three-turn
    Viterbi pass enforces angular continuity while still retaining corners and
    curved physical edges.

    Full-field specimens have no image-visible exterior.  They are rejected
    by the normalized-radius gate and therefore do not acquire an artificial
    green contour (the ``twin2.tif`` regression case).
    """
    evidence, evidence_scale = _tem_precise_sample_boundary_evidence(image)
    if evidence is None or not np.size(evidence):
        return {"valid": False, "segments": []}
    evidence = np.asarray(evidence, dtype=float).copy()
    height, width = evidence.shape
    center_x, center_y = (width-1.0)/2.0, (height-1.0)/2.0

    # The scale bar/text is an annotation, not a specimen edge.  Suppress only
    # its measured rectangle (plus a small antialiasing margin), rather than a
    # broad lower-image band that could erase a genuine sample boundary.
    if scale_bar:
        margin = max(24.0, .04*min(image.shape))
        x0 = max(0, int(math.floor((float(scale_bar.get("x0", 0))-margin) /
                                   evidence_scale)))
        x1 = min(width, int(math.ceil((float(scale_bar.get("x1", 0))+margin) /
                                      evidence_scale))+1)
        y0 = max(0, int(math.floor((float(scale_bar.get("y0", 0))-margin) /
                                   evidence_scale)))
        y1 = min(height, int(math.ceil((float(scale_bar.get("y1", 0))+margin) /
                                      evidence_scale))+1)
        if x1 > x0 and y1 > y0:
            evidence[y0:y1, x0:x1] = -1.0

    angles = np.linspace(0.0, 2.0*math.pi, 240, endpoint=False)
    radii_fraction = np.linspace(.26, .985, 100)
    local = np.zeros((len(angles), len(radii_fraction)), dtype=float)
    coordinates = np.zeros((len(angles), len(radii_fraction), 2),
                           dtype=float)
    maximum_radii = np.zeros(len(angles), dtype=float)
    for angle_index, angle in enumerate(angles):
        unit_x, unit_y = math.cos(angle), math.sin(angle)
        limits = []
        if unit_x > 1e-9:
            limits.append((width-1.0-center_x)/unit_x)
        elif unit_x < -1e-9:
            limits.append(-center_x/unit_x)
        if unit_y > 1e-9:
            limits.append((height-1.0-center_y)/unit_y)
        elif unit_y < -1e-9:
            limits.append(-center_y/unit_y)
        maximum_radius = min(value for value in limits if value > 0.0)
        maximum_radii[angle_index] = maximum_radius
        radii = radii_fraction*maximum_radius
        columns = np.clip(np.rint(center_x+radii*unit_x).astype(int),
                          0, width-1)
        rows = np.clip(np.rint(center_y+radii*unit_y).astype(int),
                       0, height-1)
        coordinates[angle_index, :, 0] = center_x+radii*unit_x
        coordinates[angle_index, :, 1] = center_y+radii*unit_y
        # A very small outward preference breaks ties between an internal
        # periodic edge and the physical outer rim; TEM evidence remains the
        # dominant term.
        local[angle_index] = evidence[rows, columns]+.025*radii_fraction

    # Repeating the circle lets the state at 360 degrees feed back into the
    # state at 0 degrees without an expensive all-start-state enumeration.
    repeated = np.vstack([local, local, local])
    parents = np.zeros(repeated.shape, dtype=np.int16)
    objective = repeated[0].copy()
    transition_cost = 2.2*np.abs(
        radii_fraction[:, None]-radii_fraction[None, :])
    for angle_index in range(1, len(repeated)):
        transitions = objective[:, None]-transition_cost
        parents[angle_index] = np.argmax(transitions, axis=0)
        objective = repeated[angle_index]+np.max(transitions, axis=0)
    choice = int(np.argmax(objective))
    choices = []
    for angle_index in range(len(repeated)-1, -1, -1):
        choices.append(choice)
        if angle_index:
            choice = int(parents[angle_index, choice])
    choices = choices[::-1][-len(angles):]
    normalized_radii = radii_fraction[np.asarray(choices, dtype=int)]
    # The twin regression is a full-field sample: its spurious radial path is
    # internal (median normalized radius about 0.35).  Real exterior outlines
    # in the two supplied specimen examples lie beyond 0.8.  Keep a margin so
    # moderately cropped specimens remain supported.
    if (float(np.median(normalized_radii)) < .62 or
            float(np.percentile(normalized_radii, 10.0)) < .42):
        return {"valid": False, "segments": [],
                "reason": "no_closed_outer_specimen_envelope"}

    points = np.asarray([
        coordinates[index, choices[index]]
        for index in range(len(angles))], dtype=float)
    # Low-pass only along the angular direction.  This removes one-pixel FFT
    # noise without imposing a polygon/ellipse or straightening true corners.
    for unused in range(5):
        points = (.10*np.roll(points, 2, axis=0)+
                  .20*np.roll(points, 1, axis=0)+.40*points+
                  .20*np.roll(points, -1, axis=0)+
                  .10*np.roll(points, -2, axis=0))
    points *= float(evidence_scale)
    points[:, 0] = np.clip(points[:, 0], 0.0, image.shape[1]-1.0)
    points[:, 1] = np.clip(points[:, 1], 0.0, image.shape[0]-1.0)
    points = np.vstack([points, points[0]])
    segments = [[float(first[0]), float(first[1]),
                 float(second[0]), float(second[1])]
                for first, second in zip(points[:-1], points[1:])]
    # Radial/Viterbi segmentation determines which side is specimen.  Reuse
    # the domain-boundary texture evidence only for the final placement so the
    # green contour follows the visible physical rim instead of a smoothed
    # radial approximation.
    segments = _refine_boundary_segments_from_tem(
        segments, evidence, evidence_scale,
        search_px=max(18.0, min(image.shape)*.035), maximum_search=56,
        offset_penalty=.0025, transition_penalty=.055,
        smoothing_iterations=1, reduction_tolerance=.06)
    label_segment = segments[len(segments)//2] if segments else None
    return {
        "valid": bool(segments),
        "segments": segments,
        "label_x": ((label_segment[0]+label_segment[2])/2.0
                    if label_segment else None),
        "label_y": ((label_segment[1]+label_segment[3])/2.0
                    if label_segment else None),
        "boundary_method": "continuous_360_tem_specimen_envelope",
        "median_normalized_radius": float(np.median(normalized_radii)),
    }


def _periodicity_sample_boundary(image, scale_bar=None):
    """Find an open sample/background edge from local FFT periodicity.

    The radial specimen detector is intentionally conservative and only
    accepts a centre-containing closed envelope.  Real TEM crops frequently
    contain a crystalline sheet that continues through two or three image
    edges, leaving only one open specimen/background interface.  In that
    geometry a closed radial contour is impossible.

    This fallback performs a Photoshop-magic-wand-like classification on a
    coarse grid after local background flattening: crystalline sample tiles
    have concentrated reciprocal-space peaks, whereas amorphous support or
    empty background has a diffuse FFT.  It is enabled only when the two
    periodicity classes are strongly separated.  Full-field twins and ordinary
    orientation domains therefore remain domains rather than being mislabeled
    as sample/background.
    """
    gray = np.asarray(image, dtype=np.float64)
    if float(np.max(gray)) > 1.5:
        gray = gray / 255.0
    height, width = gray.shape
    longer_tiles = 22
    if width >= height:
        columns = longer_tiles
        rows = max(12, int(round(longer_tiles * height / max(width, 1))))
    else:
        rows = longer_tiles
        columns = max(12, int(round(longer_tiles * width / max(height, 1))))
    x_edges = np.rint(np.linspace(0, width, columns + 1)).astype(int)
    y_edges = np.rint(np.linspace(0, height, rows + 1)).astype(int)
    periodicity = np.zeros((rows, columns), dtype=float)
    local_mean = np.zeros((rows, columns), dtype=float)
    local_texture = np.zeros((rows, columns), dtype=float)
    for row in range(rows):
        for column in range(columns):
            tile = gray[y_edges[row]:y_edges[row + 1],
                        x_edges[column]:x_edges[column + 1]].copy()
            if min(tile.shape) < 12:
                continue
            local_mean[row, column] = float(np.mean(tile))
            local_texture[row, column] = float(
                np.mean(np.abs(np.diff(tile, axis=0))) +
                np.mean(np.abs(np.diff(tile, axis=1))))
            # Remove local tone/slow illumination before measuring reciprocal
            # peak concentration.  The displayed TEM itself is never altered.
            tile -= local_mean[row, column]
            window = np.outer(np.hanning(tile.shape[0]),
                              np.hanning(tile.shape[1]))
            spectrum = np.abs(np.fft.fftshift(np.fft.fft2(tile * window)))
            yy, xx = np.mgrid[:tile.shape[0], :tile.shape[1]]
            radius = np.hypot(xx - tile.shape[1] / 2.0,
                              yy - tile.shape[0] / 2.0)
            annulus = ((radius > 3.0) &
                       (radius < min(tile.shape) * .42))
            values = spectrum[annulus]
            if not values.size:
                continue
            periodicity[row, column] = float(
                np.percentile(values, 99.5) /
                max(float(np.median(values)), 1e-9))

    values = periodicity[np.isfinite(periodicity) & (periodicity > 0.0)]
    if values.size < max(20, rows * columns // 3):
        return {"valid": False, "segments": [],
                "reason": "insufficient_local_periodicity_samples"}
    # Three classes are necessary: a crystalline sheet can contain both a
    # clean, strongly periodic interior and a darker/contaminated but still
    # genuine sample rim.  A two-class split would call that rim background
    # and draw the boundary through the middle of the specimen.
    features = np.column_stack((
        local_mean.ravel(), local_texture.ravel(),
        np.log(np.maximum(periodicity.ravel(), 1e-9))))
    feature_median = np.median(features, axis=0)
    feature_scale = (1.4826*np.median(
        np.abs(features-feature_median), axis=0)+1e-8)
    normalized = (features-feature_median)/feature_scale
    weakest = int(np.argmin(features[:, 2]))
    strongest = int(np.argmax(features[:, 2]))
    weak_distance = np.sum((normalized-normalized[weakest])**2, axis=1)
    strong_distance = np.sum((normalized-normalized[strongest])**2, axis=1)
    third = int(np.argmax(np.minimum(weak_distance, strong_distance)))
    centres = normalized[[weakest, strongest, third]].copy()
    flat_labels = np.zeros(len(features), dtype=np.int8)
    for unused in range(48):
        distances = np.sum((normalized[:, None, :]-
                            centres[None, :, :])**2, axis=2)
        flat_labels = np.argmin(distances, axis=1).astype(np.int8)
        updated = np.asarray([
            np.mean(normalized[flat_labels == index], axis=0)
            if np.any(flat_labels == index) else centres[index]
            for index in range(3)], dtype=float)
        if float(np.max(np.abs(updated-centres))) < 1e-6:
            centres = updated
            break
        centres = updated
    labels = flat_labels.reshape(rows, columns)
    cluster_periodicity = np.asarray([
        float(np.exp(np.mean(features[flat_labels == index, 2])))
        if np.any(flat_labels == index) else float("inf")
        for index in range(3)], dtype=float)
    order = np.argsort(cluster_periodicity)
    background_index = int(order[0])
    next_index, sample_index = int(order[1]), int(order[2])
    separation_ratio = float(
        cluster_periodicity[sample_index] /
        max(cluster_periodicity[background_index], 1e-9))
    rim_separation_ratio = float(
        cluster_periodicity[next_index] /
        max(cluster_periodicity[background_index], 1e-9))
    background = labels == background_index
    background_fraction = float(np.mean(background))
    # This high-specificity gate is what keeps full-field lattice twins from
    # becoming false specimen boundaries.  The supplied open-edge Kagome TEM
    # has a ratio around 2.7; the full-field twin controls are around 1.5.
    if (separation_ratio < 1.90 or rim_separation_ratio < 1.55 or
            background_fraction < .07 or
            background_fraction > .72):
        return {"valid": False, "segments": [],
                "reason": "periodicity_classes_not_sample_background",
                "periodicity_separation_ratio": separation_ratio,
                "rim_periodicity_separation_ratio": rim_separation_ratio}

    # Keep only low-periodicity regions connected to an image edge.  Enclosed
    # lattice defects and pores are sample features, not exterior background.
    border = np.zeros_like(background, dtype=bool)
    border[0, :], border[-1, :] = True, True
    border[:, 0], border[:, -1] = True, True
    exterior = np.zeros_like(background, dtype=bool)
    stack = [tuple(value) for value in np.argwhere(background & border)]
    for row, column in stack:
        exterior[row, column] = True
    while stack:
        row, column = stack.pop()
        for next_row, next_column in ((row - 1, column),
                                      (row + 1, column),
                                      (row, column - 1),
                                      (row, column + 1)):
            if (next_row < 0 or next_column < 0 or next_row >= rows or
                    next_column >= columns or
                    exterior[next_row, next_column] or
                    not background[next_row, next_column]):
                continue
            exterior[next_row, next_column] = True
            stack.append((next_row, next_column))
    exterior_fraction = float(np.mean(exterior))
    if exterior_fraction < .06:
        return {"valid": False, "segments": [],
                "reason": "no_coherent_border_background"}

    # Two conservative majority passes remove isolated grid tiles while
    # retaining bends in the physical edge.
    for unused in range(2):
        padded = np.pad(exterior.astype(np.int8), 1, mode="edge")
        votes = np.zeros_like(exterior, dtype=np.int16)
        for dy in range(3):
            for dx in range(3):
                votes += padded[dy:dy + rows, dx:dx + columns]
        exterior = votes >= 5

    segments = []
    for row in range(rows):
        y0, y1 = float(y_edges[row]), float(y_edges[row + 1])
        for column in range(columns - 1):
            if bool(exterior[row, column]) == bool(exterior[row, column + 1]):
                continue
            x = float(x_edges[column + 1])
            segments.append([x, y0, x, y1])
    for row in range(rows - 1):
        y = float(y_edges[row + 1])
        for column in range(columns):
            if bool(exterior[row, column]) == bool(exterior[row + 1, column]):
                continue
            x0, x1 = float(x_edges[column]), float(x_edges[column + 1])
            segments.append([x0, y, x1, y])
    if not segments:
        return {"valid": False, "segments": [],
                "reason": "no_open_sample_interface"}

    evidence, evidence_scale = _tem_precise_sample_boundary_evidence(gray)
    cell_size = max(width / max(columns, 1), height / max(rows, 1))
    segments = _refine_boundary_segments_from_tem(
        segments, evidence, evidence_scale,
        search_px=cell_size * .95, maximum_search=96,
        smoothing_iterations=2, reduction_tolerance=.18)
    segments, label_point = _smooth_boundary_segments(segments, iterations=1)
    label_segment = (segments[len(segments) // 2] if segments else None)
    return {
        "valid": bool(segments),
        "segments": segments,
        "label_x": (float(label_point[0]) if label_point else
                    ((label_segment[0] + label_segment[2]) / 2.0
                     if label_segment else None)),
        "label_y": (float(label_point[1]) if label_point else
                    ((label_segment[1] + label_segment[3]) / 2.0
                     if label_segment else None)),
        "boundary_method": "open_local_fft_periodicity_sample_boundary",
        "periodicity_separation_ratio": separation_ratio,
        "rim_periodicity_separation_ratio": rim_separation_ratio,
        "background_fraction": exterior_fraction,
        "confidence": "high",
    }


def _refine_boundary_segments_from_tem(segments, evidence, evidence_scale,
                                       search_px, maximum_search=72,
                                       offset_penalty=.004,
                                       transition_penalty=.045,
                                       smoothing_iterations=2,
                                       reduction_tolerance=.22):
    """Snap coarse local-FFT boundary segments onto visible TEM seams."""
    if not segments or evidence is None or not np.size(evidence):
        return segments
    height, width = evidence.shape
    search = max(2, int(round(float(search_px)/evidence_scale)))
    search = min(search, max(2, int(maximum_search)))

    def evidence_mean(x, y):
        x, y = int(round(x)), int(round(y))
        if x < 1 or y < 1 or x >= width-1 or y >= height-1:
            return 0.0
        return float(np.mean(evidence[y-1:y+2, x-1:x+2]))

    # Join grid edges first.  Refining each horizontal/vertical edge
    # independently produces cross-shaped stubs; a continuous chain instead
    # follows the visible TEM seam with one shared vertex at every junction.
    adjacency, edges = {}, set()
    for x1, y1, x2, y2 in segments:
        first = (round(float(x1), 4), round(float(y1), 4))
        second = (round(float(x2), 4), round(float(y2), 4))
        if first == second:
            continue
        edge = tuple(sorted((first, second)))
        edges.add(edge)
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)
    chains, unused = [], set(edges)
    while unused:
        degree = {}
        for edge in unused:
            for point in edge:
                degree[point] = degree.get(point, 0)+1
        endpoints = [point for point, value in degree.items() if value == 1]
        start = endpoints[0] if endpoints else next(iter(unused))[0]
        chain, current, previous = [start], start, None
        while True:
            options = [neighbor for neighbor in adjacency.get(current, ())
                       if tuple(sorted((current, neighbor))) in unused]
            if not options:
                break
            if previous is not None and len(options) > 1:
                old_angle = math.atan2(current[1]-previous[1],
                                       current[0]-previous[0])
                options.sort(key=lambda point: abs(math.atan2(
                    math.sin(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle),
                    math.cos(math.atan2(point[1]-current[1],
                                        point[0]-current[0])-old_angle))))
            neighbor = options[0]
            unused.remove(tuple(sorted((current, neighbor))))
            chain.append(neighbor)
            previous, current = current, neighbor
            if current == start:
                break
        if len(chain) >= 2:
            chains.append(chain)

    refined = []
    offsets = np.arange(-search, search+1, dtype=float)
    for raw_chain in chains:
        closed = len(raw_chain) > 2 and raw_chain[0] == raw_chain[-1]
        raw_points = raw_chain[:-1] if closed else raw_chain
        # A local-FFT grid edge may be 50--170 source pixels long.  Optimizing
        # only its endpoints inevitably draws one long straight chord across
        # a curved or angled physical boundary.  Insert closely spaced control
        # nodes first; every node is then independently snapped to the TEM
        # boundary evidence by the globally continuous dynamic programme.
        dense_points = []
        target_step = max(8.0, float(evidence_scale)*5.0)
        segment_pairs = list(zip(raw_points[:-1], raw_points[1:]))
        if closed and raw_points:
            segment_pairs.append((raw_points[-1], raw_points[0]))
        for first, second in segment_pairs:
            first = np.asarray(first, dtype=float)
            second = np.asarray(second, dtype=float)
            length = float(np.hypot(*(second-first)))
            count = max(1, int(math.ceil(length/target_step)))
            for index in range(count):
                fraction = float(index)/float(count)
                dense_points.append(tuple(
                    first*(1.0-fraction)+second*fraction))
        if not closed and raw_points:
            dense_points.append(tuple(raw_points[-1]))
        raw_points = dense_points or raw_points
        points = np.asarray([(point[0]/evidence_scale,
                              point[1]/evidence_scale)
                             for point in raw_points], dtype=float)
        if len(points) < 2:
            continue
        candidates = np.zeros((len(points), len(offsets), 2), dtype=float)
        local_scores = np.zeros((len(points), len(offsets)), dtype=float)
        for point_index, point in enumerate(points):
            previous = points[(point_index-1) % len(points)] if (
                closed or point_index > 0) else points[point_index]
            following = points[(point_index+1) % len(points)] if (
                closed or point_index+1 < len(points)) else points[point_index]
            tangent = following-previous
            norm = float(np.hypot(tangent[0], tangent[1]))
            if norm <= 1e-9:
                tangent = np.asarray([1.0, 0.0])
                norm = 1.0
            normal = np.asarray([-tangent[1]/norm, tangent[0]/norm])
            candidates[point_index] = point[None, :]+offsets[:, None]*normal
            for offset_index, candidate in enumerate(candidates[point_index]):
                local_scores[point_index, offset_index] = (
                    evidence_mean(candidate[0], candidate[1])-
                    float(offset_penalty)*abs(offsets[offset_index]))
        objective = local_scores[0].copy()
        parents = np.zeros((len(points), len(offsets)), dtype=int)
        for point_index in range(1, len(points)):
            updated = np.full(len(offsets), -1e99, dtype=float)
            for offset_index in range(len(offsets)):
                transition = objective-float(transition_penalty)*np.abs(
                    np.arange(len(offsets))-offset_index)
                parent = int(np.argmax(transition))
                updated[offset_index] = (
                    local_scores[point_index, offset_index]+transition[parent])
                parents[point_index, offset_index] = parent
            objective = updated
        chosen = [int(np.argmax(objective))]
        for point_index in range(len(points)-1, 0, -1):
            chosen.append(int(parents[point_index, chosen[-1]]))
        chosen.reverse()
        final_array = np.asarray(
            [candidates[index, choice]*evidence_scale
             for index, choice in enumerate(chosen)], dtype=float)
        # Remove only the high-frequency rectangular-grid stair-step.  This is
        # a local low-pass over densely measured nodes, not endpoint/path
        # simplification: broad bends, oblique runs and real corners remain.
        if len(final_array) >= 5:
            for unused_pass in range(4):
                old = final_array.copy()
                if closed:
                    final_array = (
                        .25*np.roll(old, 1, axis=0)+.50*old+
                        .25*np.roll(old, -1, axis=0))
                else:
                    final_array[1:-1] = (
                        .25*old[:-2]+.50*old[1:-1]+.25*old[2:])
        final_points = [tuple(point) for point in final_array]
        if closed:
            final_points.append(final_points[0])
        refined.extend([[float(first[0]), float(first[1]),
                         float(second[0]), float(second[1])]
                        for first, second in zip(final_points[:-1],
                                                 final_points[1:])])
    if refined:
        # Convert the dense evidence-following polyline into short connected
        # corner-cut vector segments.  Two conservative Chaikin passes remove
        # the local-FFT grid staircase while preserving the measured global
        # route and genuine changes of direction; unlike the former endpoint
        # simplifier this can never collapse a boundary to one straight line.
        curved, unused_label = _smooth_boundary_segments(
            refined, iterations=max(0, int(smoothing_iterations)))
        curved = curved or refined
        return _reduce_boundary_staircase(
            curved, tolerance=max(2.0, float(search_px)*float(
                reduction_tolerance)))
    return segments


def _local_square_orientation_domains(image, families, scale_bar=None,
                                      manual_edits=None):
    """Map validated lattice families to their real-space domains.

    The historical name is retained for compatibility.  A family may contain
    the two reciprocal axes of a Square lattice or the three reciprocal axes
    of a sixfold lattice.  Global FFT establishes the allowed bases;
    overlapping local FFTs then assign symmetry, spacing and orientation to
    each spatial domain.  Connected components and boundaries are emitted in
    source-image coordinates for identical raster/SVG rendering by the GUI.
    """
    manual_edits = list(manual_edits or [])
    valid_families = []
    for index, family in enumerate(families or []):
        symmetry = str(family.get("symmetry") or "Square")
        angles = (family.get("reciprocal_axis_angles_deg") or
                  family.get("axis_angles_deg") or [])
        periods = (family.get("reciprocal_axis_periods_px") or
                   family.get("axis_periods_px") or [])
        if len(periods) < 2 and family.get("lattice_constant_px"):
            periods = [family["lattice_constant_px"]] * 2
        if len(angles) < 2 or len(periods) < 2:
            continue
        vectors = []
        for angle, period in zip(angles[:2], periods[:2]):
            period = float(period)
            if period <= 0:
                vectors = []
                break
            radians = math.radians(float(angle))
            vectors.append((math.cos(radians)/period,
                            -math.sin(radians)/period))
        if len(vectors) >= 2:
            inter_angles = [float(value) for value in
                            (family.get("inter_axis_angles_deg") or [])]
            if not inter_angles:
                measured_angles = sorted(float(value) % 180.0
                                         for value in angles)
                if symmetry == "Square" and len(measured_angles) >= 2:
                    first_gap = measured_angles[1]-measured_angles[0]
                    inter_angles = [first_gap, 180.0-first_gap]
                elif len(measured_angles) >= 3:
                    inter_angles = [
                        measured_angles[1]-measured_angles[0],
                        measured_angles[2]-measured_angles[1],
                        180.0-measured_angles[2]+measured_angles[0],
                    ]
            valid_families.append({
                "orientation_index": index,
                "orientation_deg": float(family.get(
                    "orientation_deg", float(angles[0]) % 90.0)),
                "symmetry": symmetry,
                "lattice_constant_px": float(family.get(
                    "lattice_constant_px") or np.median(periods)),
                "reciprocal_axis_angles_deg": [float(value)
                                                 for value in angles],
                "reciprocal_axis_periods_px": [float(value)
                                                 for value in periods],
                "inter_axis_angles_deg": inter_angles,
                "vectors": vectors,
            })
    if not valid_families:
        return {"valid": False,
                "error": "No Square orientation is available for local FFT analysis."}
    gray = np.asarray(image, dtype=np.float64) / 255.0
    height, width = gray.shape
    patch = int(min(224, max(128, min(height, width) * .19)))
    patch -= patch % 16
    half = patch // 2
    # A quarter-patch stride creates roughly 4,900 local FFTs on a 4k TEM.
    # That does not improve the final boundary location: the coarse domain
    # grid is subsequently snapped to the full TEM texture evidence below.
    # Use half-patch overlap for large frames, retaining the denser historical
    # grid for smaller images where it is inexpensive.
    large_mixed_family_field = bool(
        max(height, width) >= 3072 and len(valid_families) > 1 and
        len(set(family["symmetry"] for family in valid_families)) > 1)
    stride = max(32, (patch // 2 if large_mixed_family_field else
                      ((3 * patch) // 4
                       if max(height, width) >= 3072 else patch // 4)))
    xs = np.arange(half, width-half+1, stride, dtype=int)
    ys = np.arange(half, height-half+1, stride, dtype=int)
    if not len(xs) or not len(ys):
        return {"valid": False,
                "error": "The image is too small for local FFT partitioning."}
    window = np.outer(np.hanning(patch), np.hanning(patch))
    scores = np.zeros((len(ys), len(xs), len(valid_families)), dtype=float)
    quality = np.zeros((len(ys), len(xs)), dtype=float)
    tile_means = np.zeros((len(ys), len(xs)), dtype=float)
    tile_textures = np.zeros((len(ys), len(xs)), dtype=float)
    expected_radii = [patch * math.hypot(*vector)
                      for family in valid_families
                      for vector in family["vectors"]]
    annulus_min = max(4.0, min(expected_radii) * .72)
    annulus_max = max(annulus_min + 4.0, max(expected_radii) * 1.28)
    grid_y, grid_x = np.mgrid[:patch, :patch]
    center = patch / 2.0
    radial = np.hypot(grid_x-center, grid_y-center)
    polar_angle = np.arctan2(-(grid_y-center), grid_x-center)
    annulus_mask = (radial >= annulus_min) & (radial <= annulus_max)
    for row, cy in enumerate(ys):
        for column, cx in enumerate(xs):
            tile = gray[cy-half:cy+half, cx-half:cx+half].copy()
            tile_mean = float(np.mean(tile))
            tile_means[row, column] = tile_mean
            tile_textures[row, column] = float(
                np.mean(np.abs(np.diff(tile, axis=0)))+
                np.mean(np.abs(np.diff(tile, axis=1))))
            tile -= tile_mean
            spectrum = np.log1p(np.abs(np.fft.fftshift(np.fft.fft2(
                tile * window))))
            background_values = spectrum[annulus_mask]
            background = float(np.median(background_values))
            spread = float(np.median(np.abs(
                background_values-background)) + 1e-7)
            for family_index, family in enumerate(valid_families):
                values = []
                for fx, fy in family["vectors"]:
                    px, py = center + fx*patch, center + fy*patch
                    values.extend((
                        _local_fft_disk_mean(spectrum, px, py),
                        _local_fft_disk_mean(
                            spectrum, 2*center-px, 2*center-py)))
                normalized_values = ((np.asarray(values, dtype=float) -
                                      background) / spread)
                # A valid point group needs energy on every defining axis.
                # Combining mean and lower quartile prevents two strong axes
                # of a sixfold pattern from being mislabeled as Square.
                family_score = float(
                    .68*np.mean(normalized_values) +
                    .32*np.percentile(normalized_values, 25))
                reciprocal_radii = [patch*math.hypot(fx, fy)
                                    for fx, fy in family["vectors"]]
                target_radius = float(np.median(reciprocal_radii))
                ring_half_width = max(1.5, target_radius*.065)
                ring = np.abs(radial-target_radius) <= ring_half_width
                ring_weights = np.maximum(spectrum-background, 0.0)*ring
                ring_total = float(np.sum(ring_weights))
                if ring_total > 1e-9:
                    harmonic_4 = abs(np.sum(
                        ring_weights*np.exp(4j*polar_angle)))/ring_total
                    harmonic_6 = abs(np.sum(
                        ring_weights*np.exp(6j*polar_angle)))/ring_total
                    own = (harmonic_4 if family["symmetry"] == "Square"
                           else harmonic_6)
                    competing = (harmonic_6 if family["symmetry"] == "Square"
                                 else harmonic_4)
                    family_score += 5.5*float(own-competing)
                scores[row, column, family_index] = family_score
            quality[row, column] = float(np.max(scores[row, column]))

    # Exclude non-sample background before assigning domains.  Border-seeded
    # background is accepted only when it is both descriptor-consistent and
    # substantially less periodic/different in tone than the enclosed region;
    # this prevents a full-field twin sample from being mistaken for a
    # foreground island merely because its two grains have different texture.
    border = np.zeros(quality.shape, dtype=bool)
    border[0, :], border[-1, :] = True, True
    border[:, 0], border[:, -1] = True, True

    def robust_center_scale(values):
        values = np.asarray(values, dtype=float)
        center_value = float(np.median(values))
        mad = float(np.median(np.abs(values-center_value)))
        scale_value = max(1.4826*mad, float(np.std(values))*.20, 1e-8)
        return center_value, scale_value

    # Model exterior background from four corners, not the entire frame edge.
    # A large specimen may touch the top/bottom/side border; using every edge
    # cell then contaminates the background model and floods through the
    # specimen (the failure seen in ``2.tif``).
    corner_rows = max(1, int(round(quality.shape[0]*.18)))
    corner_columns = max(1, int(round(quality.shape[1]*.18)))
    corner_reference = np.zeros_like(border)
    corner_reference[:corner_rows, :corner_columns] = True
    corner_reference[:corner_rows, -corner_columns:] = True
    corner_reference[-corner_rows:, :corner_columns] = True
    corner_reference[-corner_rows:, -corner_columns:] = True
    mean_center, mean_scale = robust_center_scale(
        tile_means[corner_reference])
    texture_center, texture_scale = robust_center_scale(
        tile_textures[corner_reference])
    quality_center, quality_scale = robust_center_scale(
        quality[corner_reference])
    background_like = (
        (np.abs(tile_means-mean_center) <= 1.65*mean_scale) &
        (np.abs(tile_textures-texture_center) <= 1.80*texture_scale) &
        (quality <= quality_center+2.5*quality_scale))
    background = np.zeros_like(background_like, dtype=bool)
    stack = [tuple(value) for value in np.argwhere(
        border & corner_reference & background_like)]
    for row, column in stack:
        background[row, column] = True
    while stack:
        row, column = stack.pop()
        for next_row, next_column in ((row-1, column), (row+1, column),
                                      (row, column-1), (row, column+1)):
            if (next_row < 0 or next_column < 0 or
                    next_row >= background.shape[0] or
                    next_column >= background.shape[1] or
                    background[next_row, next_column] or
                    not background_like[next_row, next_column]):
                continue
            background[next_row, next_column] = True
            stack.append((next_row, next_column))
    proposed_sample = ~background
    proposed_fraction = float(np.mean(proposed_sample))
    accept_background = False
    # A tiny complement of a border-derived descriptor cluster is commonly a
    # second grain inside a full-field sample, not the only specimen present.
    # Require a substantial enclosed specimen before excluding the exterior.
    if .28 <= proposed_fraction <= .94 and np.any(background):
        foreground_quality = float(np.median(quality[proposed_sample]))
        background_quality = float(np.median(quality[background]))
        tone_difference = abs(float(np.median(tile_means[proposed_sample]))-
                              float(np.median(tile_means[background])))
        texture_difference = abs(float(np.median(
            tile_textures[proposed_sample]))-float(np.median(
                tile_textures[background])))
        # A lattice-orientation domain can have a different FFT score without
        # being background.  Require an actual tone/texture envelope change;
        # quality alone is never sufficient to call a specimen outline.
        accept_background = bool(
            tone_difference > 1.65*mean_scale or
            texture_difference > 1.80*texture_scale)

    # A second, direct specimen-envelope estimate is essential when the
    # specimen itself touches one or more image borders.  Select pixels whose
    # mean tone or high-frequency texture differs robustly from the four-corner
    # background model, retain only the largest coherent object, then fill its
    # enclosed holes.  Full-field twins produce only scattered/partial
    # descriptor deviations and therefore fail the minimum-area gate.
    descriptor_candidate = (
        (np.abs(tile_means-mean_center) >= 2.00*mean_scale) |
        (np.abs(tile_textures-texture_center) >= 4.00*texture_scale))
    descriptor_components = []
    descriptor_seen = np.zeros_like(descriptor_candidate, dtype=bool)
    for row, column in zip(*np.where(descriptor_candidate)):
        if descriptor_seen[row, column]:
            continue
        cells, stack = [], [(int(row), int(column))]
        descriptor_seen[row, column] = True
        while stack:
            current_row, current_column = stack.pop()
            cells.append((current_row, current_column))
            for next_row, next_column in (
                    (current_row-1, current_column),
                    (current_row+1, current_column),
                    (current_row, current_column-1),
                    (current_row, current_column+1)):
                if (next_row < 0 or next_column < 0 or
                        next_row >= descriptor_candidate.shape[0] or
                        next_column >= descriptor_candidate.shape[1] or
                        descriptor_seen[next_row, next_column] or
                        not descriptor_candidate[next_row, next_column]):
                    continue
                descriptor_seen[next_row, next_column] = True
                stack.append((next_row, next_column))
        descriptor_components.append(cells)
    largest_descriptor = max(descriptor_components, key=len, default=[])
    descriptor_sample = np.zeros_like(descriptor_candidate, dtype=bool)
    for row, column in largest_descriptor:
        descriptor_sample[row, column] = True
    descriptor_fraction = float(np.mean(descriptor_sample))
    descriptor_valid = bool(.28 <= descriptor_fraction <= .94)
    if descriptor_valid:
        # Fill only holes enclosed by the retained specimen.  Background that
        # remains connected to the grid border is preserved.
        exterior = np.zeros_like(descriptor_sample, dtype=bool)
        exterior_stack = [tuple(value) for value in np.argwhere(
            border & ~descriptor_sample)]
        for row, column in exterior_stack:
            exterior[row, column] = True
        while exterior_stack:
            row, column = exterior_stack.pop()
            for next_row, next_column in ((row-1, column), (row+1, column),
                                          (row, column-1), (row, column+1)):
                if (next_row < 0 or next_column < 0 or
                        next_row >= exterior.shape[0] or
                        next_column >= exterior.shape[1] or
                        exterior[next_row, next_column] or
                        descriptor_sample[next_row, next_column]):
                    continue
                exterior[next_row, next_column] = True
                exterior_stack.append((next_row, next_column))
        descriptor_sample |= (~descriptor_sample & ~exterior)
        descriptor_fraction = float(np.mean(descriptor_sample))
    if descriptor_valid:
        sample_mask = descriptor_sample
        accept_background = True
    elif accept_background:
        sample_mask = proposed_sample
    else:
        sample_mask = np.ones_like(background, dtype=bool)

    manual_background_requested = any(str(edit.get("type") or "").startswith(
        "background_") for edit in manual_edits)
    automatic_background_rejected = False
    sample_periodicity_validation = {
        "inside_median": None, "outside_median": None,
        "separation_z": None, "same_lattice_across_boundary": False,
    }
    if np.any(sample_mask) and np.any(~sample_mask):
        inside_quality = float(np.median(quality[sample_mask]))
        outside_quality = float(np.median(quality[~sample_mask]))
        quality_mad = float(np.median(np.abs(
            quality-np.median(quality))))
        quality_scale_global = max(1.4826*quality_mad,
                                   float(np.std(quality))*.20, 1e-8)
        separation_z = ((inside_quality-outside_quality) /
                        quality_scale_global)
        # Global and local FFT must agree that an alleged exterior really
        # lacks the retained crystal lattice.  Similar lattice evidence on
        # both sides means TEM brightness/texture drift has cut one single
        # crystal into false sample/background islands.
        same_lattice_across_boundary = bool(separation_z < .80)
        sample_periodicity_validation = {
            "inside_median": inside_quality,
            "outside_median": outside_quality,
            "separation_z": float(separation_z),
            "same_lattice_across_boundary": same_lattice_across_boundary,
        }
        if (same_lattice_across_boundary and
                not manual_background_requested):
            sample_mask = np.ones_like(sample_mask, dtype=bool)
            accept_background = False
            descriptor_valid = False
            automatic_background_rejected = True
            sample_periodicity_validation["automatic_boundary_rejected"] = True

    def edit_grid_mask(edit):
        """Map a frozen manual boundary or selected result region to the grid."""
        region_polygons = edit.get("region_polygons") or []
        region_runs = edit.get("region_runs") or []
        if region_polygons or region_runs:
            output = np.zeros(sample_mask.shape, dtype=bool)
            for grid_row, center_y in enumerate(ys):
                for grid_column, center_x in enumerate(xs):
                    if any(_point_in_polygon(
                            (float(center_x), float(center_y)), polygon)
                           for polygon in region_polygons):
                        output[grid_row, grid_column] = True
                        continue
                    output[grid_row, grid_column] = any(
                        len(run) >= 4 and
                        float(run[0]) <= float(center_x) <= float(run[2]) and
                        float(run[1]) <= float(center_y) <= float(run[3])
                        for run in region_runs)
            return output
        polygon = edit.get("polygon") or []
        if len(polygon) >= 3:
            output = np.zeros(sample_mask.shape, dtype=bool)
            for grid_row, center_y in enumerate(ys):
                for grid_column, center_x in enumerate(xs):
                    output[grid_row, grid_column] = _point_in_polygon(
                        (float(center_x), float(center_y)), polygon)
            return output
        point = edit.get("point") or []
        if len(point) < 2:
            return np.zeros(sample_mask.shape, dtype=bool)
        seed_row = int(np.argmin(np.abs(ys-float(point[1]))))
        seed_column = int(np.argmin(np.abs(xs-float(point[0]))))
        # Point mode is an automatic region selection, not a one-cell brush.
        # Grow through neighboring windows whose measured tone, texture and
        # periodicity resemble the clicked seed.  This is deliberately local
        # and cannot jump across a strong specimen/domain interface.
        feature_arrays = (tile_means, tile_textures, quality)
        feature_scales = []
        for values in feature_arrays:
            mad = float(np.median(np.abs(values-np.median(values))))
            feature_scales.append(max(1.4826*mad, float(np.std(values))*.15,
                                      1e-8))
        distance = np.zeros(sample_mask.shape, dtype=float)
        weights = (1.0, 1.0, .55)
        for values, scale_value, weight in zip(
                feature_arrays, feature_scales, weights):
            distance += weight*np.abs(
                values-values[seed_row, seed_column])/scale_value
        compatible = distance <= 4.2
        output = np.zeros(sample_mask.shape, dtype=bool)
        stack = [(seed_row, seed_column)]
        output[seed_row, seed_column] = True
        while stack:
            grid_row, grid_column = stack.pop()
            for next_row, next_column in (
                    (grid_row-1, grid_column), (grid_row+1, grid_column),
                    (grid_row, grid_column-1), (grid_row, grid_column+1)):
                if (next_row < 0 or next_column < 0 or
                        next_row >= output.shape[0] or
                        next_column >= output.shape[1] or
                        output[next_row, next_column] or
                        not compatible[next_row, next_column]):
                    continue
                output[next_row, next_column] = True
                stack.append((next_row, next_column))
        # A very small seed component is expanded conservatively so a click is
        # still visible and can be undone, while never replacing the whole
        # specimen when local descriptors are nearly uniform.
        if int(np.count_nonzero(output)) < 4:
            row0, row1 = max(0, seed_row-1), min(output.shape[0], seed_row+2)
            col0, col1 = max(0, seed_column-1), min(output.shape[1], seed_column+2)
            output[row0:row1, col0:col1] = True
        return output

    def connected_component(mask, seed_row, seed_column):
        """Return the four-connected component containing the selected cell."""
        output = np.zeros(mask.shape, dtype=bool)
        if (seed_row < 0 or seed_column < 0 or
                seed_row >= mask.shape[0] or seed_column >= mask.shape[1] or
                not bool(mask[seed_row, seed_column])):
            return output
        output[seed_row, seed_column] = True
        stack = [(seed_row, seed_column)]
        while stack:
            row, column = stack.pop()
            for next_row, next_column in (
                    (row-1, column), (row+1, column),
                    (row, column-1), (row, column+1)):
                if (next_row < 0 or next_column < 0 or
                        next_row >= mask.shape[0] or
                        next_column >= mask.shape[1] or
                        output[next_row, next_column] or
                        not bool(mask[next_row, next_column])):
                    continue
                output[next_row, next_column] = True
                stack.append((next_row, next_column))
        return output

    manual_region_ids = np.zeros(sample_mask.shape, dtype=int)
    manual_deleted_mask = np.zeros(sample_mask.shape, dtype=bool)
    manual_region_number = 0
    manual_region_polygons = {}
    manual_background_polygons = []
    # Background edits are replayed in order.  Outline mode uses the rough
    # polygon as a hard prior; point mode uses the automatic descriptor grow
    # above.  The raw TEM boundary refinement below snaps both to the image.
    for edit in manual_edits:
        operation = str(edit.get("type") or "")
        region = edit_grid_mask(edit)
        if operation == "background_add":
            sample_mask[region] = False
            accept_background = True
            polygon = edit.get("polygon") or []
            if len(polygon) >= 3:
                manual_background_polygons.append([
                    [float(point[0]), float(point[1])]
                    for point in polygon])
        elif operation == "background_delete":
            point = edit.get("point") or []
            if len(point) >= 2:
                seed_row = int(np.argmin(np.abs(ys-float(point[1]))))
                seed_column = int(np.argmin(np.abs(xs-float(point[0]))))
                region = connected_component(
                    ~sample_mask, seed_row, seed_column)
                # Remove any earlier hand-drawn background whose continuous
                # component contains this click.  Its frozen vector contour
                # must disappear together with the deleted background.
                manual_background_polygons = [
                    polygon for polygon in manual_background_polygons
                    if not _point_in_polygon(
                        (float(point[0]), float(point[1])), polygon)]
            sample_mask[region] = True
        elif operation == "domain_add":
            manual_region_number += 1
            manual_region_ids[region & sample_mask] = manual_region_number
            polygon = edit.get("polygon") or []
            if len(polygon) >= 3:
                manual_region_polygons[manual_region_number] = [
                    [float(point[0]), float(point[1])]
                    for point in polygon]
        elif operation == "domain_delete":
            # Applied after the initial FFT labels exist, so a point removes
            # exactly the clicked analyzed domain rather than a descriptor-
            # similar patch in another domain.
            pass

    if len(valid_families) == 1:
        probabilities = np.ones_like(scores)
    else:
        normalized = scores - np.max(scores, axis=2, keepdims=True)
        probabilities = np.exp(np.clip(normalized, -20.0, 20.0))
        probabilities /= np.sum(probabilities, axis=2, keepdims=True)
        # Confidence-preserving spatial regularization suppresses isolated
        # one-window flips without erasing narrow real domains.
        for unused in range(4):
            updated = probabilities * 4.0
            updated[1:] += probabilities[:-1]
            updated[:-1] += probabilities[1:]
            updated[:, 1:] += probabilities[:, :-1]
            updated[:, :-1] += probabilities[:, 1:]
            probabilities = updated / 8.0
    labels = np.argmax(probabilities, axis=2).astype(int)
    labels[~sample_mask] = -1
    for manual_region in range(1, manual_region_number+1):
        region = (manual_region_ids == manual_region) & sample_mask
        if not np.any(region):
            continue
        # One user-selected region represents one physical domain.  Determine
        # its symmetry/a/orientation from the joint FFT evidence over the
        # whole selected region, then use that family consistently while the
        # raw-TEM refinement finds the precise boundary.
        mean_family_scores = np.asarray([
            float(np.mean(scores[:, :, family_index][region]))
            for family_index in range(scores.shape[2])])
        labels[region] = int(np.argmax(mean_family_scores))
    for edit in manual_edits:
        if str(edit.get("type") or "") != "domain_delete":
            continue
        region = edit_grid_mask(edit)
        if np.any(region):
            manual_deleted_mask |= region & sample_mask
            continue
        point = edit.get("point") or []
        if len(point) < 2:
            continue
        seed_row = int(np.argmin(np.abs(ys-float(point[1]))))
        seed_column = int(np.argmin(np.abs(xs-float(point[0]))))
        seed_label = int(labels[seed_row, seed_column])
        if seed_label < 0:
            continue
        visited = np.zeros(labels.shape, dtype=bool)
        stack = [(seed_row, seed_column)]
        visited[seed_row, seed_column] = True
        while stack:
            grid_row, grid_column = stack.pop()
            manual_deleted_mask[grid_row, grid_column] = True
            for next_row, next_column in (
                    (grid_row-1, grid_column), (grid_row+1, grid_column),
                    (grid_row, grid_column-1), (grid_row, grid_column+1)):
                if (next_row < 0 or next_column < 0 or
                        next_row >= labels.shape[0] or
                        next_column >= labels.shape[1] or
                        visited[next_row, next_column] or
                        int(labels[next_row, next_column]) != seed_label or
                        int(manual_region_ids[next_row, next_column]) !=
                        int(manual_region_ids[seed_row, seed_column])):
                    continue
                visited[next_row, next_column] = True
                stack.append((next_row, next_column))
    # Deleting a domain removes the classification boundary, not the sample.
    # Reassign the selected connected component outwards-in from its adjacent
    # domains.  This deterministically merges it into the surrounding domains
    # and cannot leave the former black/unanalysed hole.
    if np.any(manual_deleted_mask):
        original_labels = labels.copy()
        pending = manual_deleted_mask.copy()
        labels[pending] = -2
        manual_region_ids[pending] = 0
        while np.any(pending):
            assignments = []
            for row, column in zip(*np.where(pending)):
                neighbours = []
                for next_row, next_column in (
                        (row-1, column), (row+1, column),
                        (row, column-1), (row, column+1)):
                    if (0 <= next_row < labels.shape[0] and
                            0 <= next_column < labels.shape[1] and
                            int(labels[next_row, next_column]) >= 0 and
                            not pending[next_row, next_column]):
                        neighbours.append(int(labels[next_row, next_column]))
                if neighbours:
                    counts = np.bincount(
                        neighbours, minlength=len(valid_families))
                    assignments.append((row, column, int(np.argmax(counts))))
            if not assignments:
                # No surrounding domain means the click selected the entire
                # specimen.  Treat that as a no-op instead of inventing a new
                # domain or deleting the sample.
                labels[pending] = original_labels[pending]
                break
            for row, column, value in assignments:
                labels[row, column] = value
                pending[row, column] = False
    # A rough outline is a spatial prior only.  Its crystallographic identity
    # remains the locally measured best-scoring FFT family.  For a point edit,
    # the region grow likewise retains per-window lattice measurements.
    for unused in range(2):
        updated = labels.copy()
        for row in range(labels.shape[0]):
            for column in range(labels.shape[1]):
                if labels[row, column] < 0:
                    continue
                y0, y1 = max(0, row-1), min(labels.shape[0], row+2)
                x0, x1 = max(0, column-1), min(labels.shape[1], column+2)
                neighbours = labels[y0:y1, x0:x1].ravel()
                neighbours = neighbours[neighbours >= 0]
                if not len(neighbours):
                    continue
                counts = np.bincount(neighbours,
                                     minlength=len(valid_families))
                if int(np.max(counts)) >= 5:
                    updated[row, column] = int(np.argmax(counts))
        labels = updated
    confidence = np.max(probabilities, axis=2)

    # Connected components turn repeated/disconnected regions of the same
    # orientation into distinct domains with distinct labels and boundaries.
    component_map = np.full(labels.shape, -1, dtype=int)
    raw_components = []
    for row in range(labels.shape[0]):
        for column in range(labels.shape[1]):
            if component_map[row, column] >= 0:
                continue
            orientation = int(labels[row, column])
            if orientation < 0:
                component_map[row, column] = -2
                continue
            manual_region = int(manual_region_ids[row, column])
            component_id = len(raw_components)
            stack, cells = [(row, column)], []
            component_map[row, column] = component_id
            while stack:
                current_row, current_column = stack.pop()
                cells.append((current_row, current_column))
                for next_row, next_column in (
                        (current_row-1, current_column),
                        (current_row+1, current_column),
                        (current_row, current_column-1),
                        (current_row, current_column+1)):
                    if (next_row < 0 or next_column < 0 or
                            next_row >= labels.shape[0] or
                            next_column >= labels.shape[1] or
                            component_map[next_row, next_column] >= 0 or
                            int(labels[next_row, next_column]) != orientation or
                            int(manual_region_ids[next_row, next_column]) !=
                            manual_region):
                        continue
                    component_map[next_row, next_column] = component_id
                    stack.append((next_row, next_column))
            raw_components.append({"orientation": orientation,
                                   "manual_region": manual_region,
                                   "manual_protected": manual_region > 0,
                                   "cells": cells})
    sample_cell_count = max(1, int(np.sum(sample_mask)))
    minimum_cells = max(2, int(round(sample_cell_count * .008)))
    retained = [component for component in raw_components
                if (len(component["cells"]) >= minimum_cells or
                    component.get("manual_protected"))]
    retained.sort(key=lambda component: len(component["cells"]), reverse=True)
    domain_map = np.full(labels.shape, -1, dtype=int)
    domains = []
    for domain_index, component in enumerate(retained, 1):
        cells = component["cells"]
        orientation_index = int(component["orientation"])
        for row, column in cells:
            domain_map[row, column] = domain_index
        # Place the marker in a strong interior window, not on the boundary.
        best_cell, best_merit = cells[0], -1e99
        for row, column in cells:
            same_neighbors = sum(
                1 for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                if (0 <= row+dy < labels.shape[0] and
                    0 <= column+dx < labels.shape[1] and
                    int(labels[row+dy, column+dx]) == orientation_index))
            center_penalty = .0005 * math.hypot(
                float(xs[column])-width/2.0, float(ys[row])-height/2.0)
            merit = (same_neighbors + .35*quality[row, column] -
                     center_penalty)
            if merit > best_merit:
                best_cell, best_merit = (row, column), merit
        rows = [cell[0] for cell in cells]
        columns = [cell[1] for cell in cells]
        marker_row, marker_column = best_cell
        domains.append({
            "domain_id": domain_index,
            "orientation_index": orientation_index,
            "orientation_deg": valid_families[orientation_index][
                "orientation_deg"],
            "symmetry": valid_families[orientation_index]["symmetry"],
            "lattice_constant_px": valid_families[orientation_index][
                "lattice_constant_px"],
            "reciprocal_axis_angles_deg": valid_families[
                orientation_index]["reciprocal_axis_angles_deg"],
            "reciprocal_axis_periods_px": valid_families[
                orientation_index]["reciprocal_axis_periods_px"],
            "inter_axis_angles_deg": valid_families[
                orientation_index]["inter_axis_angles_deg"],
            "cell_count": len(cells),
            "area_fraction": float(len(cells)/sample_cell_count),
            "marker_x": float(xs[marker_column]),
            "marker_y": float(ys[marker_row]),
            "bbox": [float(xs[min(columns)]), float(ys[min(rows)]),
                     float(xs[max(columns)]), float(ys[max(rows)])],
            "mean_confidence": float(np.mean(
                [confidence[row, column] for row, column in cells])),
            "mean_local_score": float(np.mean(
                [quality[row, column] for row, column in cells])),
            "manual_protected": bool(component.get("manual_protected")),
            "manual_region": int(component.get("manual_region", 0)),
        })
        manual_region = int(component.get("manual_region", 0))
        if manual_region in manual_region_polygons:
            domains[-1]["manual_exact_polygon"] = manual_region_polygons[
                manual_region]
            domains[-1]["manual_boundary_exact"] = True

    # Each neighboring pair of retained domain components gets its own set of
    # boundary segments.  Multiple domain pairs therefore remain explicit.
    boundary_groups = {}
    sample_boundary_segments = []
    dx = float(np.median(np.diff(xs))) if len(xs) > 1 else float(stride)
    dy = float(np.median(np.diff(ys))) if len(ys) > 1 else float(stride)
    for row in range(labels.shape[0]):
        for column in range(labels.shape[1]-1):
            first, second = (int(domain_map[row, column]),
                             int(domain_map[row, column+1]))
            if first <= 0 or second <= 0 or first == second:
                continue
            if float(np.mean([quality[row, column],
                              quality[row, column+1]])) < .08:
                continue
            key = tuple(sorted((first, second)))
            x = float(xs[column]+xs[column+1])/2.0
            boundary_groups.setdefault(key, []).append(
                [x, float(ys[row])-dy/2.0, x, float(ys[row])+dy/2.0])
            
    # Sample/background interfaces are a different physical object from
    # domain/domain interfaces.  Preserve them separately so the GUI can draw
    # and label a dedicated ``sample boundary`` without implying a lattice
    # transition.  The grid edge is refined against the raw TEM below.
    if accept_background:
        for row in range(sample_mask.shape[0]):
            for column in range(sample_mask.shape[1]-1):
                if bool(sample_mask[row, column]) == bool(
                        sample_mask[row, column+1]):
                    continue
                x = float(xs[column]+xs[column+1])/2.0
                sample_boundary_segments.append(
                    [x, float(ys[row])-dy/2.0,
                     x, float(ys[row])+dy/2.0])
    for row in range(labels.shape[0]-1):
        for column in range(labels.shape[1]):
            first, second = (int(domain_map[row, column]),
                             int(domain_map[row+1, column]))
            if first <= 0 or second <= 0 or first == second:
                continue
            if float(np.mean([quality[row, column],
                              quality[row+1, column]])) < .08:
                continue
            key = tuple(sorted((first, second)))
            y = float(ys[row]+ys[row+1])/2.0
            boundary_groups.setdefault(key, []).append(
                [float(xs[column])-dx/2.0, y,
                 float(xs[column])+dx/2.0, y])
    if accept_background:
        for row in range(sample_mask.shape[0]-1):
            for column in range(sample_mask.shape[1]):
                if bool(sample_mask[row, column]) == bool(
                        sample_mask[row+1, column]):
                    continue
                y = float(ys[row]+ys[row+1])/2.0
                sample_boundary_segments.append(
                    [float(xs[column])-dx/2.0, y,
                     float(xs[column])+dx/2.0, y])
    boundary_evidence, boundary_evidence_scale = (
        _tem_texture_boundary_evidence(gray) if boundary_groups else
        (None, 1.0))
    dark_seam_evidence, dark_seam_evidence_scale = (
        _tem_dark_seam_boundary_evidence(gray) if boundary_groups else
        (None, 1.0))
    sample_evidence, sample_evidence_scale = (
        _tem_precise_sample_boundary_evidence(gray)
        if sample_boundary_segments else
        (None, 1.0))
    if sample_boundary_segments:
        sample_boundary_segments = _refine_boundary_segments_from_tem(
            sample_boundary_segments, sample_evidence,
            sample_evidence_scale, search_px=max(dx, dy)*1.35)
    sample_label_segment = (sample_boundary_segments[
        len(sample_boundary_segments)//2]
                            if sample_boundary_segments else None)
    sample_boundary = {
        "valid": bool(sample_boundary_segments),
        "segments": sample_boundary_segments,
        "label_x": ((float(sample_label_segment[0])+
                     float(sample_label_segment[2]))/2.0
                    if sample_label_segment else None),
        "label_y": ((float(sample_label_segment[1])+
                     float(sample_label_segment[3]))/2.0
                    if sample_label_segment else None),
        "boundary_method": "tem_texture_sample_background_boundary",
    }
    radial_sample_boundary = _radial_sample_boundary(gray, scale_bar)
    manual_background_edit = manual_background_requested
    if (radial_sample_boundary.get("valid") and
            not manual_background_edit):
        sample_boundary = radial_sample_boundary
    if manual_background_polygons:
        exact_segments = []
        for polygon in manual_background_polygons:
            points = polygon + [polygon[0]]
            exact_segments.extend([
                [float(first[0]), float(first[1]),
                 float(second[0]), float(second[1])]
                for first, second in zip(points[:-1], points[1:])])
        sample_boundary.setdefault("segments", []).extend(exact_segments)
        sample_boundary["valid"] = True
        sample_boundary["manual_exact_polygons"] = manual_background_polygons
        sample_boundary["boundary_method"] = (
            "automatic_plus_frozen_manual_sample_background_boundaries")
    provisional_boundaries = []
    evidence_p65 = (float(np.percentile(boundary_evidence, 65.0))
                    if boundary_evidence is not None else 0.0)
    evidence_p75 = (float(np.percentile(boundary_evidence, 75.0))
                    if boundary_evidence is not None else 0.0)
    domain_area = {int(domain["domain_id"]): float(
        domain.get("area_fraction", 0.0)) for domain in domains}
    domain_lookup = {int(domain["domain_id"]): domain for domain in domains}
    for domain_pair, segments in sorted(boundary_groups.items()):
        # Local FFT centres leave an unavoidable half-patch margin.  If an
        # open seam ends within that margin, continue it to the nearest frame
        # edge before the accurate TEM evidence trace below.  This closes the
        # two measured domain areas without altering internal/closed seams.
        first_candidate = domain_lookup[int(domain_pair[0])]
        second_candidate = domain_lookup[int(domain_pair[1])]
        manual_polygon = (first_candidate.get("manual_exact_polygon") or
                          second_candidate.get("manual_exact_polygon") or [])
        manual_exact_boundary = len(manual_polygon) >= 3
        if manual_exact_boundary:
            points = manual_polygon + [manual_polygon[0]]
            segments = [[float(first[0]), float(first[1]),
                         float(second[0]), float(second[1])]
                        for first, second in zip(points[:-1], points[1:])]
        else:
            segments = _extend_open_boundary_to_frame(
                segments, width, height,
                maximum_gap=max(dx, dy)*2.40)
        first_candidate_a = float(
            first_candidate.get("lattice_constant_px") or 0.0)
        second_candidate_a = float(
            second_candidate.get("lattice_constant_px") or 0.0)
        candidate_a_difference = (
            abs(first_candidate_a-second_candidate_a) /
            max((first_candidate_a+second_candidate_a)/2.0, 1e-9))
        candidate_same_lattice_twin = bool(
            str(first_candidate.get("symmetry") or "Unknown") ==
            str(second_candidate.get("symmetry") or "Unknown") and
            candidate_a_difference < .08)
        pair_evidence = (dark_seam_evidence
                         if candidate_same_lattice_twin
                         else boundary_evidence)
        pair_evidence_scale = (dark_seam_evidence_scale
                               if candidate_same_lattice_twin
                               else boundary_evidence_scale)
        # Local FFT supplies domain identity; the raw TEM supplies the visible
        # seam position.  No spline/corner smoothing is applied after this
        # edge-aware refinement.
        monotone_seam = (_trace_monotone_tem_seam(
            segments, pair_evidence, pair_evidence_scale,
            search_pixels=max(dx, dy)*1.05)
                         if (candidate_same_lattice_twin and
                             not manual_exact_boundary) else None)
        if monotone_seam:
            segments = monotone_seam
        elif not manual_exact_boundary:
            segments = _refine_boundary_segments_from_tem(
                segments, pair_evidence, pair_evidence_scale,
                search_px=max(dx, dy)*1.05,
                # Domain seams must remain attached to the measured
                # dark/texture discontinuity.  The former five corner-cut
                # passes plus a broad RDP tolerance visibly displaced the
                # boundary toward the local-FFT cell centre.
                maximum_search=64,
                offset_penalty=(.0015 if candidate_same_lattice_twin
                                else .004),
                transition_penalty=(.065 if candidate_same_lattice_twin
                                    else .045),
                smoothing_iterations=(2 if candidate_same_lattice_twin
                                      else 1),
                reduction_tolerance=(.045 if candidate_same_lattice_twin
                                     else .065))
        visible_values = []
        if boundary_evidence is not None:
            evidence_height, evidence_width = boundary_evidence.shape
            for x1, y1, x2, y2 in segments:
                count = max(2, int(round(max(abs(x2-x1), abs(y2-y1)) /
                                         boundary_evidence_scale))+1)
                for x, y in zip(np.linspace(x1, x2, count),
                                np.linspace(y1, y2, count)):
                    column = int(round(x/boundary_evidence_scale))
                    row = int(round(y/boundary_evidence_scale))
                    if (0 <= row < evidence_height and
                            0 <= column < evidence_width):
                        visible_values.append(float(
                            boundary_evidence[row, column]))
        mean_visible = float(np.mean(visible_values)) if visible_values else 0.0
        visible_fraction = (float(np.mean(
            np.asarray(visible_values) >= evidence_p75))
                            if visible_values else 0.0)
        minimum_area = min(domain_area.get(int(domain_pair[0]), 0.0),
                           domain_area.get(int(domain_pair[1]), 0.0))
        # A distinct domain requires a continuous visible TEM interface, not
        # merely another strong/higher-order FFT hypothesis.  Small domains
        # need stronger image evidence because isolated local-FFT flips are
        # common in noisy or slowly varying single-domain micrographs.
        visible_boundary = bool(
            (minimum_area >= .08 and visible_fraction >= .30 and
             mean_visible >= evidence_p65) or
            (minimum_area < .08 and visible_fraction >= .50 and
             mean_visible >= evidence_p75))
        first_domain = domain_lookup[int(domain_pair[0])]
        second_domain = domain_lookup[int(domain_pair[1])]
        first_symmetry = str(first_domain.get("symmetry") or "Unknown")
        second_symmetry = str(second_domain.get("symmetry") or "Unknown")
        symmetry_difference = first_symmetry != second_symmetry
        symmetry_period = (90.0 if first_symmetry == second_symmetry ==
                           "Square" else 60.0)
        orientation_difference = _circular_distance(
            float(first_domain.get("orientation_deg", 0.0)),
            float(second_domain.get("orientation_deg", 0.0)),
            symmetry_period)
        first_a = float(first_domain.get("lattice_constant_px") or 0.0)
        second_a = float(second_domain.get("lattice_constant_px") or 0.0)
        relative_a_difference = (abs(first_a-second_a) /
                                 max((first_a+second_a)/2.0, 1e-9))
        minimum_confidence = min(
            float(first_domain.get("mean_confidence", 0.0)),
            float(second_domain.get("mean_confidence", 0.0)))
        same_lattice_twin = bool(
            not symmetry_difference and relative_a_difference < .08)
        # A same-a/same-symmetry split is the most common false domain in a
        # single periodic field: a second reciprocal hypothesis can win on
        # alternating local windows even though no physical seam exists.
        # Require a continuously visible real-space interface for such twins.
        # The regression twin image has a dark, spatially coherent seam and
        # comfortably passes these gates; the Moire field in ``2.tif`` does
        # not and is therefore kept as one specimen/domain.
        if same_lattice_twin:
            visible_boundary = bool(
                visible_boundary and minimum_confidence >= .62 and
                visible_fraction >= .46)
        elif minimum_area < .08 and minimum_confidence < .40:
            # Tiny low-confidence islands in a high-resolution lattice are
            # normally a local harmonic/contrast fluctuation.  Even a strong
            # edge elsewhere in their coarse cell must not promote them to a
            # separately labelled crystallographic domain.
            visible_boundary = False
        # A weak or invisible real-space seam must not erase a genuine domain.
        # Strong, spatially coherent local-FFT evidence can independently
        # support a domain when symmetry, lattice constant or orientation is
        # clearly different.  Confidence/area gates keep legal higher-order
        # reflections in a uniform Honeycomb/Kagome image from being promoted
        # to extra domains.
        fft_domain_evidence = bool(
            minimum_area >= .04 and minimum_confidence >= .48 and (
                symmetry_difference or relative_a_difference >= .14 or
                (orientation_difference >= 7.5 and
                 minimum_confidence >= (.62 if same_lattice_twin else .48))))
        evidence_basis = []
        if visible_boundary:
            evidence_basis.append("tem_boundary")
        if fft_domain_evidence and symmetry_difference:
            evidence_basis.append("lattice_symmetry")
        if fft_domain_evidence and relative_a_difference >= .14:
            evidence_basis.append("lattice_constant")
        if fft_domain_evidence and orientation_difference >= 7.5:
            evidence_basis.append("orientation")
        supported_boundary = bool(visible_boundary or fft_domain_evidence)
        manual_boundary = bool(first_domain.get("manual_protected") or
                               second_domain.get("manual_protected"))
        if manual_boundary:
            # The user supplied a spatial prior; crystallographic properties
            # still come from the local FFT and the seam is still snapped to
            # raw TEM evidence, but an automatic confidence gate must not
            # silently merge the explicitly requested domain again.
            supported_boundary = True
            evidence_basis.append("manual_exact_boundary")
        label_segment = segments[len(segments)//2] if segments else None
        provisional_boundaries.append({
            "between_domains": list(domain_pair),
            "segments": segments,
            "visible_boundary": visible_boundary,
            "supported_boundary": supported_boundary,
            "evidence_basis": evidence_basis,
            "orientation_difference_deg": float(orientation_difference),
            "symmetry_difference": bool(symmetry_difference),
            "lattice_constant_relative_difference": float(
                relative_a_difference),
            "minimum_local_fft_confidence": float(minimum_confidence),
            "tem_boundary_mean_score": mean_visible,
            "tem_boundary_high_evidence_fraction": visible_fraction,
            "manual_boundary_exact": bool(manual_exact_boundary),
            "label_x": ((float(label_segment[0])+float(label_segment[2]))/2.0
                        if label_segment else None),
            "label_y": ((float(label_segment[1])+float(label_segment[3]))/2.0
                        if label_segment else None),
        })

    # Merge unsupported local-FFT components.  A domain may be supported by a
    # TEM-visible boundary or by strong coherent FFT evidence (orientation,
    # symmetry or lattice-constant difference).  This preserves weak-boundary
    # twins while preventing uniform Honeycomb/Kagome fields from being split
    # by legal higher-order reflections or contrast drift.
    original_domain_area = dict(domain_area)
    parent = {int(domain["domain_id"]): int(domain["domain_id"])
              for domain in domains}

    def find(domain_id):
        while parent[domain_id] != domain_id:
            parent[domain_id] = parent[parent[domain_id]]
            domain_id = parent[domain_id]
        return domain_id

    def union(first, second):
        first_root, second_root = find(first), find(second)
        if first_root == second_root:
            return
        first_size = domain_area.get(first_root, 0.0)
        second_size = domain_area.get(second_root, 0.0)
        if first_size < second_size:
            first_root, second_root = second_root, first_root
        parent[second_root] = first_root
        domain_area[first_root] = first_size+second_size

    # Connected-component labels describe spatial continuity, not a new
    # crystallographic identity.  When every retained component selects the
    # exact same globally validated FFT family, gaps caused by background,
    # low-confidence tiles or specimen edges must not create D1/D2/... with
    # identical symmetry, a, orientation and reciprocal axes.  Merge those
    # fragments before boundary support is evaluated.  True twins and mixed
    # domains have at least two orientation_index values and are untouched.
    retained_orientation_indices = set(
        int(domain.get("orientation_index", -1)) for domain in domains)
    if len(retained_orientation_indices) == 1 and len(domains) > 1:
        automatic_domains = [domain for domain in domains
                             if not domain.get("manual_protected")]
        if automatic_domains:
            anchor_domain_id = int(automatic_domains[0]["domain_id"])
            for domain in automatic_domains[1:]:
                union(anchor_domain_id, int(domain["domain_id"]))

    visibly_bounded = set()
    for boundary in provisional_boundaries:
        first, second = [int(value) for value in
                         boundary["between_domains"]]
        if (domain_lookup.get(first, {}).get("manual_protected") or
                domain_lookup.get(second, {}).get("manual_protected")):
            visibly_bounded.update((first, second))
            continue
        if boundary["supported_boundary"]:
            visibly_bounded.update((first, second))
        elif (original_domain_area.get(first, 0.0) >= .08 and
              original_domain_area.get(second, 0.0) >= .08):
            # Two substantial FFT regions without a TEM-visible interface are
            # one physical domain (the common single-domain false split).
            union(first, second)

    # A small unsupported component may be absorbed by one largest adjacent
    # domain, but never used as a bridge that accidentally joins two large
    # domains separated by a real boundary.
    unsupported_neighbors = {}
    for boundary in provisional_boundaries:
        if boundary["supported_boundary"]:
            continue
        first, second = [int(value) for value in
                         boundary["between_domains"]]
        unsupported_neighbors.setdefault(first, set()).add(second)
        unsupported_neighbors.setdefault(second, set()).add(first)
    for domain_id, area in sorted(original_domain_area.items(),
                                  key=lambda item: item[1]):
        if area >= .08 or domain_id in visibly_bounded:
            continue
        options = unsupported_neighbors.get(domain_id) or set()
        if not options:
            continue
        target = max(options, key=lambda value: original_domain_area.get(
            value, 0.0))
        union(domain_id, target)

    merged_groups = {}
    for domain in domains:
        root = find(int(domain["domain_id"]))
        merged_groups.setdefault(root, []).append(domain)
    merged_domains, root_to_new = [], {}
    ordered_groups = sorted(merged_groups.items(), key=lambda item: sum(
        int(domain.get("cell_count", 0)) for domain in item[1]), reverse=True)
    for new_id, (root, members) in enumerate(ordered_groups, 1):
        dominant = max(members, key=lambda domain: int(
            domain.get("cell_count", 0)))
        merged = dict(dominant)
        merged["domain_id"] = new_id
        merged["cell_count"] = sum(int(domain.get("cell_count", 0))
                                   for domain in members)
        merged["area_fraction"] = sum(float(domain.get(
            "area_fraction", 0.0)) for domain in members)
        merged["bbox"] = [
            min(float(domain["bbox"][0]) for domain in members),
            min(float(domain["bbox"][1]) for domain in members),
            max(float(domain["bbox"][2]) for domain in members),
            max(float(domain["bbox"][3]) for domain in members),
        ]
        merged["merged_domain_ids"] = [int(domain["domain_id"])
                                       for domain in members]
        merged_domains.append(merged)
        root_to_new[root] = new_id
    domains = merged_domains
    retained_area = sum(float(domain.get("area_fraction", 0.0))
                        for domain in domains) or 1.0
    for domain in domains:
        domain["area_fraction"] = float(
            domain.get("area_fraction", 0.0))/retained_area

    # Preserve the exact spatial cells used for the reported area fractions.
    # The UI uses these run-length rectangles for a light transparent overlay,
    # so the coloured region and the numerical percentage always describe the
    # same measurement rather than two independently inferred masks.
    final_domain_map = np.full(domain_map.shape, -1, dtype=int)
    for row in range(domain_map.shape[0]):
        for column in range(domain_map.shape[1]):
            old_domain = int(domain_map[row, column])
            if old_domain <= 0 or old_domain not in parent:
                continue
            final_domain_map[row, column] = int(
                root_to_new[find(old_domain)])
    x_bounds = np.concatenate((
        [0.0], (xs[:-1].astype(float)+xs[1:].astype(float))/2.0,
        [float(width)]))
    y_bounds = np.concatenate((
        [0.0], (ys[:-1].astype(float)+ys[1:].astype(float))/2.0,
        [float(height)]))
    area_runs = {int(domain["domain_id"]): [] for domain in domains}
    for row in range(final_domain_map.shape[0]):
        column = 0
        while column < final_domain_map.shape[1]:
            domain_id = int(final_domain_map[row, column])
            if domain_id <= 0:
                column += 1
                continue
            end = column+1
            while (end < final_domain_map.shape[1] and
                   int(final_domain_map[row, end]) == domain_id):
                end += 1
            area_runs.setdefault(domain_id, []).append([
                float(x_bounds[column]), float(y_bounds[row]),
                float(x_bounds[end]), float(y_bounds[row+1])])
            column = end
    for domain in domains:
        domain["area_runs"] = area_runs.get(
            int(domain["domain_id"]), [])
        domain["area_overlay_basis"] = (
            "local_fft_cells_used_for_area_fraction")

    combined_boundaries = {}
    for boundary in provisional_boundaries:
        if not boundary["supported_boundary"]:
            continue
        old_first, old_second = [int(value) for value in
                                 boundary["between_domains"]]
        first, second = root_to_new[find(old_first)], root_to_new[find(old_second)]
        if first == second:
            continue
        key = tuple(sorted((first, second)))
        item = combined_boundaries.setdefault(key, {
            "segments": [], "scores": [], "fractions": [],
            "visible": [], "evidence": [], "orientation": [],
            "symmetry": [], "a_difference": [], "confidence": [],
            "manual_exact": []})
        item["segments"].extend(boundary["segments"])
        item["scores"].append(boundary["tem_boundary_mean_score"])
        item["fractions"].append(
            boundary["tem_boundary_high_evidence_fraction"])
        item["visible"].append(bool(boundary["visible_boundary"]))
        item["evidence"].extend(boundary.get("evidence_basis") or [])
        item["orientation"].append(boundary[
            "orientation_difference_deg"])
        item["symmetry"].append(boundary["symmetry_difference"])
        item["a_difference"].append(boundary[
            "lattice_constant_relative_difference"])
        item["confidence"].append(boundary[
            "minimum_local_fft_confidence"])
        item["manual_exact"].append(bool(
            boundary.get("manual_boundary_exact")))
    boundaries = []
    for boundary_index, (domain_pair, item) in enumerate(
            sorted(combined_boundaries.items()), 1):
        segments = item["segments"]
        label_segment = segments[len(segments)//2] if segments else None
        boundaries.append({
            "boundary_id": boundary_index,
            "between_domains": list(domain_pair),
            "segments": segments,
            "visible_boundary": bool(any(item["visible"])),
            "supported_boundary": True,
            "evidence_basis": sorted(set(item["evidence"])),
            "orientation_difference_deg": float(max(item["orientation"])),
            "symmetry_difference": bool(any(item["symmetry"])),
            "lattice_constant_relative_difference": float(max(
                item["a_difference"])),
            "minimum_local_fft_confidence": float(max(item["confidence"])),
            "manual_boundary_exact": bool(any(item["manual_exact"])),
            "tem_boundary_mean_score": float(np.mean(item["scores"])),
            "tem_boundary_high_evidence_fraction": float(np.mean(
                item["fractions"])),
            "label_x": ((float(label_segment[0])+float(label_segment[2]))/2.0
                        if label_segment else None),
            "label_y": ((float(label_segment[1])+float(label_segment[3]))/2.0
                        if label_segment else None),
        })
    has_manual_exact_domain = any(
        domain.get("manual_boundary_exact") for domain in domains)
    if not has_manual_exact_domain:
        _precise_two_domain_area_polygons(
            domains, boundaries, width, height)
    for domain in domains:
        polygon = domain.get("manual_exact_polygon") or []
        if len(polygon) >= 3:
            domain["area_polygons"] = [polygon]
            domain["area_overlay_basis"] = "frozen_manual_domain_boundary"
    if has_manual_exact_domain:
        # Recalculate percentages on a finer integration grid.  Automatic
        # domains retain their prior classification, but every pixel inside a
        # frozen manual polygon is reassigned to that exact polygon before
        # counting.  Thus analysis and the visible overlay use the identical
        # user-drawn boundary, not the coarse local-FFT cell approximation.
        integration_step = max(1.0, max(float(width), float(height))/760.0)
        integration_x = np.arange(
            integration_step*.5, float(width), integration_step)
        integration_y = np.arange(
            integration_step*.5, float(height), integration_step)
        integration_map = np.full(
            (len(integration_y), len(integration_x)), -1, dtype=int)
        for domain in domains:
            domain_id = int(domain["domain_id"])
            for run in domain.get("area_runs") or []:
                if len(run) < 4:
                    continue
                rows = ((integration_y >= float(run[1])) &
                        (integration_y < float(run[3])))
                columns = ((integration_x >= float(run[0])) &
                           (integration_x < float(run[2])))
                integration_map[np.ix_(rows, columns)] = domain_id
        for domain in domains:
            polygon = domain.get("manual_exact_polygon") or []
            if len(polygon) < 3:
                continue
            exact_mask = _polygon_mask_at_analysis_resolution(
                polygon, width, height, maximum_side=760)
            if exact_mask is not None and exact_mask.shape == integration_map.shape:
                integration_map[exact_mask] = int(domain["domain_id"])
        for polygon in manual_background_polygons:
            exact_mask = _polygon_mask_at_analysis_resolution(
                polygon, width, height, maximum_side=760)
            if exact_mask is not None and exact_mask.shape == integration_map.shape:
                integration_map[exact_mask] = -1
        counts = {int(domain["domain_id"]): int(np.count_nonzero(
            integration_map == int(domain["domain_id"])))
                  for domain in domains}
        count_total = sum(counts.values())
        if count_total > 0:
            for domain in domains:
                domain["area_fraction"] = float(
                    counts[int(domain["domain_id"])])/float(count_total)
                domain["area_normalization_basis"] = (
                    "frozen_manual_boundary_fine_grid")
    else:
        _normalize_domain_areas_to_sample_boundary(
            domains, sample_boundary, width, height)
    return {
        "valid": True,
        "method": "overlapping_local_fft_symmetry_a_orientation_domains",
        "boundary_method": (
            "tem_boundary_plus_fft_symmetry_a_orientation_evidence"),
        "patch_px": patch, "stride_px": stride,
        "grid_x": xs.astype(int).tolist(),
        "grid_y": ys.astype(int).tolist(),
        "orientation_count": len(set(int(domain.get(
            "orientation_index", -1)) for domain in domains)),
        "family_count": len(set(int(domain.get(
            "orientation_index", -1)) for domain in domains)),
        "domain_count": len(domains),
        "boundary_count": len(boundaries),
        "sample_cell_fraction": float(np.mean(sample_mask)),
        "background_fraction": float(1.0-np.mean(sample_mask)),
        "background_excluded": bool(accept_background),
        "automatic_background_rejected": bool(
            automatic_background_rejected),
        "sample_periodicity_validation": sample_periodicity_validation,
        "sample_boundary": sample_boundary,
        "domains": domains,
        "boundaries": boundaries,
    }


def _precise_fft_assets(image, output_dir, analysis_kind="bilayer"):
    """Detect visible reciprocal peaks and preserve their complex phase.

    The detector uses the real FFT of the uploaded TEM, sub-pixel peak
    centroids, compact locally fitted ellipses, and exact Friedel partners.
    Raw arrays keep the source width:height ratio.  The GUI displays them on
    isotropic reciprocal-space axes, so FFT geometry is not stretched.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gray = np.asarray(image, dtype=np.float64) / 255.0
    height, width = gray.shape
    center_y, center_x = height // 2, width // 2
    window = np.outer(np.hanning(height), np.hanning(width))
    # Remove slow illumination gradients and apply only sub-pixel Gaussian
    # denoising before peak detection.  Peak coordinates are then refined on
    # the unfiltered magnitude below, so contrast enhancement cannot move a
    # reciprocal-lattice point.
    illumination = _gaussian_filter_fft(gray, 28.0)
    detection_gray = _gaussian_filter_fft(gray - illumination, 0.65)
    detection_spectrum = np.fft.fftshift(np.fft.fft2(
        detection_gray * window))
    raw_windowed_spectrum = np.fft.fftshift(np.fft.fft2(
        (gray - float(np.mean(gray))) * window))
    magnitude = np.abs(raw_windowed_spectrum)
    log_magnitude = np.log1p(np.abs(detection_spectrum))
    response = (_gaussian_filter_fft(log_magnitude, 1.25) -
                _gaussian_filter_fft(log_magnitude, 12.0))
    yy, xx = np.mgrid[:height, :width]
    radius = np.sqrt((xx - center_x) ** 2 + (yy - center_y) ** 2)
    # Do not discard the crystallographic first-order shell in high-resolution
    # cryo-EM images.  In those images the helix/sub-lattice reflections can
    # be much farther from the origin than the weaker pore/unit-cell shell;
    # the old 8% central exclusion removed the latter completely and forced
    # Honeycomb/Kagome a and IFFT to use a higher-order shell.  Slow background
    # is already suppressed above, so a conservative 2.5% exclusion retains
    # the unit-cell peaks without admitting the DC/illumination halo.
    minimum_radius = max(18.0, min(width, height) * 0.025)
    maximum_radius = min(width, height) * 0.46
    valid = (radius >= minimum_radius) & (radius <= maximum_radius)
    threshold = float(np.percentile(response[valid], 99.84))
    coordinates = np.argwhere(valid & (response >= threshold))
    raw = sorted(
        ({"x": int(x), "y": int(y),
          "score": float(response[y, x])} for y, x in coordinates),
        key=lambda point: point["score"], reverse=True)
    minimum_separation = max(6.0, min(width, height) * 0.0075)
    maxima = []
    for point in raw:
        if any((point["x"] - old["x"]) ** 2 +
               (point["y"] - old["y"]) ** 2 < minimum_separation ** 2
               for old in maxima):
            continue
        maxima.append(point)
        # Hundreds of paired candidates are already far more than the clear
        # crystallographic reflections used by the downstream indexer.  A
        # 4k cryo/TEM frame can otherwise contribute thousands of noise/local
        # maxima, making Friedel pairing and ellipse fitting appear hung while
        # adding no valid lattice information.  Scores are descending and the
        # FFT is centrosymmetric, so 480 maxima retain both members of every
        # clearly visible pair with a generous margin.
        if len(maxima) >= 480:
            break

    paired = []
    used = set()
    mirror_tolerance = max(5.0, minimum_separation * 0.8)
    for index, point in enumerate(maxima):
        if index in used:
            continue
        mirror_x = 2 * center_x - point["x"]
        mirror_y = 2 * center_y - point["y"]
        options = [
            (other_index,
             (other["x"] - mirror_x) ** 2 +
             (other["y"] - mirror_y) ** 2)
            for other_index, other in enumerate(maxima)
            if other_index not in used and other_index != index]
        if not options:
            continue
        other_index, distance_squared = min(options, key=lambda item: item[1])
        if distance_squared > mirror_tolerance ** 2:
            continue
        other = maxima[other_index]
        x = int(round((point["x"] + (2 * center_x - other["x"])) / 2.0))
        y = int(round((point["y"] + (2 * center_y - other["y"])) / 2.0))
        score = float((point["score"] + other["score"]) / 2.0)
        paired.extend((
            {"x": x, "y": y, "score": score},
            {"x": 2 * center_x - x, "y": 2 * center_y - y,
             "score": score}))
        used.update((index, other_index))

    # Fit a compact ellipse to the connected half-height component around
    # each actual maximum, then refine the peak to a sub-pixel centroid.
    fit_radius = max(10, int(round(min(width, height) * 0.0135)))
    for peak in paired:
        x, y = int(peak["x"]), int(peak["y"])
        if (x < fit_radius or y < fit_radius or
                x >= width - fit_radius or y >= height - fit_radius):
            peak.update({"rx": 5.0, "ry": 3.5, "angle": 0.0})
            continue
        patch = magnitude[y-fit_radius:y+fit_radius+1,
                          x-fit_radius:x+fit_radius+1]
        grid_y, grid_x = np.mgrid[-fit_radius:fit_radius+1,
                                  -fit_radius:fit_radius+1]
        annulus = grid_x * grid_x + grid_y * grid_y >= (fit_radius * .72) ** 2
        background = float(np.median(patch[annulus]))
        center_value = float(patch[fit_radius, fit_radius])
        local_threshold = background + .34 * max(
            center_value - background, 1e-12)
        binary = patch >= local_threshold
        component = np.zeros_like(binary, dtype=bool)
        stack = [(fit_radius, fit_radius)]
        while stack:
            local_y, local_x = stack.pop()
            if (local_y < 0 or local_x < 0 or
                    local_y >= binary.shape[0] or
                    local_x >= binary.shape[1] or
                    component[local_y, local_x] or
                    not binary[local_y, local_x]):
                continue
            component[local_y, local_x] = True
            stack.extend(((local_y-1, local_x), (local_y+1, local_x),
                          (local_y, local_x-1), (local_y, local_x+1)))
        if not np.any(component):
            component = grid_x * grid_x + grid_y * grid_y <= 16
        weights = np.maximum(patch - background, 0.0) * component
        total = float(np.sum(weights))
        if total <= 0:
            peak.update({"rx": 5.0, "ry": 3.5, "angle": 0.0})
            continue
        mean_x = float(np.sum(weights * grid_x) / total)
        mean_y = float(np.sum(weights * grid_y) / total)
        dx, dy = grid_x - mean_x, grid_y - mean_y
        covariance = np.asarray([
            [np.sum(weights * dx * dx) / total,
             np.sum(weights * dx * dy) / total],
            [np.sum(weights * dx * dy) / total,
             np.sum(weights * dy * dy) / total]])
        values, vectors = np.linalg.eigh(covariance)
        order = np.argsort(values)[::-1]
        values, vectors = values[order], vectors[:, order]
        maximum_rx = max(8.0, min(width, height) * .0125)
        maximum_ry = max(6.0, min(width, height) * .0087)
        rx = float(np.clip(2.25 * math.sqrt(max(values[0], 1e-9)) + 1.5,
                           3.5, maximum_rx))
        ry = float(np.clip(2.25 * math.sqrt(max(values[1], 1e-9)) + 1.5,
                           2.5, maximum_ry))
        angle = float(math.degrees(math.atan2(vectors[1, 0], vectors[0, 0])))
        peak.update({"x": round(x + mean_x, 3),
                     "y": round(y + mean_y, 3),
                     "rx": round(rx, 3), "ry": round(ry, 3),
                     "angle": round(angle, 3)})
    for index in range(0, len(paired), 2):
        first, second = paired[index], paired[index + 1]
        x = (float(first["x"]) + (2 * center_x - float(second["x"]))) / 2.0
        y = (float(first["y"]) + (2 * center_y - float(second["y"]))) / 2.0
        first["x"], first["y"] = round(x, 3), round(y, 3)
        second["x"], second["y"] = round(2 * center_x - x, 3), round(
            2 * center_y - y, 3)

    # Retain every clearly visible, locally fitted Friedel peak.  The primary
    # lattice classifier below deliberately reduces ``paired`` to one quartet,
    # but a validated second a (for example a pore/supercell repeat) needs its
    # own precise spots for the combined phase-preserving IFFT.
    all_visible_spots = [dict(point) for point in paired]

    dominant_single = None
    if analysis_kind == "single" and len(paired) >= 4:
        # Single-layer analysis reports the helix-scale lattice a.  A single
        # Square orientation contributes two orthogonal reciprocal axes, each
        # represented by an exact Friedel pair.  Select one strongest,
        # period-matched orthogonal quartet; do not reinterpret the two axes as
        # two lattice orientations and do not replace it with a weaker inner
        # supercell/repeat shell nearer the FFT origin.
        axes = []
        for pair_index in range(0, len(paired), 2):
            point = paired[pair_index]
            fx = (float(point["x"]) - center_x) / max(float(width), 1.0)
            fy = (float(point["y"]) - center_y) / max(float(height), 1.0)
            frequency = math.hypot(fx, fy)
            if frequency <= 0:
                continue
            axes.append({
                "pair_index": pair_index,
                "period": 1.0 / frequency,
                "angle": math.degrees(math.atan2(-fy, fx)) % 180.0,
                "score": float(point["score"]),
            })
        square_options = []
        for first_index, first in enumerate(axes):
            for second in axes[first_index + 1:]:
                orthogonal_error = abs(
                    _circular_distance(first["angle"], second["angle"], 180.0) -
                    90.0)
                period_error = abs(first["period"] - second["period"]) / max(
                    (first["period"] + second["period"]) / 2.0, 1e-9)
                if orthogonal_error > 12.0 or period_error > 0.20:
                    continue
                mean_period = (first["period"] + second["period"]) / 2.0
                merit = (math.sqrt(first["score"] * second["score"]) /
                         (1.0 + orthogonal_error / 5.0 + period_error * 4.0))
                square_options.append((merit, first, second))
        best = None
        if square_options:
            # Square images can contain strong lower-frequency pore/supercell
            # shells.  Keep detecting those shells, but define the primary
            # helix lattice from the highest-frequency orthogonal family that
            # remains almost as strong as the best family.  This restores the
            # established helix-a + separately indexed supercell behaviour
            # while the enlarged central search remains available to H/K.
            maximum_merit = max(option[0] for option in square_options)
            credible = [option for option in square_options
                        if option[0] >= maximum_merit * .75]
            best = min(credible, key=lambda option: (
                (option[1]["period"]+option[2]["period"])/2.0,
                -option[0]))
        if best is not None:
            unused_merit, first, second = best
            selected = []
            for axis in (first, second):
                start = int(axis["pair_index"])
                selected.extend((paired[start], paired[start + 1]))
            paired = selected
            phases = [first["angle"] % 90.0, second["angle"] % 90.0]
            dominant_single = {
                "valid": True,
                "orientation_deg": float(_circular_mean_90(
                    phases, [first["score"], second["score"]])),
                "lattice_constant_px": float(np.median(
                    [first["period"], second["period"]])),
                "axis_angles_deg": [float(first["angle"]),
                                    float(second["angle"])],
                "axis_periods_px": [float(first["period"]),
                                    float(second["period"])],
                "axis_scores": [float(first["score"]),
                                float(second["score"])],
                "peaks": [dict(point) for point in paired],
                "orthogonal_error_deg": float(abs(
                    _circular_distance(first["angle"], second["angle"], 180.0) -
                    90.0)),
                "period_mismatch_fraction": float(abs(
                    first["period"] - second["period"]) / max(
                        (first["period"] + second["period"]) / 2.0, 1e-9)),
            }
        else:
            # Never pass hundreds of unclassified maxima into IFFT.  A small
            # strongest-only diagnostic mask is safer and is flagged invalid.
            paired = paired[:min(8, len(paired))]
            dominant_single = {
                "valid": False,
                "error": "No period-consistent orthogonal first-order Square peaks were found.",
            }

    # Select every clearly resolved FFT maximum for the phase-preserving IFFT,
    # while retaining Friedel pairs as indivisible units.  This is purposely
    # broader than the one/two first-order quartets used to measure a.
    clear_visible_spots = []
    pair_scores = [float(all_visible_spots[index].get("score", 0.0))
                   for index in range(0, len(all_visible_spots), 2)]
    clear_threshold = (max(pair_scores) * .32 if pair_scores else 0.0)
    for index in range(0, len(all_visible_spots), 2):
        pair = all_visible_spots[index:index+2]
        if len(pair) < 2 or float(pair[0].get("score", 0.0)) < clear_threshold:
            continue
        for point in pair:
            marked = dict(point)
            marked["lattice_role"] = "visible_fft"
            clear_visible_spots.append(marked)
    if not clear_visible_spots:
        clear_visible_spots = [dict(point) for point in paired]

    # Detect a second Square orientation at the same a.  A twin contributes a
    # second period-matched orthogonal quartet, not merely another harmonic.
    single_lattices = ([dict(dominant_single)] if dominant_single and
                       dominant_single.get("valid") else [])
    if single_lattices:
        expected = float(dominant_single["lattice_constant_px"])
        primary_floor = .55 * min(dominant_single.get(
            "axis_scores") or [0.0])
        axes = []
        for pair_index in range(0, len(all_visible_spots), 2):
            point = all_visible_spots[pair_index]
            fx = (float(point["x"])-center_x) / max(float(width), 1.0)
            fy = (float(point["y"])-center_y) / max(float(height), 1.0)
            frequency = math.hypot(fx, fy)
            if frequency <= 0:
                continue
            period = 1.0 / frequency
            score = float(point.get("score", 0.0))
            if (abs(period-expected) / expected > .12 or
                    score < primary_floor):
                continue
            candidate = {
                "pair_index": pair_index, "period": period, "score": score,
                "angle": math.degrees(math.atan2(-fy, fx)) % 180.0}
            if any(_circular_distance(candidate["angle"], old["angle"],
                                      180.0) < 2.5 and
                   abs(candidate["period"]-old["period"]) /
                   max(old["period"], 1e-9) < .08 for old in axes):
                continue
            axes.append(candidate)
        family_candidates = []
        for first_index, first in enumerate(axes):
            for second in axes[first_index+1:]:
                orthogonal_error = abs(_circular_distance(
                    first["angle"], second["angle"], 180.0)-90.0)
                period_error = abs(first["period"]-second["period"]) / max(
                    (first["period"]+second["period"])/2.0, 1e-9)
                if orthogonal_error > 5.0 or period_error > .08:
                    continue
                orientation = float(_circular_mean_90(
                    [first["angle"] % 90.0, second["angle"] % 90.0],
                    [first["score"], second["score"]]))
                merit = math.sqrt(first["score"]*second["score"]) / (
                    1.0+orthogonal_error/3.0+period_error*5.0)
                peaks = []
                for axis in (first, second):
                    start = int(axis["pair_index"])
                    peaks.extend(dict(point) for point in
                                 all_visible_spots[start:start+2])
                family_candidates.append({
                    "valid": True, "orientation_deg": orientation,
                    "lattice_constant_px": float(np.median(
                        [first["period"], second["period"]])),
                    "axis_angles_deg": [float(first["angle"]),
                                        float(second["angle"])],
                    "axis_periods_px": [float(first["period"]),
                                        float(second["period"])],
                    "axis_scores": [float(first["score"]),
                                    float(second["score"])],
                    "orthogonal_error_deg": orthogonal_error,
                    "period_mismatch_fraction": period_error,
                    "peaks": peaks, "merit": merit})
        for family in sorted(family_candidates,
                             key=lambda value: value["merit"], reverse=True):
            if any(_circular_distance(family["orientation_deg"],
                                      old["orientation_deg"], 90.0) < 6.0
                   for old in single_lattices):
                continue
            family["lattice_role"] = "twin_orientation"
            single_lattices.append(family)
            # A polycrystalline/twinned single layer can contain more than two
            # stable Square orientations.  Keep every independently validated
            # orthogonal quartet, with a conservative cap for pathological
            # noisy FFTs.
            if len(single_lattices) >= 6:
                break
        if len(single_lattices) > 1:
            for family in single_lattices[1:]:
                twin_peaks = family.get("peaks") or []
                for point in clear_visible_spots:
                    if any((float(point["x"])-float(twin["x"]))**2 +
                           (float(point["y"])-float(twin["y"]))**2 < 16.0
                           for twin in twin_peaks):
                        point["lattice_role"] = "twin_orientation"

    # Dark, denoised FFT display. The PGM keeps the exact original aspect.
    low, high = np.percentile(log_magnitude, (72.0, 99.97))
    display = np.clip((log_magnitude - low) / max(high - low, 1e-12), 0, 1)
    display = np.power(display, 1.55)
    fft_path = output_dir / "measured_fft_high_contrast.pgm"
    write_pgm(fft_path, display * 255.0)

    # Before a supercell has been validated, a single-layer reconstruction
    # uses only crystallographically classified first-order quartets: four
    # peaks for one Square orientation, or two quartets (eight peaks) for a
    # detected twin.  Unindexed visible maxima must never enter the IFFT.
    if analysis_kind == "single" and single_lattices:
        ifft_spots = []
        for family_index, family in enumerate(single_lattices[:6]):
            for point in family.get("peaks") or []:
                marked = dict(point)
                marked["lattice_role"] = (
                    "twin_orientation" if family_index else "primary_a")
                if any((float(marked["x"])-float(old["x"]))**2 +
                       (float(marked["y"])-float(old["y"]))**2 < 4.0
                       for old in ifft_spots):
                    continue
                ifft_spots.append(marked)
    else:
        ifft_spots = [dict(point) for point in paired]
    # Reconstruct only the selected reciprocal amplitudes, preserving complex
    # phase and exact registration to the uploaded TEM.
    full_spectrum = np.fft.fftshift(np.fft.fft2(gray - float(np.mean(gray))))
    aperture_mask = np.zeros((height, width), dtype=np.float64)
    for peak in ifft_spots:
        theta = math.radians(float(peak["angle"]))
        cosine, sine = math.cos(theta), math.sin(theta)
        dx, dy = xx - float(peak["x"]), yy - float(peak["y"])
        along = cosine * dx + sine * dy
        across = -sine * dx + cosine * dy
        sigma_x = max(1.65, float(peak["rx"]) / 1.9)
        sigma_y = max(1.25, float(peak["ry"]) / 1.9)
        distance = (along * along / (sigma_x * sigma_x) +
                    across * across / (sigma_y * sigma_y))
        aperture_mask = np.maximum(aperture_mask, np.exp(-.5 * distance))
    reconstruction = np.real(np.fft.ifft2(np.fft.ifftshift(
        full_spectrum * aperture_mask)))
    limit = float(np.percentile(np.abs(reconstruction), 99.5))
    reconstruction = np.clip(
        .5 + .5 * reconstruction / max(limit, 1e-12), 0, 1)
    reconstruction = np.power(reconstruction, .92)
    reconstruction_path = output_dir / "selected_spot_ifft.pgm"
    write_pgm(reconstruction_path, reconstruction * 255.0)

    # The selected shell forms one Square quadrilateral for single-layer data
    # or two nearby quadrilaterals for a bilayer.
    def physical_radius(point):
        return math.hypot(float(point["x"]) - center_x,
                          (float(point["y"]) - center_y) *
                          width / max(float(height), 1.0))

    ordered = sorted(paired, key=physical_radius)
    first_order_periods = []
    first_order_axes = []
    shell_diagnostics = []
    if analysis_kind == "single":
        first_shell = paired[:4]
        for point in first_shell:
            fx = (float(point["x"]) - center_x) / max(float(width), 1.0)
            fy = (float(point["y"]) - center_y) / max(float(height), 1.0)
            frequency = math.hypot(fx, fy)
            if frequency > 0:
                first_order_periods.append(1.0 / frequency)
    else:
        # Build radial shells from indivisible Friedel pairs.  A Square layer
        # needs two period-matched orthogonal axes on one shell; two twisted
        # layers need four such axes on that *same* shell.  In particular, do
        # not combine the (1,0)/(0,1) shell with the sqrt(2)-higher-frequency
        # (1,1)/(1,-1) diagonal shell: that is one lattice rotated by 45 deg,
        # not evidence for a second physical layer.
        axes = []
        for pair_index in range(0, len(paired), 2):
            pair = paired[pair_index:pair_index+2]
            if len(pair) < 2:
                continue
            point = pair[0]
            fx = (float(point["x"])-center_x) / max(float(width), 1.0)
            fy = (float(point["y"])-center_y) / max(float(height), 1.0)
            frequency = math.hypot(fx, fy)
            if frequency <= 0:
                continue
            axes.append({
                "period": 1.0/frequency,
                "angle": math.degrees(math.atan2(-fy, fx)) % 180.0,
                "score": float(point.get("score", 0.0)),
                "points": [dict(value) for value in pair],
            })
        shells = []
        for axis in sorted(axes, key=lambda value: value["period"],
                           reverse=True):
            shell = next((group for group in shells
                          if abs(axis["period"]-group["period"]) /
                          max(group["period"], 1e-9) < .09), None)
            if shell is None:
                shell = {"period": axis["period"], "axes": []}
                shells.append(shell)
            shell["axes"].append(axis)
            shell["period"] = float(np.median(
                [value["period"] for value in shell["axes"]]))
        eligible = []
        for shell in shells:
            unique = []
            for axis in sorted(shell["axes"],
                               key=lambda value: value["score"], reverse=True):
                if any(_circular_distance(axis["angle"], old["angle"],
                                          180.0) < 2.5
                       for old in unique):
                    continue
                unique.append(axis)
            orthogonal = []
            for first_index, first in enumerate(unique):
                for second in unique[first_index+1:]:
                    error = abs(_circular_distance(
                        first["angle"], second["angle"], 180.0)-90.0)
                    mismatch = abs(first["period"]-second["period"]) / max(
                        (first["period"]+second["period"])/2.0, 1e-9)
                    if error <= 7.0 and mismatch <= .10:
                        orthogonal.append((
                            math.sqrt(max(first["score"], 0.0)*
                                      max(second["score"], 0.0)) /
                            (1.0+error/4.0+mismatch*5.0), first, second))
            merit = max((item[0] for item in orthogonal), default=0.0)
            shell_diagnostics.append({
                "period_px": float(shell["period"]),
                "axis_count": len(unique),
                "orthogonal_pair_count": len(orthogonal),
                "merit": float(merit),
            })
            if orthogonal:
                eligible.append((merit, shell["period"], unique, orthogonal))
        if eligible:
            maximum_merit = max(item[0] for item in eligible)
            credible = [item for item in eligible
                        if item[0] >= maximum_merit*.55]
            unused_merit, unused_period, first_order_axes, orthogonal = max(
                credible, key=lambda item: item[1])
            # Four axes are sufficient for two Square layers.  Extra axes are
            # weaker window/noise candidates and must not create more layers.
            first_order_axes = sorted(
                first_order_axes, key=lambda value: value["score"],
                reverse=True)[:4]
            first_shell = [point for axis in first_order_axes
                           for point in axis["points"]]
            first_order_periods = [axis["period"]
                                   for axis in first_order_axes]
        else:
            first_shell = []
    quadrilaterals = []
    measured_twist = None
    if analysis_kind == "single" and len(first_shell) >= 4:
        members = sorted(first_shell, key=lambda point: math.atan2(
            float(point["y"]) - center_y,
            float(point["x"]) - center_x))
        quadrilaterals.append(members)
        if len(single_lattices) > 1:
            for family in single_lattices[1:6]:
                twin_members = sorted(
                    [dict(point) for point in family.get("peaks", [])],
                    key=lambda point: math.atan2(
                        float(point["y"])-center_y,
                        float(point["x"])-center_x))
                if len(twin_members) >= 4:
                    quadrilaterals.append(twin_members[:4])
    elif len(first_order_axes) >= 4:
        phases = [axis["angle"] % 90.0 for axis in first_order_axes]
        clustered = _two_orientation_clusters(
            phases, [axis["score"] for axis in first_order_axes])
        if clustered is not None:
            centers, labels = clustered
            candidate_twist = _circular_distance(
                float(centers[0]), float(centers[1]))
            valid_families = True
            for label_index in (0, 1):
                family_axes = [axis for axis, label_value in zip(
                    first_order_axes, labels)
                               if int(label_value) == label_index]
                if len(family_axes) != 2:
                    valid_families = False
                    break
                orthogonal_error = abs(_circular_distance(
                    family_axes[0]["angle"], family_axes[1]["angle"],
                    180.0)-90.0)
                if orthogonal_error > 7.0:
                    valid_families = False
                    break
                members = [point for axis in family_axes
                           for point in axis["points"]]
                members.sort(key=lambda point: math.atan2(
                    float(point["y"]) - center_y,
                    float(point["x"]) - center_x))
                if len(members) >= 4:
                    quadrilaterals.append(members[:4])
            if (valid_families and len(quadrilaterals) == 2 and
                    .05 < candidate_twist <= 45.0):
                measured_twist = candidate_twist
            else:
                quadrilaterals = []
    if analysis_kind == "bilayer" and not quadrilaterals and \
            len(first_order_axes) >= 2:
        # One complete Square quartet is still useful evidence for lattice a,
        # but it is displayed as one lattice only and never yields twist.
        best_pair = None
        for first_index, first in enumerate(first_order_axes):
            for second in first_order_axes[first_index+1:]:
                error = abs(_circular_distance(
                    first["angle"], second["angle"], 180.0)-90.0)
                if error > 7.0:
                    continue
                merit = math.sqrt(max(first["score"], 0.0)*
                                  max(second["score"], 0.0))/(1.0+error/4.0)
                if best_pair is None or merit > best_pair[0]:
                    best_pair = (merit, first, second)
        if best_pair is not None:
            unused_merit, first, second = best_pair
            members = first["points"]+second["points"]
            members.sort(key=lambda point: math.atan2(
                float(point["y"])-center_y,
                float(point["x"])-center_x))
            quadrilaterals = [members[:4]]
            first_order_periods = [first["period"], second["period"]]
    return {
        "fft_path": str(fft_path),
        "reconstruction_path": str(reconstruction_path),
        "center": [center_x, center_y],
        "selected_spots": ifft_spots,
        "primary_selected_spots": [dict(point) for point in paired],
        "all_visible_spots": all_visible_spots,
        "selected_spot_count": len(ifft_spots),
        "first_order_quadrilaterals": quadrilaterals,
        "first_order_layer_count": len(quadrilaterals),
        "first_order_bilayer_valid": bool(
            analysis_kind == "bilayer" and len(quadrilaterals) == 2 and
            measured_twist is not None),
        "first_order_shells": shell_diagnostics,
        "first_order_lattice_constant_px": (
            float(np.median(first_order_periods))
            if first_order_periods else None),
        "full_fft_twist_angle_deg": measured_twist,
        "single_lattice": dominant_single,
        "single_lattices": single_lattices,
        "source_aspect_ratio": float(width) / max(float(height), 1.0),
    }


def _global_single_layer_lattice_candidates(
        assets, source_shape, retain_combinatorial_alternatives=False):
    """Build Square and sixfold lattice hypotheses from measured FFT peaks.

    Each hypothesis is an actual Friedel-paired set of reciprocal peaks.  The
    local-domain stage decides where it is present; no symmetry is inferred
    from a single line or from an unpaired maximum.  Honeycomb is reserved for
    a near-60/60/60 first-order family.  The user-supplied cryo-EM Kagome
    reference is represented by the experimentally distorted sixfold family,
    typically containing one 65--70 degree reciprocal-axis gap.  Ambiguous
    sixfold families remain explicitly ``Hexagonal``.
    """
    height, width = source_shape
    center_x, center_y = width/2.0, height/2.0
    spots = assets.get("all_visible_spots") or []
    axes = []
    for pair_index in range(0, len(spots)-1, 2):
        point = spots[pair_index]
        fx = (float(point["x"])-center_x)/max(float(width), 1.0)
        fy = (float(point["y"])-center_y)/max(float(height), 1.0)
        frequency = math.hypot(fx, fy)
        if frequency <= 0:
            continue
        axis = {
            "pair_index": pair_index,
            "period": 1.0/frequency,
            "angle": math.degrees(math.atan2(-fy, fx)) % 180.0,
            "score": float(point.get("score", 0.0)),
        }
        if axis["period"] < 3.0 or axis["period"] > min(width, height)*.36:
            continue
        axes.append(axis)
    if not axes:
        return []
    maximum_score = max(axis["score"] for axis in axes)
    axes = [axis for axis in axes
            if axis["score"] >= max(maximum_score*.12, 1e-9)]
    axes.sort(key=lambda value: value["score"], reverse=True)
    # Merge repeated numerical maxima describing the same reciprocal axis.
    # In bilayer mode two genuinely different layers can be separated by only
    # a few degrees.  The old 2.2-degree merge collapsed those two measured
    # axes before a layer pair could be assembled (most visibly in the
    # Honeycomb references).  Retain the resolved alternatives here; the
    # subsequent disjoint-peak and full-axis alignment checks remove noise.
    merge_angle = .55 if retain_combinatorial_alternatives else 2.2
    merge_period = .04 if retain_combinatorial_alternatives else .07
    unique_axes = []
    for axis in axes:
        if any(_circular_distance(axis["angle"], old["angle"], 180.0) <
               merge_angle
               and abs(axis["period"]-old["period"])/max(
                   old["period"], 1e-9) < merge_period
               for old in unique_axes):
            continue
        unique_axes.append(axis)
        if len(unique_axes) >= 64:
            break
    axes = unique_axes

    def periodic_mean(values, period, weights):
        phase = np.asarray(values, dtype=float)*2.0*np.pi/period
        vector = np.sum(np.asarray(weights, dtype=float)*np.exp(1j*phase))
        return float((np.angle(vector)*period/(2.0*np.pi)) % period)

    def peak_members(members):
        output = []
        for member in members:
            start = int(member["pair_index"])
            output.extend(dict(point) for point in spots[start:start+2])
        return output

    raw = []
    for first_index, first in enumerate(axes):
        for second in axes[first_index+1:]:
            angle_error = abs(_circular_distance(
                first["angle"], second["angle"], 180.0)-90.0)
            period_error = abs(first["period"]-second["period"])/max(
                (first["period"]+second["period"])/2.0, 1e-9)
            if angle_error > 6.0 or period_error > .13:
                continue
            mean_period = float(np.median(
                [first["period"], second["period"]]))
            merit = math.sqrt(first["score"]*second["score"])/(
                1.0+angle_error/3.0+period_error*6.0)
            raw.append({
                "valid": True, "symmetry": "Square",
                "orientation_deg": periodic_mean(
                    [first["angle"] % 90.0, second["angle"] % 90.0],
                    90.0, [first["score"], second["score"]]),
                "lattice_constant_px": mean_period,
                "reciprocal_axis_angles_deg": [first["angle"],
                                                 second["angle"]],
                "reciprocal_axis_periods_px": [first["period"],
                                                 second["period"]],
                "axis_angles_deg": [first["angle"], second["angle"]],
                "axis_periods_px": [first["period"], second["period"]],
                "axis_scores": [first["score"], second["score"]],
                "inter_axis_angles_deg": [90.0-angle_error,
                                           90.0+angle_error],
                "axis_pair_indices": [int(first["pair_index"]),
                                      int(second["pair_index"])],
                "peaks": peak_members((first, second)), "merit": merit,
            })

    # Three unoriented reciprocal axes span 180 degrees.  A regular sixfold
    # family has cyclic gaps 60/60/60; Kagome references are allowed the
    # experimentally observed 65--70 degree distorted gap.
    for first_index, first in enumerate(axes):
        for second_index in range(first_index+1, len(axes)):
            second = axes[second_index]
            for third in axes[second_index+1:]:
                members = (first, second, third)
                periods = [member["period"] for member in members]
                period_error = ((max(periods)-min(periods))/max(
                    float(np.mean(periods)), 1e-9))
                if period_error > .16:
                    continue
                angles = sorted(member["angle"] for member in members)
                gaps = [angles[1]-angles[0], angles[2]-angles[1],
                        180.0-angles[2]+angles[0]]
                maximum_gap_error = max(abs(gap-60.0) for gap in gaps)
                # Kagome domains in the cryo-EM references are not an ideal
                # 60/60/60 reciprocal triad: one real-space angle can reach
                # about 65--70 degrees, and finite-domain broadening moves the
                # measured reciprocal gaps farther still.  Keep that measured
                # distortion instead of forcing it into a Square pair or an
                # ideal Honeycomb triad.
                if maximum_gap_error > 15.0 or min(gaps) < 44.5:
                    continue
                if 64.0 <= max(gaps) <= 75.0 and min(gaps) >= 44.5:
                    symmetry = "Kagome"
                elif maximum_gap_error <= 5.0:
                    symmetry = "Honeycomb"
                else:
                    symmetry = "Hexagonal"
                scores = [member["score"] for member in members]
                reciprocal_phase = periodic_mean(
                    [member["angle"] % 60.0 for member in members],
                    60.0, scores)
                reciprocal_period = float(np.median(periods))
                # For triangular/hexagonal direct lattices, the shortest
                # reciprocal period equals sqrt(3)*a/2.
                lattice_constant = 2.0*reciprocal_period/math.sqrt(3.0)
                merit = float(np.prod(scores))**(1.0/3.0)/(
                    1.0+maximum_gap_error/5.0+period_error*6.0)
                raw.append({
                    "valid": True, "symmetry": symmetry,
                    "orientation_deg": (reciprocal_phase-30.0) % 60.0,
                    "lattice_constant_px": lattice_constant,
                    "reciprocal_axis_angles_deg": [member["angle"]
                                                     for member in members],
                    "reciprocal_axis_periods_px": periods,
                    "axis_angles_deg": [member["angle"]
                                         for member in members],
                    "axis_periods_px": periods, "axis_scores": scores,
                    "inter_axis_angles_deg": gaps,
                    "axis_pair_indices": [int(member["pair_index"])
                                          for member in members],
                    "peaks": peak_members(members), "merit": merit,
                })

    if not raw:
        return []
    raw.sort(key=lambda value: value["merit"], reverse=True)
    maximum_merit = raw[0]["merit"]
    candidates = []
    # Do not let numerous accidental orthogonal pairs crowd a valid sixfold
    # hypothesis out of the local-domain competition.  Keep a bounded,
    # symmetry-balanced hypothesis set; the local FFT still has to supply the
    # spatial evidence before any family is reported.
    family_limits = ({"Square": 16, "Honeycomb": 20,
                      "Kagome": 24, "Hexagonal": 8}
                     if retain_combinatorial_alternatives else
                     {"Square": 6, "Honeycomb": 3,
                      "Kagome": 3, "Hexagonal": 2})
    family_counts = {key: 0 for key in family_limits}
    for family in raw:
        if family["merit"] < maximum_merit*.10:
            continue
        symmetry = family["symmetry"]
        if family_counts.get(symmetry, 0) >= family_limits.get(symmetry, 2):
            continue
        period = 90.0 if family["symmetry"] == "Square" else 60.0
        same_members = any(
            old["symmetry"] == family["symmetry"] and
            set(old.get("axis_pair_indices") or []) ==
            set(family.get("axis_pair_indices") or [])
            for old in candidates)
        geometric_duplicate = any(
            old["symmetry"] == family["symmetry"] and
            _circular_distance(old["orientation_deg"],
                               family["orientation_deg"], period) < 4.0 and
            abs(old["lattice_constant_px"]-
                family["lattice_constant_px"])/max(
                    old["lattice_constant_px"], 1e-9) < .10
            for old in candidates)
        # Single-layer local-domain analysis still benefits from a compact
        # geometric hypothesis set.  Bilayer analysis must instead retain
        # alternate three-axis compositions: a mixed A/B triad can have an
        # orientation close to a real layer while sharing only two axes with
        # it.  Collapsing those by angle alone removes the disjoint pair that
        # is required to establish two physical layers.
        if same_members or (geometric_duplicate and
                            not retain_combinatorial_alternatives):
            continue
        family = dict(family)
        family["family_id"] = len(candidates)+1
        candidates.append(family)
        family_counts[symmetry] = family_counts.get(symmetry, 0)+1
        if len(candidates) >= sum(family_limits.values()):
            break
    return candidates


def _axis_fingerprint_alignment(first_angles, second_angles,
                                symmetry_period=60.0):
    """Align two complete unoriented reciprocal-axis families.

    Kagome's measured three-axis family is intentionally not replaced by an
    ideal 60/60/60 model.  Instead, every one-to-one axis correspondence is
    tested under one global rotation.  The residual therefore measures how
    well the *existing measured Kagome angle fingerprint* is preserved.
    """
    first = [float(value) % 180.0 for value in (first_angles or [])]
    second = [float(value) % 180.0 for value in (second_angles or [])]
    if not first or len(first) != len(second):
        return None

    def signed_offset(value, reference, period=180.0):
        return ((value-reference+period/2.0) % period)-period/2.0

    best = None
    for ordered in itertools.permutations(second):
        offsets = [signed_offset(value, reference)
                   for reference, value in zip(first, ordered)]
        phase = np.asarray(offsets, dtype=float)*2.0*np.pi/180.0
        vector = np.sum(np.exp(1j*phase))
        if abs(vector) <= 1e-12:
            continue
        rotation = float(np.angle(vector))*180.0/(2.0*np.pi)
        residuals = [_circular_distance(value, rotation, 180.0)
                     for value in offsets]
        maximum = float(max(residuals))
        rms = float(math.sqrt(np.mean(np.square(residuals))))
        twist = float(_circular_distance(0.0, rotation, symmetry_period))
        aligned = {
            "rotation_deg": float(rotation),
            "twist_angle_deg": twist,
            "maximum_residual_deg": maximum,
            "rms_residual_deg": rms,
            "axis_correspondence_deg": [float(value) for value in ordered],
        }
        rank = (maximum, rms, twist)
        if best is None or rank < best[0]:
            best = (rank, aligned)
    return best[1] if best else None


def _symmetry_aware_bilayer_fft(assets, source_shape,
                                theoretical_a_px=None,
                                theoretical_symmetries=None,
                                preferred_moire_period_px=None,
                                pixel_size_nm=None):
    """Validate two same-order Bravais families without assuming Square.

    Square contributes two reciprocal axes per layer; Honeycomb, Kagome and
    generic hexagonal lattices contribute three.  Two reported layers must
    have the same symmetry, matching ``a``, distinct peak sets and a resolved
    relative orientation.  This keeps the rejection/reporting contract
    identical for every supported symmetry.
    """
    candidates = _global_single_layer_lattice_candidates(
        assets, source_shape, retain_combinatorial_alternatives=True)
    theoretical_a_px = [float(value) for value in
                        (theoretical_a_px or [])
                        if value is not None and float(value) > 0]
    theoretical_symmetries = [str(value) for value in
                              (theoretical_symmetries or [])
                              if str(value) in
                              ("Square", "Honeycomb", "Kagome")]
    if theoretical_symmetries:
        symmetry_matches = [candidate for candidate in candidates
                            if str(candidate.get("symmetry") or "Unknown")
                            in theoretical_symmetries]
        if len(symmetry_matches) >= 2:
            candidates = symmetry_matches

    def theory_error(candidate):
        measured = float(candidate.get("lattice_constant_px") or 0.0)
        if measured <= 0 or not theoretical_a_px:
            return 0.0
        return min(abs(measured-value)/max(value, 1e-9)
                   for value in theoretical_a_px)

    def empirical_a_error(candidate):
        """Soft literature/sample prior used only as a tie breaker."""
        if theoretical_a_px or not pixel_size_nm:
            return 0.0
        symmetry = str(candidate.get("symmetry") or "Unknown")
        expected = {"Kagome": (4.0, 5.0),
                    "Honeycomb": (3.0, 4.0)}.get(symmetry)
        if not expected:
            return 0.0
        measured_nm = float(candidate.get(
            "lattice_constant_px") or 0.0)*float(pixel_size_nm)
        low, high = expected
        if low <= measured_nm <= high:
            return 0.0
        distance = low-measured_nm if measured_nm < low else measured_nm-high
        return distance/((low+high)/2.0)

    # A supplied theoretical a is a reference, not a forced answer.  Apply
    # the filter only when at least two physical layer hypotheses match it;
    # otherwise retain fully automatic analysis and expose that the reference
    # could not be used.
    matching = [candidate for candidate in candidates
                if theory_error(candidate) <= .20]
    theory_filter_applied = bool(
        theoretical_a_px and len(matching) >= 2)
    if theory_filter_applied:
        candidates = matching
    best = None
    for first_index, first in enumerate(candidates):
        symmetry = str(first.get("symmetry") or "Unknown")
        period = 90.0 if symmetry == "Square" else 60.0
        first_peaks = first.get("peaks") or []
        for second in candidates[first_index+1:]:
            if str(second.get("symmetry") or "Unknown") != symmetry:
                continue
            a1 = float(first.get("lattice_constant_px") or 0.0)
            a2 = float(second.get("lattice_constant_px") or 0.0)
            if min(a1, a2) <= 0:
                continue
            a_mismatch = abs(a1-a2)/max((a1+a2)/2.0, 1e-9)
            if a_mismatch > .12:
                continue
            alignment = None
            if symmetry == "Square":
                twist = _circular_distance(
                    float(first.get("orientation_deg") or 0.0),
                    float(second.get("orientation_deg") or 0.0), period)
                fingerprint_error = 0.0
            else:
                alignment = _axis_fingerprint_alignment(
                    first.get("axis_angles_deg"),
                    second.get("axis_angles_deg"), period)
                if not alignment:
                    continue
                # Preserve the established Kagome/Honeycomb orientation
                # definition.  The full-axis alignment validates that the two
                # distorted measured fingerprints are related by one layer
                # rotation; it must not silently redefine Kagome as an ideal
                # 60-degree lattice or replace its calibrated phase angle.
                twist = _circular_distance(
                    float(first.get("orientation_deg") or 0.0),
                    float(second.get("orientation_deg") or 0.0), period)
                fingerprint_error = float(
                    alignment["maximum_residual_deg"])
            # Below 0.35 degree the peak families are not independently
            # resolvable at the detector's angular tolerance.
            if not .35 < twist <= period/2.0:
                continue
            if symmetry != "Square":
                first_gaps = sorted(float(value) for value in
                                    (first.get("inter_axis_angles_deg") or []))
                second_gaps = sorted(float(value) for value in
                                     (second.get("inter_axis_angles_deg") or []))
                if len(first_gaps) != 3 or len(second_gaps) != 3:
                    continue
                # Kagome is deliberately not forced to 60 degrees.  Its two
                # layers must preserve the same measured distorted three-axis
                # fingerprint under an overall rotation.  Honeycomb uses the
                # same comparison and naturally remains near 60/60/60.
                fingerprint_tolerance = (8.0 if symmetry == "Kagome"
                                         else 4.0)
                if fingerprint_error > fingerprint_tolerance:
                    continue
            second_peaks = second.get("peaks") or []
            overlap = 0
            for point in first_peaks:
                if any((float(point["x"])-float(other["x"]))**2 +
                       (float(point["y"])-float(other["y"]))**2 < 1.0
                       for other in second_peaks):
                    overlap += 1
            if overlap > 0:
                continue
            denominator = math.sqrt(max(
                1e-12, a1*a1+a2*a2-2.0*a1*a2*
                math.cos(math.radians(twist))))
            predicted_period = a1*a2/denominator
            merit = (math.sqrt(
                max(float(first.get("merit", 0.0)), 0.0) *
                max(float(second.get("merit", 0.0)), 0.0)) /
                (1.0+a_mismatch*6.0+fingerprint_error/30.0))
            if theory_filter_applied:
                merit /= (1.0+5.0*(theory_error(first)+
                                    theory_error(second)))
            elif pixel_size_nm:
                # Kagome ~4--5 nm and Honeycomb ~3--4 nm are approximate,
                # never exclusion ranges.  They only break ties between
                # otherwise valid peak families; user-entered theoretical a
                # remains the higher-priority reference.
                merit /= (1.0+1.5*(empirical_a_error(first)+
                                    empirical_a_error(second)))
            if (preferred_moire_period_px and
                    float(preferred_moire_period_px) > 0):
                period_error = abs(predicted_period-
                                   float(preferred_moire_period_px))/float(
                                       preferred_moire_period_px)
                # Independent TEM repetition is corroborating evidence, not
                # a forced answer.  It only ranks otherwise valid, disjoint
                # same-symmetry FFT pairs.
                merit /= (1.0+12.0*period_error)
            if best is None or merit > best[0]:
                best = (merit, first, second, twist, fingerprint_error,
                        alignment, predicted_period)
    if best is None:
        return {"valid": False,
                "error": "Equal-order FFT peaks could not be separated reliably into two lattices of the same symmetry.",
                "candidates": candidates,
                "theoretical_a_filter_applied": theory_filter_applied,
                "theoretical_symmetries": theoretical_symmetries,
                "empirical_a_prior_applied": bool(
                    pixel_size_nm and not theoretical_a_px)}
    (unused_merit, first, second, twist, fingerprint_error,
     axis_alignment, predicted_period) = best
    selected_ids = {id(first), id(second)}
    used_peaks = list(first.get("peaks") or []) + list(
        second.get("peaks") or [])
    minimum_pair_merit = min(float(first.get("merit", 0.0)),
                             float(second.get("merit", 0.0)))
    additional_layers = []
    rejected_harmonic_layers = []
    primary_a = (float(first.get("lattice_constant_px") or 0.0) +
                 float(second.get("lattice_constant_px") or 0.0)) / 2.0
    for candidate in sorted(
            candidates, key=lambda value: float(value.get("merit", 0.0)),
            reverse=True):
        if id(candidate) in selected_ids:
            continue
        if float(candidate.get("merit", 0.0)) < minimum_pair_merit*.45:
            continue
        candidate_peaks = candidate.get("peaks") or []
        if any((float(point["x"])-float(other["x"]))**2 +
               (float(point["y"])-float(other["y"]))**2 < 16.0
               for point in candidate_peaks for other in used_peaks):
            continue
        candidate_a = float(candidate.get("lattice_constant_px") or 0.0)
        same_symmetry = (str(candidate.get("symmetry") or "Unknown") ==
                         str(first.get("symmetry") or "Unknown"))
        if same_symmetry and min(candidate_a, primary_a) > 0:
            ratio = max(candidate_a, primary_a) / min(candidate_a, primary_a)
            # A Square (h,h) diagonal shell has an apparent spacing related
            # to the first-order shell by sqrt(2).  It uses disjoint FFT spots
            # and therefore passed the old "unused peaks" test, falsely
            # turning an ordinary bilayer into a three-layer model.  Higher
            # integer-index shells are handled the same way.  A genuine third
            # layer with comparable a (ratio near 1) is not rejected.
            harmonic = next((value for value in
                             (math.sqrt(2.0), 2.0, math.sqrt(5.0))
                             if abs(ratio-value)/value <= .08), None)
            if harmonic is not None:
                rejected_harmonic_layers.append(dict(
                    candidate, harmonic_ratio=float(harmonic),
                    measured_ratio=float(ratio)))
                continue
        additional_layers.append(dict(candidate))
        used_peaks.extend(candidate_peaks)
        if len(additional_layers) >= 3:
            break
    a1 = float(first["lattice_constant_px"])
    a2 = float(second["lattice_constant_px"])
    ignored_additional_layers = list(additional_layers)
    additional_layers = []
    return {
        "valid": True,
        "symmetry": str(first.get("symmetry") or "Unknown"),
        "layers": [dict(first, layer_role="twist_pair_1"),
                   dict(second, layer_role="twist_pair_2")] + [
                       dict(layer, layer_role="additional_layer")
                       for layer in additional_layers],
        "layer_count": 2,
        "mixed_multilayer": False,
        "primary_twist_symmetry": str(
            first.get("symmetry") or "Unknown"),
        "primary_twist_layer_indices": [0, 1],
        "additional_layers": additional_layers,
        "rejected_harmonic_layers": rejected_harmonic_layers,
        "ignored_additional_layers": ignored_additional_layers,
        "theoretical_a_filter_applied": theory_filter_applied,
        "theoretical_a_px": theoretical_a_px,
        "theoretical_symmetries": theoretical_symmetries,
        "empirical_a_prior_applied": bool(
            pixel_size_nm and not theoretical_a_px and
            str(first.get("symmetry") or "Unknown") in
            ("Kagome", "Honeycomb")),
        "lattice_constant_px": float((a1+a2)/2.0),
        "twist_angle_deg": float(twist),
        "axis_fingerprint_error_deg": float(fingerprint_error),
        "axis_alignment": axis_alignment,
        "predicted_moire_period_px": float(predicted_period),
        "tem_period_reference_px": (
            float(preferred_moire_period_px)
            if preferred_moire_period_px else None),
        "candidates": candidates,
    }


def _direct_real_space_lattice(image, expected_px, scale_bar=None):
    """Estimate the short lattice spacing directly from TEM autocorrelation.

    This is deliberately independent of the reciprocal-space distance used
    for the FFT value reported in the bulk CSV.  The FFT estimate is only used
    to define a broad, physically plausible annulus; the final distance comes
    from real-space correlation peaks.
    """
    if not expected_px or expected_px < 3:
        return {"valid": False,
                "error": "A reliable FFT lattice-scale prior is unavailable."}
    height, width = image.shape
    crop_height = (max(32, int(scale_bar["y0"] - 4))
                   if scale_bar else int(height * 0.92))
    crop = np.asarray(image[:crop_height], dtype=np.float64)
    stride = max(1, int(math.ceil(max(crop.shape) / 900.0)))
    sample = crop[::stride, ::stride]
    expected = float(expected_px) / stride
    sample -= _gaussian_filter_fft(sample, max(2.0, expected * 1.5))
    sample -= float(np.mean(sample))
    window = np.outer(np.hanning(sample.shape[0]), np.hanning(sample.shape[1]))
    transformed = np.fft.fft2(sample * window)
    correlation = np.real(np.fft.fftshift(np.fft.ifft2(
        transformed * np.conjugate(transformed))))
    cy, cx = np.array(correlation.shape) // 2
    yy, xx = np.mgrid[:correlation.shape[0], :correlation.shape[1]]
    radial = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    annulus = (radial >= expected * 0.70) & (radial <= expected * 1.30)
    candidates = np.argwhere(annulus)
    ranked = sorted(
        ((float(correlation[y, x]), int(x), int(y), float(radial[y, x]))
         for y, x in candidates), reverse=True)
    chosen = []
    for score, x, y, distance in ranked:
        angle = math.degrees(math.atan2(y - cy, x - cx)) % 180.0
        if any(abs(((angle - old["angle_deg"] + 90.0) % 180.0) - 90.0) < 14.0
               for old in chosen):
            continue
        chosen.append({"distance_px": distance * stride,
                       "angle_deg": angle, "score": score})
        if len(chosen) >= 2:
            break
    if not chosen:
        return {"valid": False,
                "error": "No lattice peaks were found in the TEM real-space autocorrelation."}
    return {
        "valid": True,
        "lattice_constant_px": float(np.median(
            [item["distance_px"] for item in chosen])),
        "directions": chosen,
    }


def ocr_scale(original, scale_bar=None):
    if not original or not Path(original).is_file():
        return None, None
    tesseract = tool_executable("tesseract")
    if not tesseract:
        return None, None
    environment = os.environ.copy()
    tessdata = Path(tesseract).resolve().parent / "tessdata"
    if tessdata.is_dir():
        environment["TESSDATA_PREFIX"] = str(tessdata)
    try:
        result = subprocess.run(
            [tesseract, original, "stdout", "--psm", "11",
             "tsv"],
            check=False, text=True, capture_output=True, timeout=20,
            env=environment)
    except Exception:
        return None, None
    words = []
    for line in result.stdout.splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) < 12 or fields[0] != "5":
            continue
        try:
            x, y, width, height = map(int, fields[6:10])
            confidence = float(fields[10])
        except ValueError:
            continue
        word = fields[11].strip()
        if confidence < 0 or not word:
            continue
        if scale_bar:
            margin_x = max(40, int(scale_bar["pixel_length"] * 1.5))
            near_x = (x + width >= scale_bar["x0"] - margin_x and
                      x <= scale_bar["x1"] + margin_x)
            near_y = (y + height >= scale_bar["y0"] - 80 and
                      y <= scale_bar["y1"] + 28)
            if not (near_x and near_y):
                continue
        words.append((y, x, word))
    text = " ".join(item[2] for item in sorted(words))
    matches = re.findall(
        r"(?<![\d.])(\d+(?:[.,]\d+)?)\s*(nm|um|µm|μm)\b",
        text, flags=re.IGNORECASE)
    if not matches:
        return None, text.strip()
    value, unit = matches[-1]
    value = float(value.replace(",", "."))
    if unit.lower() != "nm":
        value *= 1000.0
    return value, text.strip()


def _circular_mean_90(values, weights=None):
    values = np.asarray(values, dtype=float) % 90.0
    weights = (np.ones(len(values), dtype=float) if weights is None else
               np.asarray(weights, dtype=float))
    phase = np.deg2rad(values * 4.0)
    vector = np.sum(weights * np.exp(1j * phase))
    return float((np.rad2deg(np.angle(vector)) / 4.0) % 90.0)


def _welch_lattice_fft(image, scale_bar=None, coarse_period_px=None):
    """Resolve the two first-order Square Bragg families from a TEM image.

    Median/mean patch spectra suppress the specimen envelope and the much
    lower-frequency moire modulation.  The closest-to-centre rule is applied
    only after sharp peaks have formed two reciprocal, orthogonal Square
    families; central background and moire satellites are therefore not
    eligible lattice peaks.
    """
    height, width = image.shape
    crop_height = (max(160, int(scale_bar["y0"] - 8))
                   if scale_bar else int(height * 0.94))
    patch = min(384, int(min(width, crop_height) * 0.46))
    patch = max(160, patch - patch % 8)
    padded = patch * 4
    step = max(48, patch // 4)
    margin_x = max(16, int(width * 0.055))
    margin_y = max(16, int(crop_height * 0.07))
    spectra = np.zeros((padded, padded), dtype=np.float64)
    accepted = 0
    window = np.outer(np.hanning(patch), np.hanning(patch))
    for y in range(margin_y, crop_height - patch + 1, step):
        for x in range(margin_x, width - patch + 1, step):
            tile = image[y:y + patch, x:x + patch]
            if tile.shape != (patch, patch) or float(np.std(tile)) < 4.0:
                continue
            tile = (tile - np.mean(tile)) * window
            spectrum = np.abs(np.fft.fftshift(
                np.fft.fft2(tile, s=(padded, padded))))
            spectrum /= float(np.median(spectrum)) + 1e-9
            spectra += np.log1p(spectrum)
            accepted += 1
    if accepted < 3:
        return {"valid": False,
                "error": "Too few valid local FFT regions were obtained from the TEM image."}
    spectra /= accepted
    center = padded // 2
    candidates = []
    # DNA-origami/SST lattice periods below four pixels are unresolved; very
    # long periods belong to the moire/envelope branch and are excluded here.
    maximum_offset = int(padded / 4.0)
    for y in range(center - maximum_offset, center + maximum_offset + 1):
        for x in range(center - maximum_offset, center + maximum_offset + 1):
            dx, dy = x - center, y - center
            radius = math.hypot(dx, dy)
            if radius <= 0:
                continue
            period = padded / radius
            maximum_lattice_period = min(80.0, patch * 0.22)
            if coarse_period_px:
                maximum_lattice_period = min(
                    maximum_lattice_period, float(coarse_period_px) * 0.45)
            if not 4.0 <= period <= maximum_lattice_period:
                continue
            value = float(spectra[y, x])
            if value < float(np.max(spectra[y - 2:y + 3,
                                            x - 2:x + 3])):
                continue
            angle = math.degrees(math.atan2(-dy, dx)) % 180.0
            candidates.append({
                "value": value, "period": period, "angle": angle,
                "x": x, "y": y, "radius": radius,
            })
    if not candidates:
        return {"valid": False,
                "error": "No first-order lattice peaks were found in the FFT."}
    maximum = max(item["value"] for item in candidates)
    strong = [item for item in candidates
              if item["value"] >= maximum * 0.72]
    # Collapse the exact Friedel partner while retaining the two nearby
    # orientations (one for each physical layer).
    collapsed = []
    for item in sorted(strong, key=lambda value: value["value"], reverse=True):
        if any(_circular_distance(item["angle"], old["angle"], 180.0) < 0.7
               for old in collapsed):
            continue
        collapsed.append(item)
    # Separate radial orders before orientation clustering. The first-order
    # Square shell is the validated shell closest to the FFT centre, i.e. the
    # shell with the largest real-space period. Without this step the strong
    # diagonal/second-order shell can masquerade as a layer rotated by 45°.
    shells = []
    for item in sorted(collapsed, key=lambda value: value["period"],
                       reverse=True):
        shell = next((group for group in shells
                      if abs(item["period"] - group["period"]) /
                      max(group["period"], 1e-9) < 0.08), None)
        if shell is None:
            shell = {"period": item["period"], "items": []}
            shells.append(shell)
        shell["items"].append(item)
        shell["period"] = float(np.median(
            [member["period"] for member in shell["items"]]))
    eligible_shells = [shell for shell in shells if len(shell["items"]) >= 4]
    if eligible_shells:
        collapsed = max(eligible_shells,
                        key=lambda shell: shell["period"])["items"]
    # A two-layer Square first-order shell contributes four unique axes after
    # Friedel collapse. Lower peaks on the same radius are window side-lobes.
    if len(collapsed) > 4:
        collapsed = sorted(collapsed,
                           key=lambda item: item["value"], reverse=True)[:4]
    if len(collapsed) < 4:
        return {"valid": False,
                "error": "First-order FFT peaks could not be separated "
                         "stably into two Square lattices."}
    phases = np.asarray([item["angle"] % 90.0 for item in collapsed])
    # The expected twist is the smallest non-zero split of two Square
    # orientation families. Search every pair of phase seeds and keep the
    # partition with two orthogonal observations per layer and least scatter.
    best = None
    for first in range(len(phases)):
        for second in range(first + 1, len(phases)):
            if _circular_distance(phases[first], phases[second]) < 0.35:
                continue
            centers = np.asarray([phases[first], phases[second]], dtype=float)
            for unused in range(12):
                distances = np.asarray([
                    [_circular_distance(value, center) for center in centers]
                    for value in phases])
                labels = np.argmin(distances, axis=1)
                updated = []
                for label in (0, 1):
                    selection = labels == label
                    if np.count_nonzero(selection) < 2:
                        break
                    updated.append(_circular_mean_90(
                        phases[selection],
                        [collapsed[index]["value"]
                         for index in np.flatnonzero(selection)]))
                if len(updated) != 2:
                    break
                centers = np.asarray(updated)
            else:
                residual = sum(
                    _circular_distance(phases[index], centers[labels[index]]) ** 2
                    for index in range(len(phases)))
                balance = abs(np.count_nonzero(labels == 0) -
                              np.count_nonzero(labels == 1))
                score = residual + balance * 0.5
                if best is None or score < best[0]:
                    best = (score, centers, labels)
    if best is None:
        return {"valid": False,
                "error": "Orientation clustering of the first-order FFT peaks failed."}
    unused_score, centers, labels = best
    maximum_phase_residual = max(
        _circular_distance(phases[index], centers[labels[index]])
        for index in range(len(phases)))
    if maximum_phase_residual > 3.0:
        return {"valid": False,
                "error": "The FFT candidate peaks could not form two orthogonal sets of first-order Square peaks.",
                "candidate_peaks": collapsed}
    layers = []
    for label in (0, 1):
        members = [collapsed[index] for index in range(len(collapsed))
                   if labels[index] == label]
        periods = [item["period"] for item in members]
        layers.append({
            "orientation_deg": float(centers[label]),
            "lattice_constant_px": float(np.median(periods)),
            "peaks": members,
        })
    twist = _circular_distance(layers[0]["orientation_deg"],
                               layers[1]["orientation_deg"])
    if twist <= 0.05 or twist > 45.0:
        return {"valid": False,
                "error": "The angle between the two FFT-derived Square lattices is unreliable."}
    a1 = layers[0]["lattice_constant_px"]
    a2 = layers[1]["lattice_constant_px"]
    denominator = math.sqrt(max(
        1e-12, a1 * a1 + a2 * a2 -
        2.0 * a1 * a2 * math.cos(math.radians(twist))))
    predicted_period = a1 * a2 / denominator
    return {
        "valid": True,
        "patch_count": accepted,
        "fft_size": padded,
        "layers": layers,
        "lattice_constant_px": float((a1 + a2) / 2.0),
        "twist_angle_deg": float(twist),
        "predicted_moire_period_px": float(predicted_period),
    }


def _long_runs(mask):
    padded = np.pad(mask.astype(np.int8), (1, 1))
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    return list(zip(starts.tolist(), ends.tolist()))


def detect_scale_bar(image):
    height, width = image.shape
    y_start = int(height * 0.55)
    roi = image[y_start:]
    low_percentile, high_percentile = np.percentile(roi, [3, 97])
    median = float(np.median(roi))
    low = min(low_percentile,
              median - 0.65 * (median - float(np.min(roi))))
    high = max(high_percentile,
               median + 0.65 * (float(np.max(roi)) - median))
    candidates = []
    for polarity, mask in (("bright", roi >= high), ("dark", roi <= low)):
        runs_by_row = []
        for local_y, row in enumerate(mask):
            runs = [(x0, x1) for x0, x1 in _long_runs(row)
            if width * 0.025 <= x1 - x0 <= width * 0.72]
            for x0, x1 in runs:
                runs_by_row.append((local_y + y_start, x0, x1))
        for y, x0, x1 in runs_by_row:
            thickness = 1
            for yy, xx0, xx1 in runs_by_row:
                if yy <= y or yy > y + max(12, int(height * 0.03)):
                    continue
                if abs(xx0 - x0) <= 3 and abs(xx1 - x1) <= 3:
                    thickness = max(thickness, yy - y + 1)
            length = x1 - x0
            # Scale bars are long, thin rectangles.  Dark particles or image
            # contamination can contain similarly long runs, but their filled
            # thickness is much larger; reject those compact blobs before
            # scoring polarity candidates.
            if thickness > max(12, int(round(length * 0.22))):
                continue
            bottom_weight = 1.0 + 0.6 * (y / max(height, 1))
            inside = image[y:min(height, y + thickness), x0:x1]
            surround_y0 = max(0, y - max(3, thickness * 2))
            surround = image[surround_y0:y, x0:x1]
            contrast = (abs(float(np.mean(inside)) - float(np.mean(surround)))
                        if surround.size else 0.0)
            if contrast < max(4.0, 0.08 * (float(np.max(roi)) -
                                           float(np.min(roi)))):
                continue
            score = length * min(thickness, 10) * bottom_weight * contrast
            candidates.append((score, x0, y, x1, y + thickness - 1,
                               polarity))
    if not candidates:
        return None
    best = max(candidates)
    unused_score, x0, y0, x1, y1, polarity = best
    # The percentile threshold locates the solid plateau of an antialiased
    # scale bar.  Expand to its visible edges so a 10-nm cryo-EM bar is not
    # shortened by several pixels (which would overestimate a by 5--8%).
    inner = image[y0:y1+1, x0:x1]
    background_parts = []
    if y0 > 0:
        background_parts.append(image[max(0, y0-8):y0,
                                      max(0, x0-4):min(width, x1+4)])
    if y1+1 < height:
        background_parts.append(image[y1+1:min(height, y1+9),
                                      max(0, x0-4):min(width, x1+4)])
    if inner.size and any(part.size for part in background_parts):
        core_value = float(np.median(inner))
        background_value = float(np.median(np.concatenate(
            [part.ravel() for part in background_parts if part.size])))
        edge_value = background_value + .22*(core_value-background_value)

        def belongs(column):
            values = image[y0:y1+1, column]
            # Rounded/slightly tilted bar ends may occupy only part of the
            # detected thickness.  A robust upper/lower quartile retains those
            # visible end pixels without following isolated image noise.
            value = float(np.percentile(
                values, 75.0 if polarity == "bright" else 25.0))
            return (value >= edge_value if polarity == "bright" else
                    value <= edge_value)

        limit = max(8, int(round((x1-x0)*.25)))
        steps = 0
        while (x0 > 0 and steps < limit and
               (belongs(x0-1) or (x0 > 1 and belongs(x0-2)))):
            x0 -= 1
            steps += 1
        steps = 0
        while (x1 < width and steps < limit and
               (belongs(x1) or (x1+1 < width and belongs(x1+1)))):
            x1 += 1
            steps += 1
    return {
        "x0": int(x0), "y0": int(y0), "x1": int(x1), "y1": int(y1),
        "pixel_length": float(x1 - x0), "polarity": polarity,
    }


def _peak_candidates(values, center, maximum=120, minimum_radius=5,
                     maximum_radius=None, nms=7):
    height, width = values.shape
    cy, cx = center
    yy, xx = np.ogrid[:height, :width]
    radius = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    if maximum_radius is None:
        maximum_radius = min(width, height) * 0.47
    masked = values.copy()
    masked[(radius < minimum_radius) | (radius > maximum_radius)] = -np.inf
    order = np.argpartition(masked.ravel(), -min(masked.size, maximum * 30))[
        -min(masked.size, maximum * 30):]
    order = order[np.argsort(masked.ravel()[order])[::-1]]
    peaks = []
    for flat in order:
        y, x = np.unravel_index(int(flat), masked.shape)
        if not np.isfinite(masked[y, x]):
            continue
        if any((x - px) ** 2 + (y - py) ** 2 < nms ** 2
               for px, py, unused in peaks):
            continue
        peaks.append((int(x), int(y), float(masked[y, x])))
        if len(peaks) >= maximum:
            break
    return peaks


def _period_groups(peaks, shape):
    height, width = shape
    cx, cy = width / 2.0, height / 2.0
    points = []
    for x, y, weight in peaks:
        fx = (x - cx) / width
        fy = (y - cy) / height
        frequency = math.hypot(fx, fy)
        if frequency > 0:
            points.append((1.0 / frequency, weight, x, y))
    points.sort(reverse=True)
    groups = []
    for point in points:
        period = point[0]
        target = next((group for group in groups
                       if abs(period - group["period"]) /
                       max(group["period"], 1e-9) < 0.10), None)
        if target is None:
            target = {"period": period, "points": []}
            groups.append(target)
        target["points"].append(point)
        target["period"] = float(np.median(
            [item[0] for item in target["points"]]))
    for group in groups:
        group["score"] = (len(group["points"]) ** 0.5 *
                          float(np.mean([item[1]
                                         for item in group["points"]])))
    return [group for group in groups if len(group["points"]) >= 2]


def _secondary_pore_lattice(groups, fft_shape, stride, source_shape,
                            helix_lattice):
    """Find an optional pore/supercell lattice sharing the helix orientation.

    A valid pore lattice must provide Friedel pairs along both axes of the
    already established single-layer Square helix lattice.  Its spacing must
    be a near-integer multiple of the helix spacing; this rejects arbitrary
    low-frequency texture, contamination and unrelated domain envelopes.
    """
    if not helix_lattice or not helix_lattice.get("valid"):
        return {"valid": False}
    targets = helix_lattice.get("axis_angles_deg") or []
    helix_period = float(helix_lattice.get("lattice_constant_px") or 0.0)
    if len(targets) < 2 or helix_period <= 0:
        return {"valid": False}
    fft_height, fft_width = fft_shape
    source_height, source_width = source_shape
    cx, cy = fft_width / 2.0, fft_height / 2.0
    candidates = []
    for group in groups:
        matched_axes = []
        for target_angle in targets[:2]:
            members = []
            for period, weight, x, y in group.get("points", []):
                angle = math.degrees(math.atan2(
                    -(y-cy) / max(float(fft_height), 1.0),
                    (x-cx) / max(float(fft_width), 1.0))) % 180.0
                if _circular_distance(angle, float(target_angle), 180.0) <= 6.0:
                    members.append({
                        "period": float(period) * stride,
                        "weight": float(weight), "x": float(x), "y": float(y),
                        "angle": float(angle)})
            if len(members) < 2:
                matched_axes = []
                break
            axis_period = float(np.median(
                [member["period"] for member in members]))
            matched_axes.append({"period": axis_period, "members": members})
        if len(matched_axes) != 2:
            continue
        periods = [axis["period"] for axis in matched_axes]
        mean_period = float(sum(periods) / 2.0)
        mismatch = abs(periods[0] - periods[1]) / max(mean_period, 1e-9)
        ratio = mean_period / helix_period
        repeat = int(round(ratio))
        if mismatch > .12 or repeat < 3 or repeat > 16:
            continue
        if abs(ratio - repeat) > .16:
            continue
        peak_points = []
        for axis in matched_axes:
            selected = sorted(
                axis["members"],
                key=lambda member: (abs(member["period"]-axis["period"]),
                                    -member["weight"]))[:2]
            for member in selected:
                peak_points.append({
                    "x": float(source_width / 2.0 +
                               (member["x"]-cx) * source_width /
                               max(float(fft_width * stride), 1.0)),
                    "y": float(source_height / 2.0 +
                               (member["y"]-cy) * source_height /
                               max(float(fft_height * stride), 1.0)),
                    "score": member["weight"]})
        peak_points.sort(key=lambda point: math.atan2(
            point["y"]-source_height/2.0,
            point["x"]-source_width/2.0))
        candidates.append({
            "valid": True,
            "lattice_constant_px": mean_period,
            "axis_periods_px": periods,
            "ratio_to_helix": ratio,
            "repeat_multiple": repeat,
            "orientation_deg": float(helix_lattice.get(
                "orientation_deg", targets[0] % 90.0)),
            "axis_angles_deg": [float(value) for value in targets[:2]],
            "period_mismatch_fraction": mismatch,
            "peaks": peak_points})
    if not candidates:
        return {"valid": False}
    return max(candidates, key=lambda item: item["lattice_constant_px"])


def _rectangular_sum(integral, height, width):
    return (integral[height:, width:] - integral[:-height, width:] -
            integral[height:, :-width] + integral[:-height, :-width])


def _moire_real_space_units(image, predicted_period_px, scale_bar=None):
    """Locate repeated TEM units near an FFT-predicted period.

    A normalized template correlation supplies candidate unit centres. The
    final period is the median of neighbour distances along two approximately
    orthogonal axes, so the FFT prediction constrains the search without
    forcing the TEM measurement to equal it.
    """
    if not predicted_period_px or not math.isfinite(predicted_period_px):
        return {"valid": False,
                "error": "No FFT-predicted Moire period is available."}
    height, width = image.shape
    crop_height = (max(160, int(scale_bar["y0"] - 8))
                   if scale_bar else int(height * 0.94))
    period = float(predicted_period_px)
    if period < 24 or period > min(width, crop_height) * 0.42:
        return {"valid": False,
                "error": "The FFT-predicted Moire period exceeds the TEM field of view."}
    size = int(round(period * 0.80))
    size = max(48, min(220, size + size % 2))
    center_x = width // 2
    center_y = int(crop_height * 0.445)
    half = size // 2
    center_x = min(max(center_x, half), width - half - 1)
    center_y = min(max(center_y, half), crop_height - half - 1)
    template = image[center_y - half:center_y + half,
                     center_x - half:center_x + half].copy()
    if template.shape != (size, size):
        return {"valid": False,
                "error": "The central TEM template extends beyond the image boundary."}
    template = ((template - np.mean(template)) *
                np.outer(np.hanning(size), np.hanning(size)))
    template -= np.mean(template)
    padded_shape = (height + size - 1, width + size - 1)
    fft_shape = (1 << (padded_shape[0] - 1).bit_length(),
                 1 << (padded_shape[1] - 1).bit_length())
    correlation = np.fft.irfft2(
        np.fft.rfft2(image, fft_shape) *
        np.fft.rfft2(template[::-1, ::-1], fft_shape), fft_shape)
    correlation = correlation[:padded_shape[0], :padded_shape[1]]
    correlation = correlation[size - 1:height, size - 1:width]
    integral = np.pad(np.cumsum(np.cumsum(image, axis=0), axis=1),
                      ((1, 0), (1, 0)))
    integral_sq = np.pad(np.cumsum(np.cumsum(image * image, axis=0), axis=1),
                         ((1, 0), (1, 0)))
    local_sum = _rectangular_sum(integral, size, size)
    local_sq = _rectangular_sum(integral_sq, size, size)
    local_variance = np.maximum(
        local_sq - local_sum * local_sum / float(size * size), 1e-9)
    normalized = correlation / np.sqrt(
        local_variance * (float(np.sum(template * template)) + 1e-9))
    working = normalized.copy()
    yy, xx = np.indices(working.shape)
    edge = max(18, int(period * 0.22))
    valid = ((xx > edge) & (xx < working.shape[1] - edge) &
             (yy > edge) & (yy < min(working.shape[0], crop_height - size) - edge))
    working[~valid] = -np.inf
    candidates = []
    suppression = max(18, int(period * 0.68))
    for unused in range(180):
        y, x = np.unravel_index(int(np.argmax(working)), working.shape)
        score = float(working[y, x])
        if not math.isfinite(score) or score < 0.18:
            break
        candidates.append({"x": float(x + half), "y": float(y + half),
                           "score": score})
        working[max(0, y - suppression):y + suppression + 1,
                max(0, x - suppression):x + suppression + 1] = -np.inf
    if len(candidates) < 6:
        return {"valid": False,
                "error": "Too few repeating-unit centers were found in the TEM image.",
                "centers": candidates}
    pairs = []
    for first in range(len(candidates)):
        for second in range(first + 1, len(candidates)):
            dx = candidates[second]["x"] - candidates[first]["x"]
            dy = candidates[second]["y"] - candidates[first]["y"]
            distance = math.hypot(dx, dy)
            if not 0.72 * period <= distance <= 1.28 * period:
                continue
            angle = math.degrees(math.atan2(dy, dx)) % 180.0
            weight = math.exp(-0.5 * ((distance - period) /
                                      max(period * 0.18, 1.0)) ** 2)
            pairs.append({"first": first, "second": second, "dx": dx,
                          "dy": dy, "distance": distance,
                          "angle": angle, "weight": weight})
    if len(pairs) < 4:
        return {"valid": False,
                "error": "The TEM unit centers do not form a two-dimensional lattice.",
                "centers": candidates}
    # Direction histogram followed by an approximately orthogonal partner.
    histogram = []
    for angle in np.arange(0.0, 180.0, 2.0):
        score = sum(pair["weight"] for pair in pairs
                    if _circular_distance(pair["angle"], angle, 180.0) < 6.0)
        histogram.append((score, float(angle)))
    first_angle = max(histogram)[1]
    orthogonal = (first_angle + 90.0) % 180.0
    second_options = [item for item in histogram
                      if _circular_distance(item[1], orthogonal, 180.0) < 20.0]
    second_angle = max(second_options)[1] if second_options else orthogonal
    selected = [pair for pair in pairs
                if min(_circular_distance(pair["angle"], first_angle, 180.0),
                       _circular_distance(pair["angle"], second_angle, 180.0)) < 9.0]
    if len(selected) < 4:
        selected = pairs
    distances = np.asarray([pair["distance"] for pair in selected])
    observed = float(np.median(distances))
    q1, q3 = np.percentile(distances, [25, 75])

    def basis_for(angle):
        members = [pair for pair in selected
                   if _circular_distance(pair["angle"], angle, 180.0) < 9.0]
        if not members:
            radians = math.radians(angle)
            return [observed * math.cos(radians),
                    observed * math.sin(radians)]
        vectors = []
        radians = math.radians(angle)
        unit_x, unit_y = math.cos(radians), math.sin(radians)
        for pair in members:
            dx, dy = pair["dx"], pair["dy"]
            if dx * unit_x + dy * unit_y < 0:
                dx, dy = -dx, -dy
            vectors.append((dx, dy))
        return [float(np.median([item[0] for item in vectors])),
                float(np.median([item[1] for item in vectors]))]

    basis_1 = basis_for(first_angle)
    basis_2 = basis_for(second_angle)
    neighbour_counts = [0] * len(candidates)
    for pair in selected:
        neighbour_counts[pair["first"]] += 1
        neighbour_counts[pair["second"]] += 1
    centers = [dict(item, neighbour_count=neighbour_counts[index])
               for index, item in enumerate(candidates)
               if neighbour_counts[index] > 0]
    representative = min(
        selected,
        key=lambda pair: (abs(pair["distance"] - observed) +
                          0.002 * math.hypot(
                              (candidates[pair["first"]]["x"] +
                               candidates[pair["second"]]["x"]) / 2 - width / 2,
                              (candidates[pair["first"]]["y"] +
                               candidates[pair["second"]]["y"]) / 2 -
                              crop_height / 2)))
    return {
        "valid": True,
        "centers": centers,
        "cell_count": len(centers),
        "basis_vectors_px": [basis_1, basis_2],
        "moire_period_px": observed,
        "period_iqr_px": [float(q1), float(q3)],
        "representative_pair": {
            "first": candidates[representative["first"]],
            "second": candidates[representative["second"]],
            "distance_px": representative["distance"],
        },
    }


def _independent_moire_period_seed(image, period_candidates_px,
                                   lattice_constant_px=None,
                                   scale_bar=None):
    """Find a real-space Moiré seed without requiring an accepted twist.

    The coarse low-frequency FFT groups only propose search scales.  A scale
    is retained only when normalized real-space matching produces a repeated
    two-dimensional set with at least two units across the short image axis.
    Trying the longest credible scale first avoids selecting its harmonic.
    """
    height, width = image.shape
    crop_height = (max(160, int(scale_bar["y0"] - 8))
                   if scale_bar else int(height * 0.94))
    minimum = max(24.0, 3.2*float(lattice_constant_px or 0.0))
    maximum = min(width, crop_height)*.42
    seeds = []
    for value in sorted((float(item) for item in
                         (period_candidates_px or [])
                         if item is not None), reverse=True):
        if not minimum <= value <= maximum:
            continue
        if any(abs(value-old)/max(old, 1e-9) < .08 for old in seeds):
            continue
        seeds.append(value)
        if len(seeds) >= 5:
            break
    for seed in seeds:
        measured = _moire_real_space_units(image, seed, scale_bar)
        autocorrelation = None
        if not measured.get("valid"):
            # A centre template can miss a Moiré cell when the image centre
            # lies on a defect or low-contrast part of the modulation.  Use a
            # global autocorrelation only to refine the search scale, then run
            # the same repeated-unit validation again at that refined scale.
            autocorrelation = _direct_real_space_lattice(
                image, seed, scale_bar)
            refined_seed = float(autocorrelation.get(
                "lattice_constant_px") or 0.0)
            if refined_seed > 0:
                measured = _moire_real_space_units(
                    image, refined_seed, scale_bar)
        if not measured.get("valid"):
            continue
        observed = float(measured.get("moire_period_px") or 0.0)
        if observed <= 0 or min(width, crop_height)/observed < 2.0:
            continue
        if int(measured.get("cell_count", 0)) < 4:
            continue
        q1, q3 = measured.get("period_iqr_px") or (observed, observed)
        if (float(q3)-float(q1))/observed > .35:
            continue
        measured = dict(measured)
        measured["independent_seed_px"] = float(seed)
        if autocorrelation and autocorrelation.get("valid"):
            measured["autocorrelation_refined_seed_px"] = float(
                autocorrelation.get("lattice_constant_px"))
        measured["method"] = "coarse_fft_seed_validated_in_real_space"
        return measured
    return {"valid": False,
            "error": "The TEM real-space image does not contain independently verifiable repeating Moire units."}


def analyze_tem(image, pixel_size_nm=None, scale_bar=None,
                analysis_kind="bilayer", output_dir=None,
                manual_edits=None, theoretical_a_nm=None,
                theoretical_symmetries=None):
    height, width = image.shape
    theoretical_a_nm = [float(value) for value in
                        (theoretical_a_nm or [])
                        if value is not None and float(value) > 0]
    theoretical_a_px = ([value/float(pixel_size_nm)
                         for value in theoretical_a_nm]
                        if pixel_size_nm else [])
    theoretical_symmetries = [str(value) for value in
                              (theoretical_symmetries or [])
                              if str(value) in
                              ("Square", "Honeycomb", "Kagome")]
    crop_height = (max(32, int(scale_bar["y0"] - 4))
                   if scale_bar else int(height * 0.90))
    crop = image[:crop_height]
    # Bound FFT cost while retaining real-space period in original pixels.
    stride = max(1, int(math.ceil(max(crop.shape) / 1024.0)))
    sample = crop[::stride, ::stride]
    sample = sample - np.mean(sample)
    window = np.outer(np.hanning(sample.shape[0]), np.hanning(sample.shape[1]))
    spectrum = np.log1p(np.abs(np.fft.fftshift(np.fft.fft2(sample * window))))
    center = (spectrum.shape[0] / 2.0, spectrum.shape[1] / 2.0)
    peaks = _peak_candidates(spectrum, center, maximum=100)
    groups = _period_groups(peaks, spectrum.shape)
    for group in groups:
        group["period"] *= stride
    plausible = [group for group in groups
                 if 4 <= group["period"] <= min(width, crop_height) * 0.45]
    strongest = sorted(plausible, key=lambda group: group["score"], reverse=True)[:8]
    strongest.sort(key=lambda group: group["period"], reverse=True)
    moire_px = strongest[0]["period"] if strongest else None
    lattice_px = strongest[-1]["period"] if len(strongest) > 1 else None
    # Avoid selecting a high harmonic as lattice: prefer a strong group whose
    # period is at least one quarter of the moire period.
    if moire_px and len(strongest) > 1:
        lattice_candidates = [group for group in strongest[1:]
                              if group["period"] >= moire_px * 0.24]
        if lattice_candidates:
            lattice_px = max(lattice_candidates,
                             key=lambda group: group["score"])["period"]
    precise_fft = (_precise_fft_assets(
        image, output_dir, analysis_kind=analysis_kind)
                   if output_dir else {})
    bilayer_seed_real_space = {"valid": False}
    single_lattice = precise_fft.get("single_lattice") or {}
    pore_lattice = (_secondary_pore_lattice(
        groups, spectrum.shape, stride, image.shape, single_lattice)
                    if analysis_kind == "single" else {"valid": False})
    # A global low-frequency shell in a twin image can be produced by the
    # domain envelope/interference of the two orientations.  It must not be
    # imposed as one shared supercell basis.  Until a supercell is independently
    # validated for each orientation, index the two helix bases separately.
    if len(precise_fft.get("single_lattices") or []) > 1:
        pore_lattice = {
            "valid": False,
            "reason": "twin_orientations_require_separate_supercell_bases",
        }
    if (analysis_kind == "single" and pore_lattice.get("valid") and
            output_dir and precise_fft.get("selected_spots")):
        # A two-a single-layer image needs both reciprocal shells in the IFFT.
        # Refine the low-frequency pore/supercell quartet on the full FFT,
        # combine it with the helix-scale quartet, and preserve the complex
        # phase so the reconstructed domains remain registered to the TEM.
        secondary_spots = _refine_fft_spot_estimates(
            image, pore_lattice.get("peaks") or [])
        for point in secondary_spots:
            point.update({
                "lattice_role": "square_supercell_first_order",
                "basis_family": "square_supercell",
                "reflection_class": "axis",
                "radial_order": 1.0,
            })
        indexed_spots = _indexed_supercell_spots(
            precise_fft, pore_lattice, image.shape)
        combined_spots = (indexed_spots or [dict(point) for point in
                          precise_fft.get("primary_selected_spots", [])])
        for point in secondary_spots:
            if any((float(point["x"])-float(old["x"]))**2 +
                   (float(point["y"])-float(old["y"]))**2 < 16.0
                   for old in combined_spots):
                continue
            combined_spots.append(point)
        if len(combined_spots) > len(precise_fft.get("selected_spots", [])):
            combined_path = Path(output_dir) / "selected_spot_ifft_two_a.pgm"
            _write_selected_spot_ifft(image, combined_spots, combined_path)
            precise_fft.setdefault("primary_selected_spots", [dict(point)
                for point in (single_lattice.get("peaks") or [])])
            precise_fft["secondary_selected_spots"] = secondary_spots
            precise_fft["indexed_reflection_spots"] = indexed_spots
            precise_fft["indexed_reflection_pair_count"] = (
                len(indexed_spots) // 2)
            precise_fft["selected_spots"] = combined_spots
            precise_fft["selected_spot_count"] = len(combined_spots)
            precise_fft["reconstruction_path"] = str(combined_path)
            precise_fft["lattice_relation"] = "helix_supercell"
            precise_fft["selection_method"] = (
                "clear_square_integer_indexed_reflections")
            precise_fft["reflection_selection"] = (
                _reflection_selection_summary(combined_spots))
            pore_lattice["refined_peaks"] = secondary_spots
            pore_lattice["relation"] = "supercell"
    if analysis_kind == "single" and single_lattice.get("valid"):
        families = precise_fft.get("single_lattices") or [single_lattice]
        layers = [{
            "orientation_deg": family["orientation_deg"],
            "lattice_constant_px": family["lattice_constant_px"],
            "peaks": family.get("peaks") or [],
            "lattice_role": ("primary" if index == 0 else
                             "twin_orientation"),
        } for index, family in enumerate(families[:6])]
        lattice_fft = {
            "valid": True,
            "method": ("single_square_multi_orientation_domains"
                       if len(layers) > 2 else
                       "single_square_twin_orientations"
                       if len(layers) > 1 else
                       "dominant_single_square_helix_first_order"),
            "measurement_scale": "helix_lattice",
            "layers": layers,
            "lattice_constant_px": single_lattice["lattice_constant_px"],
        }
        if len(layers) > 1:
            orientation_pairs = []
            for first_index in range(len(layers)):
                for second_index in range(first_index+1, len(layers)):
                    orientation_pairs.append({
                        "orientation_1_index": first_index,
                        "orientation_2_index": second_index,
                        "relative_orientation_deg": _circular_distance(
                            layers[first_index]["orientation_deg"],
                            layers[second_index]["orientation_deg"], 90.0),
                    })
            lattice_fft["twin"] = {
                "valid": True,
                "orientations_deg": [layer["orientation_deg"]
                                     for layer in layers],
                "relative_orientation_deg": _circular_distance(
                    layers[0]["orientation_deg"],
                    layers[1]["orientation_deg"], 90.0),
                "lattice_constant_px": float(np.mean(
                    [layer["lattice_constant_px"] for layer in layers])),
                "orientation_pairs": orientation_pairs,
            }
    else:
        lattice_fft = _welch_lattice_fft(image, scale_bar, moire_px)
        symmetry_bilayer = (_symmetry_aware_bilayer_fft(
            precise_fft, image.shape, theoretical_a_px,
            theoretical_symmetries, pixel_size_nm=pixel_size_nm)
                            if precise_fft else {"valid": False})
        if symmetry_bilayer.get("valid"):
            initial_period = symmetry_bilayer.get(
                "predicted_moire_period_px")
            initial_real_space = _moire_real_space_units(
                image, initial_period, scale_bar)
            if not initial_real_space.get("valid"):
                bilayer_seed_real_space = _independent_moire_period_seed(
                    image, [group.get("period") for group in strongest],
                    symmetry_bilayer.get("lattice_constant_px"), scale_bar)
                if bilayer_seed_real_space.get("valid"):
                    reranked = _symmetry_aware_bilayer_fft(
                        precise_fft, image.shape, theoretical_a_px,
                        theoretical_symmetries,
                        bilayer_seed_real_space.get("moire_period_px"),
                        pixel_size_nm)
                    if reranked.get("valid"):
                        reranked_period = float(reranked.get(
                            "predicted_moire_period_px") or 0.0)
                        measured_period = float(bilayer_seed_real_space.get(
                            "moire_period_px") or 0.0)
                        if (reranked_period > 0 and measured_period > 0 and
                                abs(reranked_period-measured_period)/
                                measured_period <= .35):
                            symmetry_bilayer = reranked
                            symmetry_bilayer[
                                "tem_period_reranking_applied"] = True
        precise_fft["bilayer_symmetry_candidates"] = (
            symmetry_bilayer.get("candidates") or [])
        if symmetry_bilayer.get("valid"):
            # Bilayer mode has an explicit two-layer ceiling.  The
            # symmetry-aware validator supplies exactly the two accepted
            # first-order peak families for display.  Welch remains the
            # higher-resolution Square estimator when it is valid; otherwise
            # the symmetry-aware result is the numerical fallback.
            if (symmetry_bilayer.get("symmetry") != "Square" or
                    not lattice_fft.get("valid")):
                lattice_fft = dict(symmetry_bilayer)
            else:
                lattice_fft = dict(lattice_fft)
                lattice_fft["layers"] = [dict(layer) for layer in
                                         symmetry_bilayer.get("layers", [])[:2]]
                lattice_fft["symmetry"] = "Square"
                lattice_fft["mixed_multilayer"] = False
                lattice_fft["detected_layer_count"] = 2
                lattice_fft["theoretical_a_filter_applied"] = (
                    symmetry_bilayer.get("theoretical_a_filter_applied", False))
            polygons = [list(layer.get("peaks") or [])
                        for layer in symmetry_bilayer.get("layers", [])[:2]]
            precise_fft.update({
                "first_order_quadrilaterals": polygons,
                "first_order_layer_count": 2,
                "first_order_bilayer_valid": True,
                "first_order_lattice_constant_px": symmetry_bilayer.get(
                    "lattice_constant_px"),
                "full_fft_twist_angle_deg": symmetry_bilayer.get(
                    "twist_angle_deg"),
                "bilayer_symmetry": symmetry_bilayer.get("symmetry"),
                "mixed_multilayer": False,
                "detected_layer_count": 2,
                "theoretical_a_filter_applied": symmetry_bilayer.get(
                    "theoretical_a_filter_applied", False),
            })

    # When no lower-frequency supercell basis supersedes the helix basis,
    # select every clearly resolved Square reflection that is integer-indexed
    # by a validated orientation.  This includes the 45-degree (h,h) shell and
    # mixed (h,k) orders, while rejecting unrelated maxima.  Single-image and
    # bulk analysis call this same worker, so they use identical selection.
    if (output_dir and precise_fft and
            not (analysis_kind == "single" and pore_lattice.get("valid"))):
        basis_families = []
        if analysis_kind == "single":
            source_families = (precise_fft.get("single_lattices") or
                               ([single_lattice]
                                if single_lattice.get("valid") else []))
            for index, family in enumerate(source_families):
                basis_families.append({
                    "family_id": "square_orientation_%d" % (index + 1),
                    "lattice_role": ("square_helix_reflection" if index == 0
                                     else "square_twin_reflection"),
                    "axis_angles_deg": family.get("axis_angles_deg") or [],
                    "axis_periods_px": family.get("axis_periods_px") or [],
                    "lattice_constant_px": family.get(
                        "lattice_constant_px"),
                    "residual_tolerance": .20,
                })
        elif lattice_fft.get("valid"):
            for index, layer in enumerate(lattice_fft.get("layers") or []):
                orientation = float(layer.get("orientation_deg", 0.0))
                period = float(layer.get("lattice_constant_px") or 0.0)
                symmetry = str(layer.get("symmetry") or
                               lattice_fft.get("symmetry") or "Square")
                axis_angles = (layer.get("reciprocal_axis_angles_deg") or
                               layer.get("axis_angles_deg") or
                               [orientation, orientation + 90.0])
                axis_periods = (layer.get("reciprocal_axis_periods_px") or
                                layer.get("axis_periods_px") or
                                [period for unused in axis_angles])
                basis_families.append({
                    "family_id": "square_layer_%d" % (index + 1),
                    "symmetry": symmetry,
                    "lattice_role": ("square_layer_reflection" if index == 0
                                     else "square_layer_2_reflection"),
                    "axis_angles_deg": axis_angles,
                    "axis_periods_px": axis_periods,
                    "residual_tolerance": .20,
                })
        indexed_spots = _indexed_square_spots(
            precise_fft, basis_families, image.shape,
            residual_tolerance=.20, score_fraction=.42)
        if len(indexed_spots) >= 4:
            indexed_path = (Path(output_dir) /
                            "selected_spot_ifft_square_indexed.pgm")
            _write_selected_spot_ifft(image, indexed_spots, indexed_path)
            precise_fft["indexed_reflection_spots"] = indexed_spots
            precise_fft["indexed_reflection_pair_count"] = (
                len(indexed_spots) // 2)
            precise_fft["selected_spots"] = indexed_spots
            precise_fft["selected_spot_count"] = len(indexed_spots)
            precise_fft["reconstruction_path"] = str(indexed_path)
            precise_fft["selection_method"] = (
                "clear_square_integer_indexed_reflections")
            precise_fft["reflection_selection"] = (
                _reflection_selection_summary(indexed_spots))
    orientation_domains = {"valid": False}
    if analysis_kind == "single" and precise_fft:
        mixed_candidates = _global_single_layer_lattice_candidates(
            precise_fft, image.shape)
        candidate_domains = (_local_square_orientation_domains(
            image, mixed_candidates, scale_bar, manual_edits)
                             if mixed_candidates else
            {"valid": False})
        # Activate heterogeneous-domain mode only when a non-Square family or
        # a genuinely different a occupies a meaningful connected area.  Pure
        # Square/twin images therefore retain the established path unchanged.
        family_area = {}
        family_score = {}
        manually_retained_indices = set()
        for domain in candidate_domains.get("domains") or []:
            index = int(domain.get("orientation_index", -1))
            family_area[index] = family_area.get(index, 0.0) + float(
                domain.get("area_fraction", 0.0))
            family_score.setdefault(index, []).append(float(
                domain.get("mean_local_score", 0.0)))
            if domain.get("manual_protected"):
                manually_retained_indices.add(index)
        significant = []
        for index, family in enumerate(mixed_candidates):
            area = family_area.get(index, 0.0)
            score = float(np.mean(family_score.get(index, [0.0])))
            if ((area >= .055 and score >= .20) or
                    index in manually_retained_indices):
                significant.append((index, family, area, score))
        non_square = [item for item in significant
                      if item[1].get("symmetry") != "Square"]
        different_a = (not pore_lattice.get("valid") and any(
            abs(first[1]["lattice_constant_px"]-
                second[1]["lattice_constant_px"])/max(
                    first[1]["lattice_constant_px"], 1e-9) >= .14
            for first_index, first in enumerate(significant)
            for second in significant[first_index+1:]))
        use_mixed = bool(candidate_domains.get("valid") and
                         (non_square or different_a or manual_edits))
        if use_mixed:
            orientation_domains = candidate_domains
            retained_indices = [item[0] for item in significant]
            retained_layers = []
            for index in retained_indices:
                family = mixed_candidates[index]
                retained_layers.append({
                    "symmetry": family.get("symmetry", "Square"),
                    "orientation_deg": family.get("orientation_deg"),
                    "lattice_constant_px": family.get(
                        "lattice_constant_px"),
                    "reciprocal_axis_angles_deg": family.get(
                        "reciprocal_axis_angles_deg") or [],
                    "reciprocal_axis_periods_px": family.get(
                        "reciprocal_axis_periods_px") or [],
                    "inter_axis_angles_deg": family.get(
                        "inter_axis_angles_deg") or [],
                    "peaks": family.get("peaks") or [],
                    "domain_area_fraction": family_area.get(index, 0.0),
                })
            total_area = sum(layer["domain_area_fraction"]
                             for layer in retained_layers) or 1.0
            lattice_fft.update({
                "valid": True,
                "method": "single_heterogeneous_symmetry_a_domains",
                "measurement_scale": "domain_specific_lattice",
                "layers": retained_layers,
                "lattice_constant_px": float(sum(
                    layer["lattice_constant_px"]*
                    layer["domain_area_fraction"]
                    for layer in retained_layers)/total_area),
                "heterogeneous_domains": True,
                "twin": {"valid": False,
                          "reason": "heterogeneous_symmetry_or_a"},
            })
            pore_lattice = {
                "valid": False,
                "reason": "domain_specific_lattices_not_global_supercell",
            }
            precise_fft["mixed_lattice_families"] = retained_layers
            domain_spots = []
            for layer in retained_layers:
                for point in layer.get("peaks") or []:
                    if any((float(point["x"])-float(old["x"]))**2 +
                           (float(point["y"])-float(old["y"]))**2 < 4.0
                           for old in domain_spots):
                        continue
                    marked = dict(point)
                    marked["lattice_role"] = "domain_first_order"
                    marked["basis_family"] = layer["symmetry"]
                    domain_spots.append(marked)
            # Apply the same integer reciprocal-index validation used for
            # Square diagonal/mixed-order peaks to every measured Bravais
            # basis.  Honeycomb and Kagome therefore retain their own legal
            # oblique higher orders rather than borrowing Square's 45° rule.
            domain_bases = []
            for index, layer in enumerate(retained_layers):
                domain_bases.append({
                    "family_id": "domain_%d_%s" % (
                        index+1, str(layer.get("symmetry") or "lattice")),
                    "symmetry": str(layer.get("symmetry") or "Square"),
                    "lattice_role": "%s_domain_reflection" % str(
                        layer.get("symmetry") or "lattice").lower(),
                    "axis_angles_deg": layer.get(
                        "reciprocal_axis_angles_deg") or [],
                    "axis_periods_px": layer.get(
                        "reciprocal_axis_periods_px") or [],
                    "lattice_constant_px": layer.get(
                        "lattice_constant_px"),
                    "residual_tolerance": .19,
                })
            indexed_domain_spots = _indexed_square_spots(
                precise_fft, domain_bases, image.shape,
                residual_tolerance=.19, score_fraction=.42)
            if len(indexed_domain_spots) >= len(domain_spots):
                domain_spots = indexed_domain_spots
            if output_dir and domain_spots:
                path = Path(output_dir) / "selected_spot_ifft_domains.pgm"
                _write_selected_spot_ifft(image, domain_spots, path)
                precise_fft["selected_spots"] = domain_spots
                precise_fft["selected_spot_count"] = len(domain_spots)
                precise_fft["reconstruction_path"] = str(path)
                precise_fft["selection_method"] = (
                    "domain_specific_integer_indexed_reciprocal_reflections")
                precise_fft["reflection_selection"] = (
                    _reflection_selection_summary(domain_spots))
        else:
            orientation_domains = _local_square_orientation_domains(
                image, precise_fft.get("single_lattices") or
                [single_lattice], scale_bar, manual_edits)
        lattice_fft["orientation_domains"] = orientation_domains
    if analysis_kind == "single":
        # Sample/background segmentation is independent of whether the global
        # FFT produced a complete lattice model.  This matters for cropped
        # sheets with an open specimen edge: they still need a visible sample
        # boundary even when the first-order peaks are incomplete.
        current_boundary = (orientation_domains.get("sample_boundary") or
                            {})
        manual_background_edit = any(str(edit.get("type") or "").startswith(
            "background_") for edit in (manual_edits or []))
        if (not current_boundary.get("valid") and
                not manual_background_edit and
                not orientation_domains.get(
                    "automatic_background_rejected")):
            independent_boundary = _radial_sample_boundary(image, scale_bar)
            if not independent_boundary.get("valid"):
                independent_boundary = _periodicity_sample_boundary(
                    image, scale_bar)
            if independent_boundary.get("valid"):
                orientation_domains = dict(orientation_domains)
                orientation_domains["sample_boundary"] = independent_boundary
                orientation_domains["background_excluded"] = True
        lattice_fft["orientation_domains"] = orientation_domains
    elif analysis_kind == "bilayer":
        # Bilayer analysis is a global lattice/FFT/Moire measurement.  Sample
        # and orientation-domain segmentation belongs to the single-layer
        # workflow and must not affect either the result or its rendering.
        orientation_domains = {
            "valid": False,
            "disabled": True,
            "reason": "bilayer_global_lattice_analysis",
            "domains": [],
            "boundaries": [],
            "sample_boundary": {"valid": False, "disabled": True},
            "background_excluded": False,
        }
        lattice_fft["orientation_domains"] = orientation_domains
    if precise_fft.get("first_order_lattice_constant_px"):
        lattice_px = precise_fft["first_order_lattice_constant_px"]
    fft_twist_reliable = bool(
        analysis_kind == "bilayer" and
        precise_fft.get("first_order_bilayer_valid") and
        lattice_fft.get("valid") and
        lattice_fft.get("twist_angle_deg") is not None)
    fft_twist_reason = (
        "" if fft_twist_reliable else
        "The two sets of equal-order twist peaks in the FFT are too close, or the evidence is insufficient for reliable separation.")
    if analysis_kind == "bilayer" and not fft_twist_reliable:
        # Do not leak a rejected candidate angle/period through the detailed
        # JSON or UI fallbacks.  Retain only the shell-derived lattice a.
        candidate_families = (precise_fft.get(
            "bilayer_symmetry_candidates") or [])
        best_family = max(candidate_families,
                          key=lambda value: float(value.get("merit", 0.0)),
                          default={})
        measured_a = (precise_fft.get("first_order_lattice_constant_px") or
                      best_family.get("lattice_constant_px") or
                      lattice_fft.get("lattice_constant_px") or lattice_px)
        measured_symmetry = str(best_family.get("symmetry") or
                                lattice_fft.get("symmetry") or "Unknown")
        lattice_fft = dict(lattice_fft)
        lattice_fft.update({
            "valid": False,
            "error": fft_twist_reason,
            "symmetry": measured_symmetry,
            "layers": ([{"symmetry": measured_symmetry,
                         "lattice_constant_px": float(measured_a)}]
                       if measured_a else []),
            "lattice_constant_px": (float(measured_a)
                                     if measured_a else None),
        })
        lattice_fft.pop("twist_angle_deg", None)
        lattice_fft.pop("predicted_moire_period_px", None)
    mixed_multilayer = bool(lattice_fft.get("mixed_multilayer"))
    real_space = (_moire_real_space_units(
        image, lattice_fft.get("predicted_moire_period_px"), scale_bar)
        if fft_twist_reliable and not mixed_multilayer else
        {"valid": False,
         "error": (
             "Moire period is not calculated for single-layer analysis."
             if analysis_kind == "single" else
             "The TEM period in a mixed multilayer image cannot be assigned uniquely to one layer pair."
             if mixed_multilayer else fft_twist_reason)})
    if lattice_fft.get("valid"):
        lattice_px = lattice_fft["lattice_constant_px"]
    fft_lattice_px = (precise_fft.get("first_order_lattice_constant_px") or
                      lattice_px)
    tem_lattice = _direct_real_space_lattice(
        image, fft_lattice_px, scale_bar)
    tem_lattice_px = tem_lattice.get("lattice_constant_px")
    if real_space.get("valid"):
        moire_px = real_space["moire_period_px"]
    elif analysis_kind == "bilayer":
        moire_px = None
    moire_units_in_short_axis = (
        float(min(width, crop_height)) / float(moire_px)
        if moire_px else 0.0)
    tem_period_reliable = bool(
        fft_twist_reliable and real_space.get("valid") and
        moire_units_in_short_axis >= 2.0 and
        int(real_space.get("cell_count", 0)) >= 4)
    if analysis_kind == "single":
        tem_period_reason = (
            "Moire period is not calculated for single-layer analysis.")
    elif mixed_multilayer:
        tem_period_reason = (
            "A mixed multilayer image was detected; the TEM period cannot be "
            "assigned uniquely to one layer pair. Only reliable FFT-derived "
            "twist values for same-symmetry layer pairs are reported.")
    elif not fft_twist_reliable:
        tem_period_reason = (
            "The TEM field of view contains fewer than two reliable Moire "
            "units, so the TEM twist cannot be determined.")
    elif not real_space.get("valid"):
        tem_period_reason = real_space.get(
            "error", "No stable Moire period was identified in TEM real space.")
    elif moire_units_in_short_axis < 2.0:
        tem_period_reason = (
            "The TEM field of view spans fewer than two Moire units along "
            "its short axis, so TEM twist and period cannot be determined "
            "reliably.")
    elif int(real_space.get("cell_count", 0)) < 4:
        tem_period_reason = (
            "Too few usable Moire units are present in the TEM image for a "
            "reliable period fit.")
    else:
        tem_period_reason = ""
    result = {
        "analysis_kind": analysis_kind,
        "theoretical_a_nm": theoretical_a_nm,
        "theoretical_symmetries": theoretical_symmetries,
        "theoretical_a_filter_applied": bool(
            precise_fft.get("theoretical_a_filter_applied")),
        "manual_domain_edits": list(manual_edits or []),
        "moire_period_px": moire_px,
        "lattice_constant_px": lattice_px,
        "fft_peak_count": len(peaks),
        "period_candidates_px": [round(group["period"], 4)
                                 for group in strongest],
        "lattice_fft": lattice_fft,
        "tem_lattice_real_space": tem_lattice,
        "tem_lattice_constant_px": tem_lattice_px,
        "fft_lattice_constant_px": fft_lattice_px,
        "moire_real_space": real_space,
        "tem_period_reliable": tem_period_reliable,
        "tem_period_reliability_reason": tem_period_reason,
        "fft_twist_reliable": fft_twist_reliable,
        "fft_twist_reliability_reason": fft_twist_reason,
        "mixed_multilayer": mixed_multilayer,
        "detected_layer_count": int(
            lattice_fft.get("layer_count") or
            len(lattice_fft.get("layers") or [])),
        "moire_units_in_short_axis": moire_units_in_short_axis,
        "fft_assets": precise_fft,
        "pore_lattice": pore_lattice,
        "orientation_domains": orientation_domains,
    }
    if analysis_kind == "single":
        result["moire_period_px"] = None
        result["single_layer_orientations_deg"] = [
            float(layer.get("orientation_deg"))
            for layer in lattice_fft.get("layers", [])
            if layer.get("orientation_deg") is not None]
    if fft_twist_reliable:
        predicted = lattice_fft["predicted_moire_period_px"]
        result["fft_twist_angle_deg"] = lattice_fft["twist_angle_deg"]
        result["fft_predicted_moire_period_px"] = predicted
        if real_space.get("valid"):
            result["moire_consistency_percent"] = float(
                100.0 * (real_space["moire_period_px"] - predicted) /
                predicted)
    if pixel_size_nm:
        result["moire_period_nm"] = (result.get("moire_period_px") * pixel_size_nm
                                     if result.get("moire_period_px") else None)
        result["lattice_constant_nm"] = (lattice_px * pixel_size_nm
                                          if lattice_px else None)
        result["tem_lattice_constant_nm"] = (
            tem_lattice_px * pixel_size_nm if tem_lattice_px else None)
        result["fft_lattice_constant_nm"] = (
            fft_lattice_px * pixel_size_nm if fft_lattice_px else None)
        if pore_lattice.get("valid"):
            result["pore_lattice_constant_nm"] = (
                pore_lattice["lattice_constant_px"] * pixel_size_nm)
    return result


def _circular_distance(a, b, period=90.0):
    distance = abs(a - b) % period
    return min(distance, period - distance)


def _two_orientation_clusters(angles, weights):
    if len(angles) < 4:
        return None
    angles = np.asarray(angles, dtype=float) % 90.0
    weights = np.asarray(weights, dtype=float)
    # Initialize from the strongest separated pair, then circular k-means.
    first = int(np.argmax(weights))
    separated = np.asarray([_circular_distance(value, angles[first])
                            for value in angles])
    second = int(np.argmax(separated * (weights / max(weights.max(), 1e-9))))
    centers = np.asarray([angles[first], angles[second]], dtype=float)
    for unused in range(30):
        distances = np.asarray([
            [_circular_distance(value, center) for center in centers]
            for value in angles])
        labels = np.argmin(distances, axis=1)
        updated = []
        for label in (0, 1):
            selected = labels == label
            if not np.any(selected):
                return None
            phase = np.deg2rad(angles[selected] * 4.0)
            vector = np.sum(weights[selected] * np.exp(1j * phase))
            updated.append((np.rad2deg(np.angle(vector)) / 4.0) % 90.0)
        updated = np.asarray(updated)
        if max(_circular_distance(updated[i], centers[i])
               for i in (0, 1)) < 1e-4:
            centers = updated
            break
        centers = updated
    return centers, labels


def analyze_fft(image):
    height, width = image.shape
    center = (height / 2.0, width / 2.0)
    contrast = np.maximum(image - np.median(image),
                          np.median(image) - image)
    peaks = _peak_candidates(contrast, center, maximum=80,
                             minimum_radius=max(6, min(width, height) * 0.03),
                             maximum_radius=min(width, height) * 0.45,
                             nms=max(5, int(min(width, height) * 0.012)))
    cx, cy = width / 2.0, height / 2.0
    # Use the dominant radial shell to avoid mixing higher-order reflections.
    radii = np.asarray([math.hypot(x - cx, y - cy) for x, y, unused in peaks])
    if len(radii) < 4:
        return {"twist_angle_deg": None, "peaks": [],
                "error": "Too few lattice peaks were identified in the FFT."}
    bins = max(12, int(min(width, height) / 18))
    histogram, edges = np.histogram(radii, bins=bins,
                                    weights=[p[2] for p in peaks])
    shell = int(np.argmax(histogram))
    low, high = edges[shell], edges[shell + 1]
    selected = [(x, y, weight) for (x, y, weight), radius in zip(peaks, radii)
                if low <= radius <= high]
    if len(selected) < 4:
        selected = peaks[:min(24, len(peaks))]
    angles = [math.degrees(math.atan2(-(y - cy), x - cx)) % 90.0
              for x, y, unused in selected]
    weights = [max(weight, 1e-9) for unused, unused_y, weight in selected]
    clustered = _two_orientation_clusters(angles, weights)
    if clustered is None:
        return {"twist_angle_deg": None, "peaks": [],
                "error": "FFT lattice orientations could not be separated stably into two groups."}
    centers, labels = clustered
    twist = _circular_distance(float(centers[0]), float(centers[1]))
    marked = [
        {"x": int(point[0]), "y": int(point[1]), "layer": int(label)}
        for point, label in zip(selected, labels)]
    return {
        "twist_angle_deg": float(twist),
        "orientation_1_deg": float(centers[0]),
        "orientation_2_deg": float(centers[1]),
        "center": [cx, cy],
        "peaks": marked,
        "peak_count": len(marked),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("tem", "fft"))
    parser.add_argument("pgm")
    parser.add_argument("--original")
    parser.add_argument("--analysis-kind", choices=("single", "bilayer"),
                        default="bilayer")
    parser.add_argument("--output-dir")
    parser.add_argument("--domain-edits")
    parser.add_argument("--theoretical-a-nm", action="append", type=float,
                        default=[])
    parser.add_argument("--theoretical-symmetry", action="append",
                        choices=("Square", "Honeycomb", "Kagome"),
                        default=[])
    parser.add_argument("--pixel-size-nm", type=float)
    parser.add_argument("--scale-value-nm", type=float)
    args = parser.parse_args()
    image = read_pgm(args.pgm)
    if args.mode == "tem":
        manual_edits = []
        if args.domain_edits:
            manual_edits = json.loads(Path(args.domain_edits).read_text(
                encoding="utf-8"))
        bar = detect_scale_bar(image)
        scale_nm, ocr_text = ocr_scale(args.original, bar)
        if args.scale_value_nm and args.scale_value_nm > 0:
            scale_nm = float(args.scale_value_nm)
        pixel_size = (float(args.pixel_size_nm)
                      if args.pixel_size_nm and args.pixel_size_nm > 0 else
                      scale_nm / bar["pixel_length"]
                      if bar and scale_nm else None)
        result = {
            "mode": "tem", "width": image.shape[1], "height": image.shape[0],
            "scale_bar": bar, "scale_value_nm": scale_nm,
            "ocr_text": ocr_text, "pixel_size_nm": pixel_size,
        }
        result.update(analyze_tem(
            image, pixel_size, bar, args.analysis_kind, args.output_dir,
            manual_edits, args.theoretical_a_nm,
            args.theoretical_symmetry))
    else:
        result = {"mode": "fft", "width": image.shape[1],
                  "height": image.shape[0]}
        result.update(analyze_fft(image))
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()

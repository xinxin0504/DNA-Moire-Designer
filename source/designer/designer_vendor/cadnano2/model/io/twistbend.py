"""Planning helpers for the interactive Twist and Bend editor.

The module deliberately contains no Qt code.  It converts a set of helix
regions and deformation parameters into integer insertion/deletion edits.  A
view may preview the returned geometry before applying those edits to cadnano.
"""

from __future__ import division

import math

from .indelanalysis import optimize_bend_plan_curvature
from .missinghelixcalibration import (
    NODES as MISSING_HELIX_NODES,
    S8_R4X4C_BRANCH_RMSE,
    S8_R4X4C_CALIBRATION_VERSION,
    S8_R4X4C_EXTENDED_PROFILE,
    S8_R4X4C_FAILED_INDEL_POINTS,
    S8_R4X4C_TWIST_POINTS,
    S8_R4X4C_VALIDATED_INDEL_RANGE,
    VALIDATION_RMSE as MISSING_HELIX_RMSE)


RISE_NM = 0.34
TARGET_PITCH = 10.44
CENTER_SPACING_NM = 2.8
STIFFNESS_SPACING_NM = 2.0
DSDNA_STRETCH_PN = 1100.0
DSDNA_BEND_PN_NM2 = 230.0
TWIST_DOMAIN_MAX_INDEL = 3

# Legacy lattice-wide Twist calibration and current Bend calibration.  The
# affine Twist entries are retained only for old-file provenance; new Add and
# Remove Twist predictions use the W x L table below.  Bend still maps the
# elastic prediction onto the corresponding SNUPI result:
#
#     calibrated = intercept + slope * uncalibrated
#
# Twist was fitted from 13 converged systems per lattice; Bend from 4 systems
# per lattice.  The raw predictions remain available in every result so the
# calibration is auditable and can be replaced without losing the underlying
# physical-model output.
SNUPI_CALIBRATION_VERSION = \
    '2026-08-22-WL-topology-S8-R4x4C-nominal-angle'
SNUPI_CALIBRATION = {
    'twist': {
        'square': {
            'intercept': -0.08908547710384841,
            'slope': 0.7801627015505553,
            'rmse': 0.0642047928338592,
            'count': 13,
        },
        'honeycomb': {
            'intercept': -0.1670919226913426,
            'slope': 0.7600602565841517,
            'rmse': 0.1250751350261012,
            'count': 13,
        },
    },
    'bend': {
        'square': {
            'intercept': 0.0,
            'slope': 0.8679254814050458,
            'rmse': 1.5425177622409517,
            'count': 4,
        },
        'honeycomb': {
            'intercept': 0.0,
            'slope': 1.0052559269975638,
            'rmse': 0.811662494514633,
            'count': 4,
        },
    },
}

# Coarse-grained Twist calibration for solid, regular, single W x L bundle
# cross-sections.  Each tuple is
#
#     (native elastic x0, non-negative native anchor T0,
#      deletion slope k-, insertion slope k+, split branches)
#
# and is evaluated as T0 + k-*(x-x0) below native, or T0 + k+*(x-x0)
# above native.  The parameters are interpolated in W and L; the final Twist
# itself is never interpolated.  These models intentionally do not apply to
# pores, incomplete sections, separated modules or ambiguous cross-sections.
TWIST_CROSS_SECTION_CALIBRATION = {
    'honeycomb': {
        (4, 1): (.900000, .440703, 1.564749, 1.685614, True),
        (8, 1): (.399857, .375734, 1.456099, 1.070584, True),
        (12, 1): (.314508, .341243, 1.461079, .943853, True),
        (4, 2): (.624229, 0.0, 1.101283, 1.019071, True),
        (8, 2): (.341437, 0.0, .863044, 1.081851, True),
        (12, 2): (.260927, 0.0, .838808, .939565, True),
        (4, 3): (.403168, .153935, .940880, .940880, False),
        (8, 3): (.299611, .088676, .668059, .668059, False),
        (4, 4): (.353844, 0.0, .821023, .821023, False),
        (8, 4): (.264216, 0.0, .500344, .500344, False),
        (6, 6): (.237096, .005845, .551209, .551209, False),
        (10, 6): (.112159, 0.0, .450682, .450682, False),
        (8, 8): (.094809, 0.0, .504727, .504727, False),
    },
    'square': {
        (4, 1): (.825093, .982838, 1.515103, 1.295900, True),
        (8, 1): (.550824, .795328, 1.385637, 1.278554, True),
        (12, 1): (.252395, .753354, 2.758262, 3.441981, True),
        (4, 2): (.711564, .654682, 1.135997, 1.076154, True),
        (8, 2): (.368969, .508484, 1.450367, 1.204797, True),
        (12, 2): (.133676, .472648, 4.640838, 2.525380, True),
        (4, 3): (.627925, .351324, .981723, .981723, False),
        (8, 3): (.251544, .279715, 1.497573, 1.497573, False),
        (4, 4): (.556618, .297149, .844081, .844081, False),
        (8, 4): (.158020, .180127, 2.723007, 2.723007, False),
        (6, 6): (.155166, .086511, 2.365168, 2.365168, False),
        (10, 6): (.077606, .081905, 2.352061, 2.352061, False),
        (8, 8): (.077359, .047478, 1.975104, 1.975104, False),
    },
}

# Source-data anchors from Kube et al., Nat Commun 2020, Fig. 2e.  The
# publication's sign is inverted here so positive means the right-handed
# rotation used by this editor.  The published Square source-data table lists
# the 6x6 and 8x8 magnitudes as 0.078125 and 0.15625 deg/base respectively.
# That isolated reversal is inconsistent with the expected monotonic decrease
# in global twist as the cross-section becomes stiffer, and it also conflicts
# with our SNUPI calibration trend.  The application therefore deliberately
# uses the physically monotonic order (6x6=0.15625, 8x8=0.078125).  Keep this
# documented model correction distinct from a transcription of the raw table.
_SQUARE_J = (31.0, 527.0, 2695.0, 8545.0)
_SQUARE_TWIST_1067 = (.9375, .5625, .15625, .078125)
_SQUARE_TWIST_1033 = (.265625, 0.0, .046875, 0.0)
_HONEYCOMB_J = (85.0, 277.0, 4027.0, 6003.0, 8557.0, 57000.0)
_HONEYCOMB_TWIST_1050 = (.9, .428571, .238095, .138889, .120482, 0.0)


class TwistBendError(ValueError):
    pass


def _lattice_name(value):
    """Normalize a lattice name or helical pitch to a calibration key."""
    if isinstance(value, str):
        return ('honeycomb' if value.lower().startswith(('h', 'hon'))
                else 'square')
    return 'honeycomb' if float(value) < 10.58 else 'square'


def _snupi_calibrate(kind, lattice, value, preserve_zero=False):
    """Apply the final affine calibration while retaining a physical zero."""
    value = float(value)
    model = SNUPI_CALIBRATION[kind][_lattice_name(lattice)]
    if preserve_zero and abs(value) <= 1e-12:
        return 0.0
    return model['intercept'] + model['slope'] * value


def _cluster_axis(values, tolerance):
    """Cluster nearly equal physical coordinates without changing order."""
    groups = []
    for value in sorted(float(item) for item in values):
        if not groups or abs(value-sum(groups[-1])/len(groups[-1])) > tolerance:
            groups.append([value])
        else:
            groups[-1].append(value)
    return [sum(group)/float(len(group)) for group in groups]


def _nearest_axis(value, centers):
    return min(range(len(centers)), key=lambda index: abs(value-centers[index]))


def classify_regular_cross_section(units, lattice):
    """Recognize one complete solid W x L cross-section.

    The classifier uses the real 2.8-nm mechanical coordinates, not helix
    numbering.  It accepts a complete rectangular Square section or complete
    staggered Honeycomb rows, and rejects holes, gaps, separated modules and
    ambiguous spacing.  W is the horizontal row width and L the row count in
    the left cross-section view.
    """
    lattice = _lattice_name(lattice)
    # Helices generated for calibration are often axially staggered while
    # retaining the same physical dsDNA length.  That is a valid mechanical
    # cross-section.  Reject only genuinely unequal/truncated axial coverage,
    # not equal-length helices whose coordinate endpoints differ.
    native_counts = [int(unit.get(
        'native_bases', unit['end']-unit['start']+1)) for unit in units]
    coverage_compatible = (not native_counts or
                           len(set(native_counts)) == 1)
    if (any(not unit.get('coverage_complete', True) for unit in units) and
            not coverage_compatible):
        return None
    lattice_points = [unit.get('lattice_coord') for unit in units]
    if lattice_points and all(point is not None for point in lattice_points):
        lattice_points = [(int(point[0]), int(point[1]))
                          for point in lattice_points]
        rows = sorted(set(point[0] for point in lattice_points))
        columns = sorted(set(point[1] for point in lattice_points))
        complete = set((row, column) for row in rows for column in columns)
        consecutive_rows = (not rows or rows == list(range(rows[0], rows[-1]+1)))
        consecutive_columns = (not columns or
            columns == list(range(columns[0], columns[-1]+1)))
        if (len(lattice_points) == len(set(lattice_points)) and
                set(lattice_points) == complete and consecutive_rows and
                consecutive_columns):
            return {'width': len(columns), 'layers': len(rows),
                    'spacing_nm': CENTER_SPACING_NM, 'lattice': lattice,
                    'solid_regular': True}
        return None
    points = [(float(unit['coord'][0]), float(unit['coord'][1]))
              for unit in units]
    if len(points) < 2 or len(set((round(x, 6), round(y, 6))
                                 for x, y in points)) != len(points):
        return None
    distances = []
    for index, first in enumerate(points):
        for second in points[index+1:]:
            distance = math.hypot(first[0]-second[0], first[1]-second[1])
            if distance > 1e-6:
                distances.append(distance)
    if not distances:
        return None
    spacing = min(distances)
    tolerance = max(.08, spacing*.18)
    ys = _cluster_axis([point[1] for point in points], tolerance)
    rows = [[] for unused in ys]
    for x, y in points:
        rows[_nearest_axis(y, ys)].append(x)
    widths = [len(row) for row in rows]
    if not widths or min(widths) < 1 or len(set(widths)) != 1:
        return None
    width, layers = widths[0], len(rows)
    if width*layers != len(points):
        return None
    row_starts = []
    horizontal = []
    for row in rows:
        ordered = sorted(row)
        row_starts.append(ordered[0])
        gaps = [ordered[index+1]-ordered[index]
                for index in range(len(ordered)-1)]
        horizontal.extend(gaps)
        if gaps and (max(gaps)-min(gaps) > tolerance or
                     abs(sum(gaps)/len(gaps)-spacing) > tolerance):
            return None
    if layers > 1:
        vertical = [ys[index+1]-ys[index] for index in range(layers-1)]
        if max(vertical)-min(vertical) > tolerance:
            return None
    if width > 1 and horizontal:
        row_step = sum(horizontal)/float(len(horizontal))
    else:
        row_step = spacing
    if lattice == 'square':
        # Rows must align and both lattice spacings must be the same.
        if max(row_starts)-min(row_starts) > tolerance:
            return None
        if layers > 1 and abs(sum(vertical)/len(vertical)-row_step) > tolerance:
            return None
    else:
        # Honeycomb rows alternate by approximately half a helix spacing and
        # are separated by sqrt(3)/2 of that spacing.
        if layers > 1:
            vertical_mean = sum(vertical)/len(vertical)
            if abs(vertical_mean-row_step*math.sqrt(3.0)/2.0) > tolerance:
                return None
            offsets = [((value-row_starts[0])/row_step) % 1.0
                       for value in row_starts]
            normalized = [min(value, 1.0-value) for value in offsets]
            if any(min(abs(value), abs(value-.5)) > .12
                   for value in normalized):
                return None
    return {'width': int(width), 'layers': int(layers),
            'spacing_nm': float(row_step), 'lattice': lattice,
            'solid_regular': True}


def _linear_profile(first, second, fraction):
    fraction = max(0.0, min(1.0, float(fraction)))
    values = [first[index] + fraction*(second[index]-first[index])
              for index in range(4)]
    return tuple(values) + (bool(first[4] or second[4]),)


def _profile_at_width(nodes, layer, width):
    available = sorted((node_width, profile)
                       for (node_width, node_layer), profile in nodes.items()
                       if node_layer == layer)
    if not available:
        return None
    for node_width, profile in available:
        if node_width == width:
            return profile
    lower = [item for item in available if item[0] < width]
    upper = [item for item in available if item[0] > width]
    if lower and upper:
        left, right = lower[-1], upper[0]
        return _linear_profile(left[1], right[1],
                               (width-left[0])/float(right[0]-left[0]))
    return (lower[-1] if lower else upper[0])[1]


def interpolate_twist_cross_section(lattice, width, layers):
    """Interpolate x0, T0, k- and k+; never interpolate final Twist."""
    lattice = _lattice_name(lattice)
    width, layers = int(width), int(layers)
    nodes = TWIST_CROSS_SECTION_CALIBRATION[lattice]
    if (width, layers) in nodes:
        return {'profile': nodes[(width, layers)], 'exact': True,
                'extrapolated': False}
    layer_values = sorted(set(layer for unused, layer in nodes))
    if layers in layer_values:
        widths = sorted(node_width for node_width, node_layer in nodes
                        if node_layer == layers)
        return {'profile': _profile_at_width(nodes, layers, width),
                'exact': False,
                'extrapolated': width < min(widths) or width > max(widths)}
    lower = [layer for layer in layer_values if layer < layers]
    upper = [layer for layer in layer_values if layer > layers]
    if lower and upper:
        low_layer, high_layer = lower[-1], upper[0]
        low = _profile_at_width(nodes, low_layer, width)
        high = _profile_at_width(nodes, high_layer, width)
        profile = _linear_profile(
            low, high, (layers-low_layer)/float(high_layer-low_layer))
        return {'profile': profile, 'exact': False, 'extrapolated': False}
    # Beyond the calibrated boundary, continue the nearest W/L slopes.  The
    # native anchor tends to zero for thick sections; Honeycomb reaches zero
    # immediately and Square approaches it monotonically over two extra rows.
    boundary_layer = (lower[-1] if lower else upper[0])
    profile = list(_profile_at_width(nodes, boundary_layer, width))
    if layers > max(layer_values):
        if lattice == 'honeycomb':
            profile[1] = 0.0
        else:
            fade = max(0.0, 1.0-(layers-max(layer_values))/2.0)
            profile[1] = max(0.0, profile[1]*fade)
    return {'profile': tuple(profile), 'exact': False,
            'extrapolated': True}


def _calibrate_cross_section_twist(lattice, raw_value, section):
    interpolation = interpolate_twist_cross_section(
        lattice, section['width'], section['layers'])
    x0, t0, k_minus, k_plus, split = interpolation['profile']
    delta = float(raw_value)-x0
    value = t0 + (k_minus if delta < 0 else k_plus)*delta
    return value, {
        'calibration_x0': x0, 'calibration_T0': t0,
        'calibration_k_deletion': k_minus,
        'calibration_k_insertion': k_plus,
        'calibration_split_branches': split,
        'calibration_exact_node': interpolation['exact'],
        'calibration_extrapolated': interpolation['extrapolated']}


def calibrate_saved_twist_prediction(prediction, lattice=None):
    """Upgrade one prediction stored by a pre-calibration JSON task."""
    result = dict(prediction or {})
    if not result or result.get('calibration_version') == SNUPI_CALIBRATION_VERSION:
        return result
    lattice = _lattice_name(lattice or result.get('lattice', 'square'))
    connectivity = float(result.get('connectivity_fraction', 1.0))
    old_value = float(result.get('uncalibrated_twist_per_base_deg',
                                 result.get('twist_per_base_deg', 0.0)))
    connected_value = (old_value / connectivity
                       if connectivity > 1e-12 else 0.0)
    width = result.get('cross_section_width')
    layers = result.get('cross_section_layers')
    applicable = bool(result.get('cross_section_solid_regular') and
                      width and layers and connectivity >= .999)
    metadata = {}
    if applicable:
        calibrated, metadata = _calibrate_cross_section_twist(
            lattice, connected_value,
            {'width': int(width), 'layers': int(layers)})
        calibrated *= connectivity
    elif (result.get('irregular_calibrated') and
          result.get('irregular_profile') and connectivity >= .999):
        profile = tuple(float(value) for value in
                        result['irregular_profile'])
        neighbors = tuple(result.get('irregular_neighbors') or ())
        exact_s8 = bool(
            result.get('irregular_calibration_exact_node') and neighbors and
            neighbors[0][0] == 'S8-R4x4C')
        if exact_s8:
            profile = S8_R4X4C_EXTENDED_PROFILE
        x0, t0, k_minus, k_plus = profile
        calibrated = (t0 + (k_minus if connected_value < x0 else k_plus) *
                      (connected_value-x0)) * connectivity
        applicable = True
        metadata.update({
            'irregular_profile': profile,
            'irregular_calibrated': True,
        })
        if exact_s8:
            metadata.update({
                'irregular_calibration_version':
                    S8_R4X4C_CALIBRATION_VERSION,
                'irregular_calibration_valid_indel_range':
                    S8_R4X4C_VALIDATED_INDEL_RANGE,
                'irregular_calibration_failed_indel_points':
                    S8_R4X4C_FAILED_INDEL_POINTS,
                'irregular_calibration_points': S8_R4X4C_TWIST_POINTS,
                'irregular_calibration_branch_rmse':
                    dict(S8_R4X4C_BRANCH_RMSE),
                'irregular_failed_points_excluded_from_fit': True,
            })
    else:
        # Old task records do not contain enough geometry to prove that they
        # belong to the new solid W x L calibration domain.
        calibrated = old_value
    # Twist rates are calibrated per nominal design base.  Existing indels
    # remain part of the effective-pitch calculation, but must not enlarge the
    # integration span a second time when converting deg/base to total angle.
    # Pre-2026-08-22 records stored the actual length in ``length_bp``; their
    # existing mean indel load is sufficient to recover the nominal span.
    stored_length = float(result.get('length_bp', 0.0))
    actual_length = float(result.get('actual_length_bp', stored_length))
    if 'nominal_length_bp' in result:
        length = float(result['nominal_length_bp'])
    elif 'mean_indel_per_helix' in result:
        length = stored_length-float(result['mean_indel_per_helix'])
    else:
        length = stored_length
    result.update({
        'lattice': lattice,
        'uncalibrated_twist_per_base_deg': old_value,
        'twist_per_base_deg': calibrated,
        'total_twist_deg': calibrated * length,
        'length_bp': length,
        'nominal_length_bp': length,
        'actual_length_bp': actual_length,
        'twist_integration_length_basis': 'nominal',
        'handedness': ('右手' if calibrated > 1e-9 else
                       '左手' if calibrated < -1e-9 else '近似无扭转'),
        'calibration_version': SNUPI_CALIBRATION_VERSION,
        'calibrated': applicable,
        'cross_section_calibration_applicable': applicable,
    })
    result.update(metadata)
    return result


def calibrate_saved_bend_prediction(prediction, lattice='square'):
    """Upgrade one elastic Bend prediction stored by an older JSON task."""
    result = dict(prediction or {})
    if not result or result.get('calibration_version') == SNUPI_CALIBRATION_VERSION:
        return result
    lattice = _lattice_name(lattice or result.get('lattice', 'square'))
    raw = float(result.get('uncalibrated_angle_degrees',
                           result.get('angle_degrees', 0.0)))
    angle = _snupi_calibrate('bend', lattice, raw, preserve_zero=True)
    length_nm = float(result.get('length_nm', 0.0))
    curvature = (math.radians(angle) / length_nm
                 if length_nm > 1e-12 else 0.0)
    model = SNUPI_CALIBRATION['bend'][lattice]
    result.update({
        'lattice': lattice,
        'uncalibrated_angle_degrees': raw,
        'angle_degrees': angle,
        'curvature_per_nm': curvature,
        'radius_nm': (1.0 / curvature if curvature > 1e-12 else None),
        'calibration_version': SNUPI_CALIBRATION_VERSION,
        'calibration_intercept': model['intercept'],
        'calibration_slope': model['slope'],
        'calibrated': True,
    })
    return result


def calibrate_saved_plan(plan, lattice='square'):
    """Upgrade saved predictions and preview angles without changing edits."""
    result = dict(plan or {})
    lattice = _lattice_name(lattice)
    for key in ('twist_before_prediction', 'twist_prediction'):
        if result.get(key):
            result[key] = calibrate_saved_twist_prediction(
                result[key], lattice)
    if result.get('elastic_prediction'):
        result['elastic_prediction'] = calibrate_saved_bend_prediction(
            result['elastic_prediction'], lattice)
        result['radius_nm'] = result['elastic_prediction']['radius_nm']
    preview = dict(result.get('preview_transform') or {})
    if preview:
        if result.get('kind') == 'bend' and result.get('elastic_prediction'):
            preview['angle'] = result['elastic_prediction']['angle_degrees']
        elif (result.get('twist_prediction') and
              result.get('twist_before_prediction')):
            preview['angle'] = (
                result['twist_prediction']['total_twist_deg'] -
                result['twist_before_prediction']['total_twist_deg'])
        result['preview_transform'] = preview
    return result


def inclusive_length(start, end, existing=None):
    """Return actual nucleotides in an inclusive lattice-index interval."""
    start, end = sorted((int(start), int(end)))
    value = end - start + 1
    for idx, length in (existing or {}).items():
        if start <= int(idx) <= end:
            value += int(length)
    return value


def validate_regions(regions):
    """Reject overlapping tasks that would prescribe incompatible edits."""
    normalized = []
    for region_index, region in enumerate(regions):
        start, end = sorted((int(region['start']), int(region['end'])))
        if start == end:
            raise TwistBendError('区域 %d 至少需要包含 2 个碱基位置。' %
                                 (region_index + 1))
        helices = tuple(sorted(set(int(v) for v in region.get('helices', ()))))
        if not helices:
            raise TwistBendError('区域 %d 尚未选择 helix。' %
                                 (region_index + 1))
        current = dict(region)
        current.update(start=start, end=end, helices=helices)
        # Tasks are an ordered virtual-design pipeline.  Overlap is therefore
        # intentional: a later Remove/Add Twist or Bend starts from the
        # predicted result of every earlier task.  Individual edit sites are
        # still protected by the planner's forbidden-base checks.
        normalized.append(current)
    return normalized


def _even_positions(start, end, count, allowed):
    """Choose ``count`` allowed positions as evenly as possible."""
    allowed = sorted(set(int(v) for v in allowed if start <= int(v) <= end))
    if count <= 0:
        return []
    if len(allowed) < count:
        raise TwistBendError(
            'base %d–%d 内只有 %d 个安全 indel 位点，但需要 %d 个。' %
            (start, end, len(allowed), count))
    chosen = []
    remaining = list(allowed)
    for rank in range(count):
        target = start + (rank + 1) * (end - start) / float(count + 1)
        candidate = min(remaining, key=lambda value: (abs(value - target), value))
        chosen.append(candidate)
        remaining.remove(candidate)
    return sorted(chosen)


def _allocate_integer(total, weights, capacities=None):
    """Balanced largest-remainder allocation with optional capacities."""
    total = int(max(0, total))
    if not weights:
        return []
    weights = [max(0.0, float(value)) for value in weights]
    denominator = sum(weights) or float(len(weights))
    raw = [total * (value if sum(weights) else 1.0) / denominator
           for value in weights]
    values = [int(math.floor(value)) for value in raw]
    capacities = ([int(value) for value in capacities]
                  if capacities is not None else [total] * len(weights))
    values = [min(value, capacities[index])
              for index, value in enumerate(values)]
    remaining = total - sum(values)
    order = sorted(range(len(values)),
                   key=lambda index: (raw[index] - values[index],
                                      weights[index], -index), reverse=True)
    while remaining > 0:
        progressed = False
        for index in order:
            if values[index] < capacities[index]:
                values[index] += 1
                remaining -= 1
                progressed = True
                if remaining == 0:
                    break
        if not progressed:
            raise TwistBendError('安全位点不足，无法分配所需 indel。')
    return values


def _allocate_evenly_across_units(total, capacities):
    """Spread integer edits over the selected helices, including fractions.

    For example, 0.1 edit/helix over 16 helices rounds to two edits and
    places them near the quarter and three-quarter positions, rather than on
    the first two helices.  Additional layers are distributed the same way.
    """
    total = int(max(0, total))
    capacities = [max(0, int(value)) for value in capacities]
    values = [0] * len(capacities)
    if total > sum(capacities):
        raise TwistBendError('安全位点不足，无法分配所需 indel。')
    remaining = total
    while remaining:
        available = [index for index, capacity in enumerate(capacities)
                     if values[index] < capacity]
        take = min(remaining, len(available))
        chosen = []
        for rank in range(take):
            position = int((rank + .5) * len(available) / float(take))
            position = min(len(available)-1, position)
            index = available[position]
            if index not in chosen:
                chosen.append(index)
        # Rounding above is normally unique; retain a deterministic fallback.
        if len(chosen) < take:
            chosen.extend(index for index in available if index not in chosen)
        for index in chosen[:take]:
            values[index] += 1
            remaining -= 1
    return values


def _twist_domain_size(lattice_pitch):
    """Return the native domain size used by Add/Remove Twist placement."""
    return 7 if _lattice_name(lattice_pitch) == 'honeycomb' else 8


def _edit_adjustments_by_helix(edits):
    """Group already-selected edit deltas for final-domain load accounting."""
    grouped = {}
    for edit in edits or ():
        grouped.setdefault(int(edit['helix']), []).append(
            (int(edit['idx']), int(edit['length'])))
    return grouped


def _twist_domain_context(unit, edit_length, domain_size, adjustments=()):
    """Build allowed sites and signed final-load capacities per native domain.

    Existing insertions/deletions and any selected removal edits are included
    before new sites are allocated.  A deletion is additionally forbidden on
    a staple oligo that is already at the 21-nt hard minimum.
    """
    start, end = int(unit['start']), int(unit['end'])
    domain_ids = list(range(start // domain_size, end // domain_size + 1))
    loads = dict((domain, 0) for domain in domain_ids)
    for idx, length in unit.get('existing', {}).items():
        domain = int(idx) // domain_size
        if domain in loads:
            loads[domain] += int(length)
    for idx, length in adjustments:
        domain = int(idx) // domain_size
        if domain in loads:
            loads[domain] += int(length)

    protected = (set(unit.get('deletion_protected', ()))
                 if int(edit_length) < 0 else set())
    candidates = dict((domain, []) for domain in domain_ids)
    for idx in unit['allowed']:
        idx = int(idx)
        if idx in protected:
            continue
        domain = idx // domain_size
        if domain in candidates:
            candidates[domain].append(idx)

    capacities = {}
    for domain in domain_ids:
        load = loads[domain]
        signed_room = (TWIST_DOMAIN_MAX_INDEL-load if edit_length > 0
                       else load+TWIST_DOMAIN_MAX_INDEL)
        capacities[domain] = min(len(candidates[domain]),
                                 max(0, int(signed_room)))
    return domain_ids, loads, candidates, capacities


def _twist_partition_assignments(total, start, end, domain_size,
                                 domain_ids, capacities):
    """Create one axial target per equal partition, then enforce domain caps.

    Axial uniformity is defined by the requested edit count ``total``: the
    selected region is split into ``total`` equal bins and each bin contributes
    one target at its centre.  Absolute 7/8-bp domains are capacity/safety
    buckets only; they do not define the primary spatial distribution.
    """
    total = int(max(0, total))
    if total == 0:
        return []
    if total > sum(capacities.values()):
        raise TwistBendError(
            '7/8-bp domain 的安全容量不足：需要 %d 个 indel，最多只能放置 %d 个。' %
            (total, sum(capacities.values())))

    width = float(int(end)-int(start)+1)
    assignments = []
    for rank in range(total):
        # Treat integer base coordinates as base centres, so the continuous
        # selected interval spans start-0.5 through end+0.5.
        target = int(start)-.5 + (rank+.5)*width/float(total)
        target_base = int(math.floor(target))
        domain = target_base // int(domain_size)
        domain = min(domain_ids, key=lambda value: (
            abs(value-domain), value))
        partition_start = int(start) + int(math.floor(
            rank*width/float(total)))
        partition_end = int(start) + int(math.floor(
            (rank+1)*width/float(total))) - 1
        assignments.append({'target': target, 'domain': domain,
                            'partition': rank,
                            'partition_start': partition_start,
                            'partition_end': partition_end})
    return assignments


def _choose_twist_partition_sites(candidates, assignments, capacities,
                                  allow_repeated_sites=False):
    """Resolve equal-partition targets under absolute-domain constraints.

    For twist, insertion and deletion contributions are both domain based.
    No artificial cross-helix phase offset is applied.  Selection first uses
    the intersection of the target's equal axial partition and absolute
    7/8-bp domain, then another domain inside that partition.  Only when the
    entire partition lacks a safe site may a nearest-site repair leave it.
    """
    assignments = sorted(assignments, key=lambda item: item['target'])
    if not assignments:
        return []
    remaining = dict((int(domain), sorted(set(int(value) for value in values)))
                     for domain, values in candidates.items())
    room = dict((int(domain), int(value))
                for domain, value in capacities.items())
    chosen = []
    unresolved = []

    def options_for(assignment, inside_partition):
        values = []
        for domain, sites in remaining.items():
            if room.get(domain, 0) <= 0:
                continue
            for site in sites:
                inside = (assignment['partition_start'] <= site <=
                          assignment['partition_end'])
                if inside != inside_partition:
                    continue
                values.append((domain, site))
        return values

    def reserve(assignment, options):
        preferred = int(assignment['domain'])
        domain, site = min(options, key=lambda item: (
            item[0] != preferred,
            abs(float(item[1])-float(assignment['target'])),
            abs(item[0]-preferred), item[1]))
        if not allow_repeated_sites:
            remaining[domain].remove(site)
        room[domain] -= 1
        chosen.append({'idx': site, 'domain': domain,
                       'assignment': assignment})

    # Resolve every bin that has an in-bin safe site before any fallback can
    # consume a site belonging to a later, disjoint equal partition.
    for assignment in assignments:
        options = options_for(assignment, True)
        if options:
            reserve(assignment, options)
        else:
            unresolved.append(assignment)
    for assignment in unresolved:
        options = options_for(assignment, False)
        if not options:
            raise TwistBendError('安全 domain 容量不足，无法修复 indel 分配。')
        reserve(assignment, options)
    return sorted(chosen, key=lambda item: item['idx'])


def equal_partition_indel_sites(total, start, end, domain_size,
                                candidates, capacities,
                                allow_repeated_sites=False):
    """Return safe indel sites using the shared equal-partition rule.

    ``candidates`` and ``capacities`` are keyed by absolute native-domain
    number.  The requested interval is first divided into ``total`` disjoint
    axial partitions.  One target is created at the centre of every
    partition, then resolved to the nearest safe site inside both that
    partition and its preferred absolute 7/8-bp domain.  Existing safety and
    per-domain capacity constraints remain authoritative; leaving a partition
    is only a last-resort repair when it contains no usable site at all.

    This public helper is shared by Add/Remove Twist, Add Bending and Curved
    Design so those tools cannot silently drift back to different axial
    distribution rules.
    """
    domain_ids = sorted(int(value) for value in candidates)
    if not domain_ids:
        if int(total) <= 0:
            return []
        raise TwistBendError('选区内没有可用的 7/8-bp domain。')
    assignments = _twist_partition_assignments(
        total, start, end, domain_size, domain_ids, capacities)
    return _choose_twist_partition_sites(
        candidates, assignments, capacities,
        allow_repeated_sites=allow_repeated_sites)


def _twist_unit_capacity(unit, edit_length, domain_size, adjustments=()):
    return sum(_twist_domain_context(
        unit, edit_length, domain_size, adjustments)[3].values())


def _domain_aware_twist_edits(units, counts, edit_length, domain_size,
                              prior_edits=()):
    """Place Add/Remove Twist edits using equal partitions and domain quotas."""
    adjustments = _edit_adjustments_by_helix(prior_edits)
    edits = []
    quota_report = {}
    final_load_report = {}
    maximum_final_load = 0
    for unit_index, (unit, count) in enumerate(zip(units, counts)):
        number = int(unit['helix'])
        domain_ids, loads, candidates, capacities = _twist_domain_context(
            unit, edit_length, domain_size, adjustments.get(number, ()))
        quotas = dict((domain, 0) for domain in domain_ids)
        selected = equal_partition_indel_sites(
            count, unit['start'], unit['end'], domain_size,
            candidates, capacities)
        for item in selected:
            idx, domain = int(item['idx']), int(item['domain'])
            assignment = item['assignment']
            quotas[domain] += 1
            edits.append({'helix': number, 'idx': idx,
                          'length': int(edit_length),
                          'partition': int(assignment['partition']),
                          'partition_count': int(count),
                          'ideal_base': float(assignment['target']),
                          'partition_start': int(
                              assignment['partition_start']),
                          'partition_end': int(
                              assignment['partition_end']),
                          'target_domain': int(assignment['domain']),
                          'domain': domain,
                          'domain_size': int(domain_size)})
            loads[domain] += int(edit_length)
        for domain in domain_ids:
            if abs(loads[domain]) > TWIST_DOMAIN_MAX_INDEL:
                raise TwistBendError(
                    'helix %d 的 domain %d 最终 indel 负载为 %+d，超过 ±%d。' %
                    (number, domain, loads[domain],
                     TWIST_DOMAIN_MAX_INDEL))
            maximum_final_load = max(maximum_final_load, abs(loads[domain]))
        quota_report[number] = dict(quotas)
        final_load_report[number] = dict(loads)
    return edits, {
        'domain_size_bp': int(domain_size),
        'maximum_indel_per_domain_allowed': TWIST_DOMAIN_MAX_INDEL,
        'maximum_final_indel_per_domain_observed': maximum_final_load,
        'domain_indel_quotas': quota_report,
        'final_domain_indel_loads': final_load_report,
        'indel_distribution_method':
            ('equal-partition-and-native-domain-intersection-first; '
             'nearest-safe-site repair; no forced stagger')}


def _domain_aware_signed_edits(units, signed_counts, domain_size):
    """Place mixed-sign bending edits with the same rule used by Twist."""
    edits = []
    quota_report = {}
    final_load_report = {}
    maximum_final_load = 0
    for unit, signed_count in zip(units, signed_counts):
        signed_count = int(signed_count)
        number = int(unit['helix'])
        if signed_count == 0:
            domain_ids, loads, unused_candidates, unused_capacities = \
                _twist_domain_context(unit, 1, domain_size)
            quota_report[number] = dict((domain, 0)
                                        for domain in domain_ids)
            final_load_report[number] = dict(loads)
            maximum_final_load = max(
                [maximum_final_load] + [abs(value)
                                        for value in loads.values()])
            continue
        edit_length = 1 if signed_count > 0 else -1
        placed, metadata = _domain_aware_twist_edits(
            [unit], [abs(signed_count)], edit_length, domain_size)
        edits.extend(placed)
        quota_report[number] = metadata['domain_indel_quotas'][number]
        final_load_report[number] = metadata[
            'final_domain_indel_loads'][number]
        maximum_final_load = max(
            maximum_final_load,
            int(metadata['maximum_final_indel_per_domain_observed']))
    return edits, {
        'domain_size_bp': int(domain_size),
        'maximum_indel_per_domain_allowed': TWIST_DOMAIN_MAX_INDEL,
        'maximum_final_indel_per_domain_observed': maximum_final_load,
        'domain_indel_quotas': quota_report,
        'final_domain_indel_loads': final_load_report,
        'indel_distribution_method':
            ('equal-partition-and-native-domain-intersection-first; '
             'nearest-safe-site repair; no forced stagger')}


def _region_units(region, helix_data):
    units = []
    for number in region['helices']:
        data = helix_data.get(number)
        if data is None:
            raise TwistBendError('找不到 helix %d。' % number)
        existing = data.get('insertions', {})
        scaffold_intervals = list(data.get('scaffold_intervals', ()))
        staple_intervals = list(data.get('staple_intervals', ()))
        paired = None
        if scaffold_intervals or staple_intervals:
            scaffold = set()
            staple = set()
            for low, high in scaffold_intervals:
                scaffold.update(range(max(region['start'], int(low)),
                                      min(region['end'], int(high))+1))
            for low, high in staple_intervals:
                staple.update(range(max(region['start'], int(low)),
                                   min(region['end'], int(high))+1))
            paired = scaffold & staple
            # Scaffold-only or staple-only helices do not contribute to the
            # double-stranded mechanical cross-section and receive no indel.
            if not paired:
                continue
        base_indices = (sorted(paired) if paired is not None else
                        list(range(region['start'], region['end']+1)))
        base_set = set(base_indices)
        actual = sum(max(0, 1+int(existing.get(idx, 0)))
                     for idx in base_indices)
        existing_delta = sum(int(existing.get(idx, 0))
                             for idx in base_indices)
        allowed = [idx for idx in base_indices
                   if region['start'] < idx < region['end']
                   if idx not in data.get('forbidden', set()) and
                   idx not in existing]
        if paired is not None:
            coverage_complete = all(idx in paired for idx in
                                    range(region['start'], region['end']+1))
        else:
            # Synthetic callers and legacy saved tasks may not provide strand
            # coverage. Geometry tests remain valid, while the live editor
            # always supplies these intervals.
            coverage_complete = True
        units.append({'helix': number, 'start': region['start'],
                      'end': region['end'], 'actual': actual,
                      'native_bases': len(base_indices),
                      'existing_delta': existing_delta,
                      'existing': dict((int(idx), int(length))
                                       for idx, length in existing.items()
                                       if idx in base_set),
                      'allowed': allowed, 'coord': data['coord'],
                      'deletion_protected': set(
                          int(idx) for idx in
                          data.get('deletion_protected', set())
                          if idx in base_set),
                      'lattice_coord': data.get('lattice_coord'),
                      'coverage_complete': coverage_complete})
    if not units:
        raise TwistBendError('选区内没有可用于 Twist/Bending 的双链 DNA。')
    return units


def _square_irregular_descriptor(region, units, helix_data):
    """Describe one Square outer frame with missing double-stranded helices."""
    selected_points = []
    for number in region['helices']:
        item = helix_data.get(number, {})
        point = item.get('lattice_coord')
        if point is not None:
            selected_points.append((int(point[0]), int(point[1])))
    native_counts = [int(unit.get(
        'native_bases', unit['end']-unit['start']+1)) for unit in units]
    # Equal dsDNA lengths with staggered endpoints are the normal input form
    # of the missing-helix calibration set.  A genuinely shortened helix still
    # invalidates the node, preserving the partial-coverage safety rule.
    if native_counts and len(set(native_counts)) != 1:
        return None
    present_raw = [unit.get('lattice_coord') for unit in units
                   if unit.get('lattice_coord') is not None]
    if not selected_points or not present_raw:
        return None
    all_rows = [point[0] for point in selected_points]
    all_cols = [point[1] for point in selected_points]
    row0, col0 = min(all_rows), min(all_cols)
    row_span = max(all_rows)-row0+1
    col_span = max(all_cols)-col0+1
    if row_span != col_span or row_span not in (4, 6, 8):
        return None
    size = row_span
    present = set((int(point[0])-row0, int(point[1])-col0)
                  for point in present_raw)
    expected = set((row, col) for row in range(size)
                   for col in range(size))
    if not present <= expected:
        return None
    missing = sorted(expected-present)
    if not missing:
        return None
    center = (size-1)/2.0
    mx = sum(point[0] for point in missing)/float(len(missing))
    my = sum(point[1] for point in missing)/float(len(missing))
    eccentricity = (math.hypot(mx-center, my-center) /
                    max(math.sqrt(2.0)*center, 1e-12))
    degrees, depths = [], []
    boundary = 0
    for row, col in present:
        degree = sum((row+dr, col+dc) in present
                     for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)))
        degrees.append(degree)
        if degree < 4:
            boundary += 1
        depths.append(min(abs(row-r)+abs(col-c)
                          for r in range(-1, size+1)
                          for c in range(-1, size+1)
                          if (r, c) not in present))
    actual_j = _cross_section_polar_moment(units)
    solid_units = [{'coord': (row*CENTER_SPACING_NM,
                              col*CENTER_SPACING_NM)}
                   for row in range(size) for col in range(size)]
    solid_j = _cross_section_polar_moment(solid_units)
    descriptor = (
        len(missing)/float(size*size), eccentricity,
        actual_j/max(solid_j, 1e-12), boundary/float(len(present)),
        sum(degrees)/float(len(degrees)),
        sum(depths)/float(len(depths)))
    rows = [point[0] for point in missing]
    cols = [point[1] for point in missing]
    gaps = (min(rows), min(cols), size-1-max(rows), size-1-max(cols))
    one_cell_wall = any(value == 1 for value in gaps if value > 0)
    return {'size': size, 'descriptor': descriptor,
            'missing_cells': missing, 'one_cell_wall': one_cell_wall}


def _calibrate_square_irregular(region, units, helix_data, raw_value):
    topology = _square_irregular_descriptor(region, units, helix_data)
    if topology is None:
        return None
    size = topology['size']
    candidates = [node for node in MISSING_HELIX_NODES if node[1] == size]
    if not candidates:
        return None
    target = topology['descriptor']
    scales = []
    for index in range(len(target)):
        values = [node[3][index] for node in candidates]
        scales.append(max(max(values)-min(values), 1e-9))
    ranked = []
    for node in candidates:
        distance = math.sqrt(sum(
            ((target[index]-node[3][index])/scales[index])**2
            for index in range(len(target))))
        ranked.append((distance, node))
    ranked.sort(key=lambda item: item[0])
    exact = ranked[0][0] < 1e-6
    # One-cell-wall S6/S8 nodes are retained for exact recognition but are not
    # used to interpolate otherwise stable closed sections.
    excluded = {'S6-R3x3C', 'S6-R4x4C', 'S8-R5x5C', 'S8-R6x6C'}
    pool = ranked if exact else [item for item in ranked
                                 if item[1][0] not in excluded]
    chosen = pool[:min(4, len(pool))]
    if exact:
        # A directly simulated topology is not an interpolation problem.
        # Mixing it with nearby nodes moved its native anchor and could even
        # change the predicted handedness near zero.
        chosen = [ranked[0]]
        weights = [1.0]
    else:
        weights = [(distance+.15)**-2 for distance, unused in chosen]
        total = sum(weights)
        weights = [value/total for value in weights]
    profile = tuple(sum(weights[index]*chosen[index][1][4][field]
                        for index in range(len(chosen)))
                    for field in range(4))
    x0, t0, k_minus, k_plus = profile
    coefficient = k_minus if raw_value < x0 else k_plus
    value = t0+coefficient*(raw_value-x0)
    nearest_distance = chosen[0][0]
    extrapolated = nearest_distance > 1.5
    low_confidence_wall = bool(size >= 6 and topology['one_cell_wall'])
    metadata = {
        'irregular_calibrated': True,
        'irregular_calibration_exact_node': exact,
        'irregular_calibration_extrapolated': extrapolated,
        'irregular_one_cell_wall': low_confidence_wall,
        'irregular_outer_size': size,
        'irregular_descriptor': target,
        'irregular_profile': profile,
        'irregular_neighbors': tuple(
            (item[1][0], weights[index], item[0])
            for index, item in enumerate(chosen)),
        'irregular_validation_rmse': MISSING_HELIX_RMSE}
    if exact and chosen[0][1][0] == 'S8-R4x4C':
        metadata.update({
            'irregular_calibration_version':
                S8_R4X4C_CALIBRATION_VERSION,
            'irregular_calibration_valid_indel_range':
                S8_R4X4C_VALIDATED_INDEL_RANGE,
            'irregular_calibration_failed_indel_points':
                S8_R4X4C_FAILED_INDEL_POINTS,
            'irregular_calibration_points': S8_R4X4C_TWIST_POINTS,
            'irregular_calibration_branch_rmse':
                dict(S8_R4X4C_BRANCH_RMSE),
            'irregular_failed_points_excluded_from_fit': True,
        })
    return value, metadata


def _log_interpolate(x, xs, ys):
    """Interpolate against log(J), clamping outside measured stiffness."""
    x = max(1e-12, float(x))
    if x <= xs[0]:
        return float(ys[0])
    if x >= xs[-1]:
        return float(ys[-1])
    lx = math.log(x)
    for left in range(len(xs)-1):
        if xs[left] <= x <= xs[left+1]:
            fraction = ((lx-math.log(xs[left])) /
                        (math.log(xs[left+1])-math.log(xs[left])))
            return ys[left] + fraction * (ys[left+1]-ys[left])
    return float(ys[-1])


def _cross_section_polar_moment(units):
    """Return the selected bundle's polar second moment in nm^4.

    cadnano geometry is displayed at 2.8 nm helix-centre spacing.  The
    calibration paper used an approximately 2.0 nm packed-cylinder spacing,
    so only this stiffness calculation is rescaled; the visible geometry is
    deliberately unchanged.
    """
    if not units:
        return 0.0
    scale = STIFFNESS_SPACING_NM / CENTER_SPACING_NM
    points = [(unit['coord'][0]*scale, unit['coord'][1]*scale)
              for unit in units]
    cx = sum(point[0] for point in points) / float(len(points))
    cy = sum(point[1] for point in points) / float(len(points))
    radius = 1.0
    intrinsic = math.pi * radius**4 / 2.0
    area = math.pi * radius**2
    return sum(intrinsic + area*((x-cx)**2 + (y-cy)**2)
               for x, y in points)


def _connectivity_fraction(region, helix_data):
    """Fraction of selected helices in the largest crossover component."""
    selected = set(int(value) for value in region['helices'])
    if len(selected) <= 1:
        return 1.0
    graph = dict((number, set()) for number in selected)
    for number in selected:
        for idx, partner in helix_data[number].get('crossovers', ()):
            partner = int(partner)
            if (region['start'] <= int(idx) <= region['end'] and
                    partner in selected and partner != number):
                graph[number].add(partner)
                graph[partner].add(number)
    largest = 1
    unseen = set(selected)
    while unseen:
        stack = [unseen.pop()]
        component = set(stack)
        while stack:
            current = stack.pop()
            for neighbor in graph[current]:
                if neighbor not in component:
                    component.add(neighbor)
                    unseen.discard(neighbor)
                    stack.append(neighbor)
        largest = max(largest, len(component))
    return largest / float(len(selected))


def estimate_global_twist(region, helix_data, lattice_pitch,
                          extra_edits=None, extra_base_delta=None):
    """Estimate physical global twist for a selected bundle region.

    The result combines effective helical pitch, selected cross-section polar
    moment, and crossover connectivity.  For one solid, complete and regular
    W x L cross-section it applies the final native-anchored coarse-grained
    calibration. Square outer frames with missing helices use the frozen
    topology-aware missing-helix database; geometries outside either
    calibrated domain retain the raw elastic prediction and are labelled as
    uncalibrated.
    """
    region = validate_regions([region])[0]
    units = _region_units(region, helix_data)
    turns = sum(int(unit.get(
        'native_bases', region['end']-region['start']+1))
        for unit in units) / float(lattice_pitch)
    current_bases = sum(unit['actual'] for unit in units)
    if extra_base_delta is None:
        extra_base_delta = sum(int(edit['length']) for edit in
                               (extra_edits or ())
                               if int(edit['helix']) in region['helices'] and
                               region['start'] <= int(edit['idx']) <= region['end'])
    final_bases = current_bases + int(extra_base_delta or 0)
    effective_pitch = final_bases / float(turns)
    polar_moment = _cross_section_polar_moment(units)
    mechanical_region = dict(region, helices=[unit['helix'] for unit in units])
    connectivity = _connectivity_fraction(mechanical_region, helix_data)
    lattice = _lattice_name(lattice_pitch)
    if lattice == 'square':
        at_1033 = _log_interpolate(polar_moment, _SQUARE_J,
                                   _SQUARE_TWIST_1033)
        at_1067 = _log_interpolate(polar_moment, _SQUARE_J,
                                   _SQUARE_TWIST_1067)
        fraction = (effective_pitch-10.33) / (10.67-10.33)
        raw_twist_per_base = at_1033 + fraction*(at_1067-at_1033)
        in_pitch_range = 10.2 <= effective_pitch <= 10.8
        in_j_range = _SQUARE_J[0] <= polar_moment <= _SQUARE_J[-1]
    else:
        reference = _log_interpolate(polar_moment, _HONEYCOMB_J,
                                     _HONEYCOMB_TWIST_1050)
        # Pitch sensitivity is directly constrained at J≈4027 and tends to
        # zero for the thickest measured bundle.  The thin-bundle endpoint is
        # anchored by the reported ~0.9°/bp 6HB rotation at 10.5 bp/turn.
        sensitivity = _log_interpolate(
            polar_moment, (85.0, 4027.0, 57000.0), (1.8, 1.15238, 0.0))
        raw_twist_per_base = reference + sensitivity*(effective_pitch-10.5)
        in_pitch_range = 10.0 <= effective_pitch <= 10.75
        in_j_range = _HONEYCOMB_J[0] <= polar_moment <= _HONEYCOMB_J[-1]
    section = classify_regular_cross_section(units, lattice)
    uncalibrated_twist_per_base = raw_twist_per_base * connectivity
    calibration_metadata = {}
    calibration_applicable = bool(section and connectivity >= .999)
    irregular_result = None
    if calibration_applicable:
        calibrated_connected, calibration_metadata = (
            _calibrate_cross_section_twist(lattice, raw_twist_per_base,
                                           section))
        twist_per_base = calibrated_connected * connectivity
    elif lattice == 'square' and connectivity >= .999:
        irregular_result = _calibrate_square_irregular(
            region, units, helix_data, raw_twist_per_base)
        if irregular_result is not None:
            calibrated_connected, calibration_metadata = irregular_result
            twist_per_base = calibrated_connected * connectivity
            calibration_applicable = True
        else:
            twist_per_base = uncalibrated_twist_per_base
    else:
        twist_per_base = uncalibrated_twist_per_base
    native_bases = sum(int(unit.get(
        'native_bases', region['end']-region['start']+1)) for unit in units)
    mean_indel_per_helix = (
        final_bases-native_bases) / float(len(units))
    calibrated_domain = calibration_metadata.get(
        'irregular_calibration_valid_indel_range')
    calibration_domain_exceeded = bool(
        calibrated_domain and not (
            float(calibrated_domain[0])-1e-9 <= mean_indel_per_helix <=
            float(calibrated_domain[1])+1e-9))
    if calibrated_domain:
        calibration_metadata['irregular_calibration_domain_exceeded'] = \
            calibration_domain_exceeded
    # The effective pitch above intentionally uses the indel-adjusted
    # ``final_bases``.  The calibrated twist rate, however, is expressed per
    # nominal design base (the same convention used by the SNUPI calibration
    # measurements).  Integrating it over the indel-adjusted length would
    # count each insertion/deletion twice: once through effective_pitch and a
    # second time through the angle span.
    nominal_length = native_bases / float(len(units))
    actual_length = final_bases / float(len(units))
    total_twist = twist_per_base * nominal_length
    confidence = '中等'
    reasons = []
    if connectivity < .999:
        confidence = '低'
        reasons.append('选区 crossover 未完全连通（最大连通 %.0f%%）' %
                       (connectivity*100.0))
    if not in_pitch_range:
        confidence = '低'
        reasons.append('有效螺距超出文献标定区间')
    if not in_j_range:
        confidence = '低'
        reasons.append('截面刚度超出文献标定区间')
    if len(units) < 2:
        confidence = '低'
        reasons.append('单 helix 不适用整体束扭转标定')
    if lattice == 'honeycomb' and confidence != '低':
        reasons.append('Honeycomb 对 crossover 拓扑较敏感')
    if calibration_applicable and irregular_result is None:
        node = ('精确节点' if calibration_metadata['calibration_exact_node']
                else 'W/L 参数插值')
        confidence = ('高' if calibration_metadata['calibration_exact_node']
                      else '中高')
        if calibration_metadata['calibration_extrapolated']:
            node = '校准边界外推'
            confidence = '低'
        reasons.append(
            '已应用实心规则截面 W%d×L%d 的粗粒化校正 %s（%s）' %
            (section['width'], section['layers'],
             SNUPI_CALIBRATION_VERSION, node))
    elif irregular_result is not None:
        node = ('直接模拟拓扑节点' if
                calibration_metadata['irregular_calibration_exact_node']
                else '拓扑筛选参数插值')
        confidence = ('高' if
                      calibration_metadata['irregular_calibration_exact_node']
                      else '中等')
        if calibration_metadata['irregular_calibration_extrapolated']:
            confidence = '低'
            node = '拓扑校准域外推'
        if calibration_metadata['irregular_one_cell_wall']:
            confidence = '低'
            reasons.append('S6/S8 单 helix 壁厚截面柔性高，粗粒化模拟稳定性不足')
        if calibration_domain_exceeded:
            confidence = '低'
            reasons.append(
                'S8-R4x4C平均indel/helix超出收敛校准域−10到+6；'
                '+8和+10 insertion未收敛且未参与拟合')
        reasons.append(
            '已应用 Square 非实心截面 S%d 的%s（验证 RMSE %.3f°/base）' %
            (calibration_metadata['irregular_outer_size'], node,
             calibration_metadata['irregular_validation_rmse']))
    else:
        confidence = '低'
        if section is None:
            reasons.append('没有适用的实心或 Square 非实心校正节点；保留原始弹性预测')
        elif connectivity < .999:
            reasons.append('截面未完全连通；保留原始弹性预测')
    result = {
        'lattice': lattice, 'effective_pitch': effective_pitch,
        'polar_moment_nm4': polar_moment,
        'connectivity_fraction': connectivity,
        'uncalibrated_twist_per_base_deg': uncalibrated_twist_per_base,
        'twist_per_base_deg': twist_per_base,
        'total_twist_deg': total_twist,
        # ``length_bp`` is retained as the downstream integration length for
        # saved tasks, Remove/Add Twist, the 3D preview and reports.
        'length_bp': nominal_length,
        'nominal_length_bp': nominal_length,
        'actual_length_bp': actual_length,
        'twist_integration_length_basis': 'nominal',
        'mean_indel_per_helix': mean_indel_per_helix,
        'handedness': ('右手' if total_twist > 1e-9 else
                       '左手' if total_twist < -1e-9 else '近似无扭转'),
        'confidence': confidence,
        'calibration_version': SNUPI_CALIBRATION_VERSION,
        'calibrated': calibration_applicable,
        'cross_section_calibration_applicable': calibration_applicable,
        'cross_section_width': section['width'] if section else None,
        'cross_section_layers': section['layers'] if section else None,
        'cross_section_solid_regular': bool(section),
        'note': '；'.join(reasons)}
    result.update(calibration_metadata)
    return result


def estimate_elastic_bend(region, helix_data, per_helix_counts,
                          lattice='square'):
    """Predict equilibrium bend from indel eigenstrain using elastic rods.

    This is a compact CanDo-inspired cross-sectional energy minimization. Each
    dsDNA contributes the published axial and bending stiffness; crossover
    connectivity controls how much axial mismatch is transmitted through the
    bundle.  The raw elastic angle is then mapped to SNUPI Static by the final
    lattice-specific calibration.  It is intentionally fast and is not a
    replacement for full FE.
    """
    region = validate_regions([region])[0]
    units = _region_units(region, helix_data)
    length_nm = (sum(unit['actual'] for unit in units) /
                 float(len(units)) * RISE_NM)
    cx = sum(unit['coord'][0] for unit in units)/float(len(units))
    cy = sum(unit['coord'][1] for unit in units)/float(len(units))
    rows = []
    for unit in units:
        delta_nm = (float(per_helix_counts.get(unit['helix'], 0)) *
                    RISE_NM)
        rows.append((unit['coord'][0]-cx, unit['coord'][1]-cy,
                     delta_nm/max(RISE_NM, length_nm)))
    mean_strain = sum(row[2] for row in rows)/float(len(rows))
    connectivity = _connectivity_fraction(region, helix_data)
    axial = DSDNA_STRETCH_PN * max(.05, connectivity)
    bending = DSDNA_BEND_PN_NM2 * len(rows)
    sxx = sum(row[0]*row[0] for row in rows)
    syy = sum(row[1]*row[1] for row in rows)
    sxy = sum(row[0]*row[1] for row in rows)
    bx = -axial*sum(row[0]*(row[2]-mean_strain) for row in rows)
    by = -axial*sum(row[1]*(row[2]-mean_strain) for row in rows)
    a11, a22, a12 = axial*sxx+bending, axial*syy+bending, axial*sxy
    determinant = a11*a22-a12*a12
    if abs(determinant) < 1e-12:
        kx, ky = 0.0, 0.0
    else:
        kx = (bx*a22-by*a12)/determinant
        ky = (by*a11-bx*a12)/determinant
    curvature = math.sqrt(kx*kx+ky*ky)
    raw_angle = math.degrees(curvature*length_nm)
    lattice = _lattice_name(lattice)
    angle = _snupi_calibrate('bend', lattice, raw_angle,
                             preserve_zero=True)
    curvature = math.radians(angle) / length_nm if length_nm > 1e-12 else 0.0
    direction = (math.degrees(math.atan2(-ky, kx)) % 360.0
                 if curvature > 1e-12 else 0.0)
    confidence = '中等'
    calibration = SNUPI_CALIBRATION['bend'][lattice]
    reasons = ['快速弹性杆模型 + SNUPI Static 校正 %s（n=%d，RMSE %.3f°）' %
               (SNUPI_CALIBRATION_VERSION, calibration['count'],
                calibration['rmse'])]
    if connectivity < .999:
        confidence = '低'
        reasons.append('crossover 最大连通分量 %.0f%%' %
                       (connectivity*100.0))
    if angle > 120.0:
        confidence = '低'
        reasons.append('超过文献中几何/弹性模型较可靠的约120°范围')
    return {'angle_degrees': angle,
            'uncalibrated_angle_degrees': raw_angle,
            'curvature_per_nm': curvature,
            'radius_nm': (1.0/curvature if curvature > 1e-12 else None),
            'direction_degrees': direction, 'length_nm': length_nm,
            'connectivity_fraction': connectivity,
            'confidence': confidence, 'note': '；'.join(reasons),
            'lattice': lattice,
            'calibration_version': SNUPI_CALIBRATION_VERSION,
            'calibration_intercept': calibration['intercept'],
            'calibration_slope': calibration['slope'],
            'calibrated': True,
            'stretch_stiffness_pn': DSDNA_STRETCH_PN,
            'bend_stiffness_pn_nm2': DSDNA_BEND_PN_NM2}


def estimate_existing_bend(region, helix_data, lattice='square'):
    """Infer one net elastic bend from indels already present in a JSON.

    This is intentionally opt-in.  Existing indels may encode twist, length
    compensation or local repair rather than a designed bend, so the result
    is an approximate net curvature over the selected region—not recovered
    design intent.
    """
    region = validate_regions([region])[0]
    units = _region_units(region, helix_data)
    counts = dict((unit['helix'], unit['existing_delta']) for unit in units)
    result = estimate_elastic_bend(region, helix_data, counts, lattice)
    differential = (max(counts.values())-min(counts.values())
                    if counts else 0)
    result.update({
        'per_helix_existing_delta': counts,
        'differential_indel_span': differential,
        'source': 'existing_indels',
        'approximate': True})
    prefix = ('从已有 indel 反推净曲率；无法区分 Bend、Twist、长度补偿或局部修复。')
    if differential == 0:
        prefix += ' 各 helix 的净 indel 相同，因此未检测到截面差异弯曲。'
    result['note'] = prefix + ' ' + result['note']
    if result['confidence'] != '低':
        result['confidence'] = '中等'
    return result


def _edits_for_counts(units, counts, length):
    edits = []
    for unit, count in zip(units, counts):
        for idx in _even_positions(unit['start'], unit['end'], count,
                                   unit['allowed']):
            edits.append({'helix': unit['helix'], 'idx': idx,
                          'length': int(length)})
    return edits


def _preferred_existing_removals(units, direction, target_magnitude):
    """Choose whole existing indels to undo before creating new ones.

    ``direction`` is the desired signed base change: -1 removes right-handed
    twist and +1 removes left-handed twist.  Existing insertions therefore
    qualify for -1, while existing deletions qualify for +1.  cadnano stores
    one insertion object per site, so an existing multi-base insertion is
    removed as one indivisible edit instead of being silently overwritten.
    """
    target_magnitude = max(0, int(target_magnitude))
    candidates = []
    for unit_index, unit in enumerate(units):
        for idx, original_length in sorted(unit.get('existing', {}).items()):
            original_length = int(original_length)
            delta = -original_length
            if delta and (1 if delta > 0 else -1) == int(direction):
                candidates.append((abs(delta), unit_index, int(idx),
                                   original_length))

    # Bounded subset selection maximizes the amount handled by undoing
    # existing indels without overshooting the desired correction.  Typical
    # regions contain few existing indels; limiting sums to the target keeps
    # the live Apply calculation compact even for long designs.
    states = {0: ()}
    for magnitude, unit_index, idx, original_length in candidates:
        for subtotal, chosen in sorted(list(states.items()), reverse=True):
            updated = subtotal + magnitude
            if updated <= target_magnitude and updated not in states:
                states[updated] = chosen + ((magnitude, unit_index, idx,
                                             original_length),)
    removed_magnitude = max(states)
    chosen = states[removed_magnitude]
    edits = []
    per_helix_delta = dict((unit['helix'], 0) for unit in units)
    for unused_magnitude, unit_index, idx, original_length in chosen:
        unit = units[unit_index]
        delta = -int(original_length)
        edits.append({'helix': unit['helix'], 'idx': idx,
                      'length': delta, 'operation': 'remove_existing',
                      'original_length': int(original_length)})
        per_helix_delta[unit['helix']] += delta
    return edits, per_helix_delta, removed_magnitude


def plan_remove_twist(region, helix_data, lattice_pitch,
                      measured_twist_per_base=None):
    """Plan bidirectional indel edits that minimize macroscopic twist.

    A right-handed baseline first undoes existing insertions, then adds any
    still-required deletions.  A left-handed baseline performs the symmetric
    operation: undo existing deletions first, then add insertions.
    """
    region = validate_regions([region])[0]
    units = _region_units(region, helix_data)
    domain_size = _twist_domain_size(lattice_pitch)
    # The lattice fixes a continuous amount of rotation over the selected
    # index span.  Do not round that rotation to a whole number of turns:
    # doing so creates discontinuities where extending the region by one base
    # can suddenly *reduce* the required number of deletions.  Existing indels
    # change the current physical base count, but not the native lattice
    # rotation represented by the selected indices.
    turns = [int(unit.get(
        'native_bases', region['end']-region['start']+1)) /
        float(lattice_pitch) for unit in units]
    current = sum(unit['actual'] for unit in units)
    # Estimate the signed integer base change that brings the empirical global
    # twist closest to zero.  Positive is insertion, negative is deletion.
    # A local finite difference supplies the first estimate; the full-range
    # sign-crossing search below handles asymmetric calibrated branches and
    # integer rounding at interpolation boundaries.
    baseline = estimate_global_twist(region, helix_data, lattice_pitch,
                                     extra_base_delta=0)
    baseline_value = (baseline['twist_per_base_deg']
                      if measured_twist_per_base is None else
                      float(measured_twist_per_base))
    direction = (-1 if baseline_value > 1e-12 else
                 1 if baseline_value < -1e-12 else 0)
    layer_delta = direction*max(1, len(units))
    one_layer = estimate_global_twist(
        region, helix_data, lattice_pitch,
        extra_base_delta=layer_delta) if direction else baseline
    slope = ((one_layer['twist_per_base_deg']-
              baseline['twist_per_base_deg']) / float(layer_delta)
             if direction else 0.0)
    estimated_delta = (int(round(-baseline_value / slope))
                       if abs(slope) > 1e-12 else 0)
    if direction and estimated_delta*direction < 0:
        estimated_delta = 0
    matching_existing = sum(
        abs(int(length)) for unit in units
        for length in unit.get('existing', {}).values()
        if int(length) and
        (1 if -int(length) > 0 else -1) == direction)
    new_capacity = sum(len(unit['allowed']) for unit in units)
    capacity = matching_existing + new_capacity
    estimated_magnitude = max(0, min(capacity, abs(estimated_delta)))
    candidate_magnitudes = set([0, estimated_magnitude, capacity])
    # The calibrated response is locally linear, so only inspect a compact
    # integer neighborhood.  This keeps live updates responsive for 64+ helices.
    radius = max(4, min(12, len(units)//4 + 2))
    candidate_magnitudes.update(
        max(0, min(capacity, estimated_magnitude+offset))
        for offset in range(-radius, radius+1))

    residual_cache = {}
    def residual_for_delta(delta):
        delta = int(delta)
        if delta not in residual_cache:
            model_value = estimate_global_twist(
                region, helix_data, lattice_pitch,
                extra_base_delta=delta)['twist_per_base_deg']
            residual_cache[delta] = (
                baseline_value + model_value-baseline['twist_per_base_deg'])
        return residual_cache[delta]

    # The native-anchored calibration is piecewise linear at the native pitch,
    # and its deletion/insertion slopes can differ.  A local derivative alone
    # may therefore miss the correct integer solution when the correction
    # crosses that anchor.  Locate the sign crossing over the full allowed
    # range, then inspect its compact integer neighborhood.
    if direction and capacity:
        end_residual = residual_for_delta(direction*capacity)
        if baseline_value*end_residual <= 0:
            low, high = 0, capacity
            while high-low > 1:
                middle = (low+high)//2
                middle_residual = residual_for_delta(direction*middle)
                if baseline_value*middle_residual > 0:
                    low = middle
                else:
                    high = middle
            root_magnitude = min(
                (low, high), key=lambda magnitude:
                abs(residual_for_delta(direction*magnitude)))
            candidate_magnitudes.update(
                max(0, min(capacity, root_magnitude+offset))
                for offset in range(-radius, radius+1))

    # Whole existing indels are indivisible.  Keep only magnitudes that can be
    # completed with the available new one-base edits.
    feasible = []
    removal_cache = {}
    placement_cache = {}
    for magnitude in candidate_magnitudes:
        removal = _preferred_existing_removals(
            units, direction, magnitude) if direction else ([], {}, 0)
        removal_cache[magnitude] = removal
        remaining_for_candidate = magnitude-removal[2]
        grouped = _edit_adjustments_by_helix(removal[0])
        capacities = [_twist_unit_capacity(
            unit, direction, domain_size,
            grouped.get(int(unit['helix']), ())) for unit in units]
        if remaining_for_candidate <= sum(capacities):
            try:
                candidate_counts = _allocate_evenly_across_units(
                    remaining_for_candidate, capacities)
                candidate_edits, candidate_metadata = \
                    _domain_aware_twist_edits(
                        units, candidate_counts, direction, domain_size,
                        prior_edits=removal[0])
            except TwistBendError:
                continue
            feasible.append(magnitude)
            placement_cache[magnitude] = (
                candidate_counts, candidate_edits, candidate_metadata)
    if not feasible:
        raise TwistBendError(
            '现有 indel、±3/domain 硬限制和安全位点共同约束下，无法生成 Remove Twist 方案。')
    selected_magnitude = min(feasible, key=lambda magnitude: (
        abs(residual_for_delta(direction*magnitude)), magnitude))
    existing_edits, per_helix_delta, removed_magnitude = removal_cache[
        selected_magnitude]
    remaining = selected_magnitude-removed_magnitude
    counts, new_edits, domain_metadata = placement_cache[selected_magnitude]
    for edit in new_edits:
        edit['operation'] = 'add'
        per_helix_delta[edit['helix']] += int(edit['length'])
    edits = existing_edits + new_edits
    actual_delta = sum(int(edit['length']) for edit in edits)
    final_bases = current + actual_delta
    achieved = final_bases / float(sum(turns))
    prediction = estimate_global_twist(region, helix_data, lattice_pitch,
                                       extra_base_delta=actual_delta)
    if measured_twist_per_base is not None:
        prediction = dict(prediction)
        residual = residual_for_delta(actual_delta)
        prediction['model_twist_per_base_deg'] = prediction[
            'twist_per_base_deg']
        prediction['twist_per_base_deg'] = residual
        prediction['total_twist_deg'] = residual*prediction['length_bp']
        prediction['handedness'] = (
            '右手' if residual > 1e-9 else
            '左手' if residual < -1e-9 else '近似无扭转')
        prediction['measurement_anchored'] = True
        prediction['measured_baseline_twist_per_base_deg'] = baseline_value
        prediction['note'] = (
            '基线来自上传模拟文件；indel 响应增量来自当前弹性/校准模型。' +
            (('；' + prediction.get('note', ''))
             if prediction.get('note') else ''))
    result = {'kind': 'remove_twist',
            'edits': edits,
            'turns': sum(turns), 'actual_bases': current,
            'achieved_pitch': achieved,
            'target_pitch': TARGET_PITCH,
            'net_base_delta': actual_delta,
            'removed_insertions': sum(
                edit.get('operation') == 'remove_existing' and
                edit.get('original_length', 0) > 0 for edit in edits),
            'removed_deletions': sum(
                edit.get('operation') == 'remove_existing' and
                edit.get('original_length', 0) < 0 for edit in edits),
            'added_insertions': sum(
                edit.get('operation') != 'remove_existing' and
                edit['length'] > 0 for edit in edits),
            'added_deletions': sum(
                edit.get('operation') != 'remove_existing' and
                edit['length'] < 0 for edit in edits),
            'baseline_twist_per_base_deg': baseline_value,
            'baseline_source': ('simulation_measurement'
                                if measured_twist_per_base is not None
                                else 'elastic_prediction'),
            'twist_prediction': prediction,
            # Preserve the historical sign convention: positive means the
            # deletion-like correction used for a right-handed baseline.
            'per_helix_counts': dict(
                (unit['helix'], -per_helix_delta.get(unit['helix'], 0))
                for unit in units)}
    result.update(domain_metadata)
    return result


def plan_add_twist(region, helix_data, lattice_pitch, angle_degrees,
                   handedness='right', indels_per_helix=None):
    """Plan a cumulative relative twist using balanced integer indels."""
    region = validate_regions([region])[0]
    units = _region_units(region, helix_data)
    domain_size = _twist_domain_size(lattice_pitch)
    sign = 1 if str(handedness).lower().startswith('r') else -1
    requested = float(angle_degrees) * sign
    # In this editor a positive/right-handed extra rotation is produced by
    # adding bases; a negative/left-handed rotation uses deletions.
    # One deleted/inserted nucleotide changes the ensemble mismatch by
    # 1/TARGET_PITCH turn.  The visible bundle rotation is the mean mismatch
    # across the selected helices.
    if indels_per_helix is None:
        count = int(round(abs(requested) / 360.0 * len(units) *
                          TARGET_PITCH))
    else:
        requested_average = float(indels_per_helix)
        requested = abs(float(angle_degrees)) * (1 if requested_average >= 0 else -1)
        count = int(round(abs(requested_average) * len(units)))
    edit_length = 1 if requested >= 0 else -1
    capacities = [_twist_unit_capacity(
        unit, edit_length, domain_size) for unit in units]
    counts = _allocate_evenly_across_units(count, capacities)
    edits, domain_metadata = _domain_aware_twist_edits(
        units, counts, edit_length, domain_size)
    achieved = (sum(counts) / float(max(1, len(units)) * TARGET_PITCH) *
                360.0)
    if edit_length < 0:
        achieved *= -1
    turns = sum(max(1, unit['actual'] / float(lattice_pitch))
                for unit in units)
    prediction = estimate_global_twist(region, helix_data, lattice_pitch,
                                       extra_edits=edits)
    result = {'kind': 'add_twist', 'edits': edits,
            'requested_angle': requested, 'achieved_angle': achieved,
            'degrees_per_100bp': achieved * 100.0 /
                max(1.0, sum(unit.get(
                    'native_bases', region['end']-region['start']+1)
                    for unit in units) / float(len(units))),
            'lattice_turns': turns,
            'twist_prediction': prediction,
            'per_helix_counts': dict((unit['helix'], count)
                                     for unit, count in zip(units, counts))}
    result.update(domain_metadata)
    return result


def plan_bend(region, helix_data, angle_degrees, direction_degrees,
              material_twist_degrees=0.0, incremental=False,
              elastic_compensation=False, lattice='square'):
    """Plan symmetric bend indels using actual lattice coordinates."""
    region = validate_regions([region])[0]
    units = _region_units(region, helix_data)
    theta = math.radians(float(angle_degrees))
    direction = math.radians(float(direction_degrees))
    # The selector uses Qt screen coordinates. Therefore 90° (arrow upward)
    # corresponds to negative lattice Y, not positive mathematical Y.
    ux, uy = math.cos(direction), -math.sin(direction)
    # If the input bundle already carries material twist, the helix occupying
    # the inside/outside of a global bend changes continuously along the
    # region.  Average the rotating cross-section projection; this is the
    # design-space counterpart of measuring twist in a transported Bishop
    # frame rather than mistaking bend-frame rotation for material twist.
    material_twist = math.radians(float(material_twist_degrees))
    samples = 33
    projections = []
    for unit in units:
        x, y = unit['coord']
        values = []
        for rank in range(samples):
            phase = material_twist * rank / float(samples-1)
            ca, sa = math.cos(phase), math.sin(phase)
            rotated_x, rotated_y = x*ca-y*sa, x*sa+y*ca
            values.append(rotated_x*ux + rotated_y*uy)
        projections.append(sum(values)/float(len(values)))
    neutral = (max(projections) + min(projections)) / 2.0
    span_actual = sum(unit['actual'] for unit in units) / float(len(units))
    inner_extent = max(projections) - neutral
    maximum = 720.0
    if abs(theta) > 1e-9:
        radius = span_actual * RISE_NM / abs(theta)
        minimum_radius = inner_extent + 0.7
        if radius <= minimum_radius:
            maximum = math.degrees(
                span_actual * RISE_NM / max(0.7, minimum_radius)) * 0.92
            raise TwistBendError(
                '该选区的弯曲角过大，会使内侧半径反转。当前最大安全角约为 %.1f°；'
                '请减小角度或增加选区长度。' % maximum)
    def counts_for(design_angle):
        design_theta = math.radians(float(design_angle))
        raw = [-design_theta * (value-neutral) / RISE_NM
               for value in projections]
        return [int(round(value)) -
                (0 if incremental else unit['existing_delta'])
                for value, unit in zip(raw, units)]

    design_angle = float(angle_degrees)
    signed_counts = counts_for(design_angle)
    elastic = estimate_elastic_bend(
        region, helix_data,
        dict((unit['helix'], count)
             for unit, count in zip(units, signed_counts)), lattice)
    if elastic_compensation and abs(float(angle_degrees)) > 1e-9:
        # Integer indels and intrinsic dsDNA bending stiffness make the first
        # geometric plan undershoot. Iteratively pre-compensate the design
        # strain so the elastic equilibrium approaches the requested angle.
        best = (abs(elastic['angle_degrees']-abs(float(angle_degrees))),
                design_angle, list(signed_counts), dict(elastic))
        for unused in range(6):
            achieved = elastic['angle_degrees']
            if achieved <= 1e-9:
                break
            ratio = abs(float(angle_degrees))/achieved
            if abs(ratio-1.0) < .01:
                break
            design_angle = min(maximum*.995,
                               max(0.0, design_angle*ratio))
            signed_counts = counts_for(design_angle)
            elastic = estimate_elastic_bend(
                region, helix_data,
                dict((unit['helix'], count)
                     for unit, count in zip(units, signed_counts)), lattice)
            candidate = (abs(elastic['angle_degrees']-
                             abs(float(angle_degrees))),
                         design_angle, list(signed_counts), dict(elastic))
            if candidate[0] < best[0]:
                best = candidate
        unused_error, design_angle, signed_counts, elastic = best
    # Bending used to choose globally even allowed indices with
    # ``_even_positions``.  It now uses the exact same equal-partition,
    # absolute-domain and nearest-safe-site resolver as Add/Remove Twist.
    # Existing indels remain immutable baseline strain and all prior
    # crossover/nick/deletion protections are retained by the shared context.
    domain_size = 7 if _lattice_name(lattice) == 'honeycomb' else 8
    edits, domain_metadata = _domain_aware_signed_edits(
        units, signed_counts, domain_size)
    radius = elastic['radius_nm']
    result = {'kind': 'bend', 'edits': edits,
            'requested_angle': float(angle_degrees),
            'geometric_design_angle': design_angle,
            'elastic_prediction': elastic,
            'direction': float(direction_degrees) % 360.0,
            'material_twist_degrees': float(material_twist_degrees),
            'radius_nm': radius, 'neutral_axis': neutral,
            'lattice': _lattice_name(lattice),
            'per_helix_counts': dict((unit['helix'], count)
                                     for unit, count in zip(units,
                                                            signed_counts))}
    result.update(domain_metadata)
    optimize_bend_plan_curvature(region, helix_data, result,
                                 maximum_passes=3)
    return result


def merge_plans(plans):
    """Combine plans and reject conflicting edits at the same base."""
    edits = {}
    for plan_index, plan in enumerate(plans):
        for edit in plan.get('edits', ()):
            key = (int(edit['helix']), int(edit['idx']))
            old = edits.get(key)
            if old is not None and old['length'] != int(edit['length']):
                raise TwistBendError(
                    '任务 %d 在 helix %d[%d] 同时要求 insertion 和 deletion。' %
                    (plan_index + 1, key[0], key[1]))
            edits[key] = dict(edit)
    return sorted(edits.values(), key=lambda item: (item['helix'], item['idx']))

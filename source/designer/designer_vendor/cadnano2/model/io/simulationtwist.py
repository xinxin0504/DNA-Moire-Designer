"""Measure material twist and centreline bend from labelled simulation XYZ.

The supported XYZ representation contains one contiguous coordinate block per
double-stranded helix.  The first line is the total point count; every later
line is ``label x y z``.  Block order follows the active caDNAno helix order.
Nominal base indices are reconstructed from the current design, including
insertions and deletions.  Material twist is measured in a parallel-transport
(Bishop) frame, so centreline bending is not counted as twist.
"""

from __future__ import division

import math
import os
import shlex


class SimulationAnalysisError(ValueError):
    pass


def _add(a, b):
    return tuple(a[index] + b[index] for index in range(3))


def _sub(a, b):
    return tuple(a[index] - b[index] for index in range(3))


def _scale(value, factor):
    return tuple(component * factor for component in value)


def _dot(a, b):
    return sum(a[index] * b[index] for index in range(3))


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _norm(value):
    return math.sqrt(max(0.0, _dot(value, value)))


def _unit(value, fallback=(1.0, 0.0, 0.0)):
    length = _norm(value)
    return (_scale(value, 1.0/length) if length > 1e-12 else fallback)


def _mean(values):
    if not values:
        raise SimulationAnalysisError('没有可计算的坐标。')
    return tuple(sum(value[index] for value in values)/float(len(values))
                 for index in range(3))


def _rotate_minimal(vector, first_axis, second_axis):
    """Parallel-transport ``vector`` through the minimal tangent rotation."""
    first_axis, second_axis = _unit(first_axis), _unit(second_axis)
    cosine = max(-1.0, min(1.0, _dot(first_axis, second_axis)))
    rotation_axis = _cross(first_axis, second_axis)
    sine = _norm(rotation_axis)
    if sine <= 1e-12:
        if cosine >= 0.0:
            return vector
        seed = (1.0, 0.0, 0.0)
        if abs(_dot(seed, first_axis)) > .9:
            seed = (0.0, 1.0, 0.0)
        rotation_axis = _unit(_cross(first_axis, seed))
        return _sub(_scale(rotation_axis, 2.0*_dot(rotation_axis, vector)),
                    vector)
    rotation_axis = _scale(rotation_axis, 1.0/sine)
    # Rodrigues formula, with cos(theta)=cosine and sin(theta)=sine.
    return _add(_add(_scale(vector, cosine),
                     _scale(_cross(rotation_axis, vector), sine)),
                _scale(rotation_axis,
                       _dot(rotation_axis, vector)*(1.0-cosine)))


def _unwrap(values):
    result = []
    for value in values:
        if result:
            while value-result[-1] > math.pi:
                value -= 2.0*math.pi
            while value-result[-1] < -math.pi:
                value += 2.0*math.pi
        result.append(value)
    return result


def _linear_fit(x_values, y_values):
    x_mean = sum(x_values)/float(len(x_values))
    y_mean = sum(y_values)/float(len(y_values))
    denominator = sum((value-x_mean)**2 for value in x_values)
    if denominator <= 1e-12:
        raise SimulationAnalysisError('选区过短，无法拟合 twist。')
    slope = sum((x-x_mean)*(y-y_mean)
                for x, y in zip(x_values, y_values))/denominator
    intercept = y_mean-slope*x_mean
    residual = [y-(slope*x+intercept)
                for x, y in zip(x_values, y_values)]
    rms = math.sqrt(sum(value*value for value in residual)/len(residual))
    return slope, intercept, rms


def read_labelled_xyz(path):
    if os.path.splitext(path)[1].lower() != '.xyz':
        raise SimulationAnalysisError(
            '当前精确映射首先支持带 helix 分块标签的 XYZ。PDB/CIF 或 '
            'DAT+TOP 需要同时提供 nucleotide mapping，不能仅凭原子/链数量猜测。')
    blocks, labels = [], []
    current_label = None
    with open(path, 'r', encoding='utf-8', errors='replace') as source:
        try:
            declared = int(source.readline().strip())
        except (TypeError, ValueError):
            raise SimulationAnalysisError('XYZ 第一行必须是坐标点总数。')
        for line in source:
            fields = line.split()
            if len(fields) < 4:
                continue
            try:
                point = tuple(float(value) for value in fields[-3:])
            except ValueError:
                continue
            label = fields[0]
            if label != current_label:
                labels.append(label)
                blocks.append([])
                current_label = label
            blocks[-1].append(point)
    count = sum(len(block) for block in blocks)
    if count != declared:
        raise SimulationAnalysisError(
            'XYZ 点数不一致：文件声明 %d，实际读取 %d。' %
            (declared, count))
    if not blocks:
        raise SimulationAnalysisError('XYZ 中没有可读取的 helix 坐标块。')
    return labels, blocks


def _covered_indices(intervals):
    covered = set()
    for low, high in intervals:
        covered.update(range(int(low), int(high)+1))
    return covered


def _double_stranded_mapping(unit):
    paired = (_covered_indices(unit.get('scaffold_intervals', ())) &
              _covered_indices(unit.get('staple_intervals', ())))
    mapping = []
    for base in sorted(paired):
        multiplicity = max(0, 1+int(unit.get('insertions', {}).get(base, 0)))
        for sub_index in range(multiplicity):
            mapping.append(base + sub_index/float(max(1, multiplicity)))
    return mapping


def _interpolate(parameters, points, target):
    if target < parameters[0]-1e-9 or target > parameters[-1]+1e-9:
        raise SimulationAnalysisError(
            '选区 base %g 超出模拟文件的双链坐标覆盖范围。' % target)
    for index, value in enumerate(parameters):
        if abs(value-target) <= 1e-9:
            return points[index]
        if value > target:
            first, second = parameters[index-1], value
            weight = (target-first)/float(second-first)
            return _add(_scale(points[index-1], 1.0-weight),
                        _scale(points[index], weight))
    return points[-1]


def _map_xyz_blocks(helix_data, labels, blocks, selected):
    active = []
    for number, unit in helix_data.items():
        mapping = _double_stranded_mapping(unit)
        if mapping:
            active.append((int(number), mapping))
    selected_active = [item for item in active if item[0] in set(selected)]
    if len(blocks) == len(active):
        targets = active
    elif len(blocks) == len(selected_active):
        targets = selected_active
    else:
        raise SimulationAnalysisError(
            'XYZ 含 %d 个 helix 块；当前设计有 %d 个双链 helix，选区有 %d 个。'
            '无法唯一映射。' % (len(blocks), len(active),
                              len(selected_active)))
    mapped = {}
    for (number, parameters), label, points in zip(targets, labels, blocks):
        if len(parameters) != len(points):
            raise SimulationAnalysisError(
                'helix %d（XYZ 标签 %s）双链碱基数不匹配：设计 %d，XYZ %d。' %
                (number, label, len(parameters), len(points)))
        mapped[number] = (parameters, points, label)
    missing = sorted(set(selected)-set(mapped))
    if missing:
        raise SimulationAnalysisError(
            '模拟文件缺少选中 helix：%s。' %
            ', '.join(str(value) for value in missing))
    return mapped


def _coordinate_scale_to_nm(mapped, source_type):
    """Return the physical nm represented by one input coordinate unit."""
    if source_type in ('cadnano_pdb', 'cadnano_mmcif'):
        # PDB/mmCIF readers already convert Angstrom coordinates to nm.
        return 1.0
    if source_type == 'oxdna_dat_top':
        # Standard oxDNA length unit.
        return .8518
    steps = []
    for unused_number, (unused_parameters, points, unused_label) in \
            mapped.items():
        for first, second in zip(points, points[1:]):
            distance = _norm(_sub(second, first))
            if distance > 1e-9:
                steps.append(distance)
    if not steps:
        return 1.0
    steps.sort()
    median = steps[len(steps)//2]
    # Labelled SNUPI XYZ files are normally Angstrom (about 3.4 per bp).
    # Accept an already-nm XYZ as well instead of relying on its filename.
    return .1 if median > 1.5 else 1.0


def _scaled_mapping(mapped, scale):
    if abs(float(scale)-1.0) <= 1e-12:
        return mapped
    return dict((number, (parameters,
                          [_scale(point, scale) for point in points], label))
                for number, (parameters, points, label) in mapped.items())


def _cross_section_gradient(coords, values):
    """Fit one affine indel field and return its cross-section span."""
    count = float(len(coords))
    mean_x = sum(point[0] for point in coords)/count
    mean_y = sum(point[1] for point in coords)/count
    mean_value = sum(values)/count
    centered = [(point[0]-mean_x, point[1]-mean_y)
                for point in coords]
    shifted = [float(value)-mean_value for value in values]
    cxx = sum(point[0]*point[0] for point in centered)
    cxy = sum(point[0]*point[1] for point in centered)
    cyy = sum(point[1]*point[1] for point in centered)
    bx = sum(point[0]*value for point, value in zip(centered, shifted))
    by = sum(point[1]*value for point, value in zip(centered, shifted))
    # A tiny ridge also supports one-row/one-column sections without changing
    # the well-conditioned 2-D solution in any meaningful way.
    ridge = max(1e-12, (cxx+cyy)*1e-10)
    cxx, cyy = cxx+ridge, cyy+ridge
    determinant = cxx*cyy-cxy*cxy
    gx = (bx*cyy-by*cxy)/determinant
    gy = (by*cxx-bx*cxy)/determinant
    fitted = [gx*point[0]+gy*point[1] for point in centered]
    span = max(fitted)-min(fitted) if fitted else 0.0
    total_variance = sum(value*value for value in shifted)
    residual_variance = sum((value-fit)**2
                            for value, fit in zip(shifted, fitted))
    explained = (max(0.0, 1.0-residual_variance/total_variance)
                 if total_variance > 1e-12 else 1.0)
    return {'gradient': (gx, gy), 'span_bases': span,
            'explained_fraction': explained,
            'residual_rms_bases': math.sqrt(
                residual_variance/max(1.0, count))}


def classify_design_bending(helix_data, selected_helices, start, end):
    """Classify programmed Bend from the JSON indel field.

    A common-mode indel profile changes pitch/length but cannot bend the
    bundle.  Bend is the cross-sectional differential component.  Mixed-sign
    indels are treated as explicit Bend intent; a one-sign pattern must also
    exhibit a resolvable cross-section gradient.  Nonuniform but radially or
    checkerboard-symmetric patterns are labelled complex rather than silently
    interpreted as a single Bend direction.
    """
    selected = sorted(set(int(value) for value in selected_helices))
    start, end = sorted((int(start), int(end)))
    events = {}
    positive = negative = 0
    for column, number in enumerate(selected):
        for idx, length in helix_data[number].get('insertions', {}).items():
            idx, length = int(idx), int(length)
            if not start <= idx <= end or not length:
                continue
            events.setdefault(idx, []).append((column, length))
            positive += int(length > 0)
            negative += int(length < 0)
    if not events:
        return {'classification': 'none', 'apply_bishop_correction': False,
                'reason': '选区没有 insertion/deletion；无设计性 bending。',
                'positive_sites': 0, 'negative_sites': 0,
                'maximum_gradient_span_bases': 0.0,
                'maximum_gradient_explained_fraction': 1.0}

    coords = [helix_data[number]['coord'] for number in selected]
    cumulative = [0] * len(selected)
    all_common_mode = True
    maximum_span = 0.0
    maximum_explained = 0.0
    maximum_residual = 0.0
    for idx in sorted(events):
        for column, length in events[idx]:
            cumulative[column] += length
        common = len(set(cumulative)) == 1
        all_common_mode = all_common_mode and common
        if common:
            continue
        fit = _cross_section_gradient(coords, cumulative)
        if fit['span_bases'] > maximum_span:
            maximum_span = fit['span_bases']
            maximum_explained = fit['explained_fraction']
            maximum_residual = fit['residual_rms_bases']

    uniform_total = len(set(cumulative)) == 1
    if all_common_mode:
        classification = 'none'
        apply = False
        reason = ('所有 helix 的 indel 数量与轴向分布完全相同；这是均匀螺距/'
                  '长度变化，不产生设计性 bending。')
    elif not (positive and negative) and uniform_total:
        classification = 'none'
        apply = False
        reason = ('所有 helix 的同符号 indel 总量相同；轴向位点虽可错开，但'
                  '不存在整体截面长度梯度，因此不判定为全局 bending。')
    elif positive and negative:
        classification = 'bending'
        apply = True
        reason = ('选区同时包含 insertion 和 deletion，存在截面差分长度；'
                  '按设计性 bending 处理。')
    elif maximum_span >= .5 and maximum_explained >= .25:
        classification = 'bending'
        apply = True
        reason = ('单一符号 indel 在截面上明显不对称（梯度跨度 %.3f base）；'
                  '按设计性 bending 处理。' % maximum_span)
    else:
        classification = 'complex'
        apply = False
        reason = ('indel 分布不均匀，但不能由单一截面梯度可靠解释；标记为复杂'
                  '局部变形，不自动执行 Bishop bending 扣除。')
    return {
        'classification': classification,
        'apply_bishop_correction': apply,
        'reason': reason, 'positive_sites': positive,
        'negative_sites': negative,
        'maximum_gradient_span_bases': maximum_span,
        'maximum_gradient_explained_fraction': maximum_explained,
        'maximum_gradient_residual_rms_bases': maximum_residual}


def _principal_axis(points):
    """Best-fit line direction through a sequence of 3-D points."""
    center = _mean(points)
    shifted = [_sub(point, center) for point in points]
    covariance = [[sum(point[row]*point[column] for point in shifted)
                   for column in range(3)] for row in range(3)]
    axis = _unit(_sub(points[-1], points[0]), (0.0, 0.0, 1.0))
    for unused in range(32):
        updated = tuple(sum(covariance[row][column]*axis[column]
                            for column in range(3)) for row in range(3))
        axis = _unit(updated, axis)
    if _dot(axis, _sub(points[-1], points[0])) < 0.0:
        axis = _scale(axis, -1.0)
    return axis


def _fit_cross_section_phases(bases, planes, centers, reference,
                              first_axes, second_axes):
    phases = []
    for plane, center, first_axis, second_axis in zip(
            planes, centers, first_axes, second_axes):
        projected = []
        for point in plane:
            vector = _sub(point, center)
            projected.append((_dot(vector, first_axis),
                              _dot(vector, second_axis)))
        a = sum(q[0]*p[0]+q[1]*p[1]
                for q, p in zip(reference, projected))
        b = sum(q[0]*p[1]-q[1]*p[0]
                for q, p in zip(reference, projected))
        phases.append(math.atan2(b, a))
    degrees = [math.degrees(value) for value in _unwrap(phases)]
    slope, unused_intercept, fit_rms = _linear_fit(bases, degrees)
    return float(slope), float(fit_rms)


def _analyze_mapped(path, source_type, mapping_text, mapped, helix_data,
                    selected_helices, start, end):
    """Measure a uniquely mapped set of helix centreline coordinates."""
    coordinate_scale_nm = _coordinate_scale_to_nm(mapped, source_type)
    mapped = _scaled_mapping(mapped, coordinate_scale_nm)
    selected = sorted(set(int(value) for value in selected_helices))
    if len(selected) < 2:
        raise SimulationAnalysisError('精确 Twist 至少需要两个 helix。')
    start, end = sorted((int(start), int(end)))
    if end-start < 2:
        raise SimulationAnalysisError('精确 Twist 至少需要三个 base 平面。')
    bases = list(range(start, end+1))
    planes = []
    for base in bases:
        planes.append([
            _interpolate(mapped[number][0], mapped[number][1], float(base))
            for number in selected])
    centers = [_mean(plane) for plane in planes]
    # SNUPI centre lines contain base-scale corrugation and thermal/static
    # relaxation noise.  Differentiating adjacent planes converts that small
    # zig-zag into a rapidly rotating frame and falsely subtracts genuine
    # material Twist.  Estimate the macroscopic tangent over an axial window;
    # the transported frame still follows real Bend, while suppressing the
    # high-frequency component that is not bundle-centreline curvature.
    tangent_half_window = max(2, min(24, len(centers)//6))
    tangents = []
    for index in range(len(centers)):
        low = max(0, index-tangent_half_window)
        high = min(len(centers)-1, index+tangent_half_window)
        delta = _sub(centers[high], centers[low])
        tangents.append(_unit(delta, tangents[-1] if tangents else
                              (1.0, 0.0, 0.0)))

    reference = []
    cx = sum(helix_data[number]['coord'][0] for number in selected)/len(selected)
    cy = sum(helix_data[number]['coord'][1] for number in selected)/len(selected)
    for number in selected:
        reference.append((helix_data[number]['coord'][0]-cx,
                          helix_data[number]['coord'][1]-cy))

    first_vectors = [_sub(point, centers[0]) for point in planes[0]]
    global_axis = _principal_axis(centers)
    global_first = None
    for vector in first_vectors:
        projected = _sub(vector, _scale(global_axis,
                                        _dot(vector, global_axis)))
        if _norm(projected) > 1e-9:
            global_first = _unit(projected)
            break
    if global_first is None:
        raise SimulationAnalysisError('首个截面退化，无法建立截面坐标架。')
    global_second = _unit(_cross(global_axis, global_first))
    global_slope, global_fit_rms = _fit_cross_section_phases(
        bases, planes, centers, reference,
        [global_first]*len(bases), [global_second]*len(bases))

    bishop_first = None
    for vector in first_vectors:
        projected = _sub(vector, _scale(tangents[0],
                                        _dot(vector, tangents[0])))
        if _norm(projected) > 1e-9:
            bishop_first = _unit(projected)
            break
    if bishop_first is None:
        raise SimulationAnalysisError('首个截面退化，无法建立 Bishop frame。')
    bishop_second = _unit(_cross(tangents[0], bishop_first))
    bishop_first_axes, bishop_second_axes = [], []
    transported_first = bishop_first
    transported_second = bishop_second
    previous_tangent = tangents[0]
    for index, tangent in enumerate(tangents):
        if index:
            transported_first = _rotate_minimal(
                transported_first, previous_tangent, tangent)
            transported_first = _unit(_sub(
                transported_first,
                _scale(tangent, _dot(transported_first, tangent))))
            transported_second = _unit(_cross(tangent, transported_first))
        bishop_first_axes.append(transported_first)
        bishop_second_axes.append(transported_second)
        previous_tangent = tangent
    bishop_slope, bishop_fit_rms = _fit_cross_section_phases(
        bases, planes, centers, reference,
        bishop_first_axes, bishop_second_axes)
    design_bend = classify_design_bending(
        helix_data, selected, start, end)
    if design_bend['apply_bishop_correction']:
        slope, fit_rms = bishop_slope, bishop_fit_rms
        method = ('signed all-plane fit in a smoothed parallel-transport '
                  'Bishop frame; enabled by the JSON indel bend classifier')
    else:
        slope, fit_rms = global_slope, global_fit_rms
        method = ('signed all-plane fit about the global best-fit axis; '
                  'Bishop bend subtraction disabled by the JSON indel '
                  'classifier')
    first_window = max(2, min(8, len(tangents)//4))
    first_tangent = _unit(_mean(tangents[:first_window]))
    last_tangent = _unit(_mean(tangents[-first_window:]))
    bend_angle = math.degrees(math.acos(max(-1.0, min(
        1.0, _dot(first_tangent, last_tangent)))))
    tangent_change = _sub(last_tangent, first_tangent)
    bend_direction = (math.degrees(math.atan2(
        _dot(tangent_change, bishop_second),
        _dot(tangent_change, bishop_first))) % 360.0
        if bend_angle > 1e-9 else 0.0)
    arc_length = sum(_norm(_sub(centers[index+1], centers[index]))
                     for index in range(len(centers)-1))
    radius = (arc_length/math.radians(bend_angle)
              if bend_angle > 1e-9 else None)
    return {
        'source_path': os.path.abspath(path),
        'source_type': source_type,
        'helices': selected, 'start': start, 'end': end,
        'twist_per_base_deg': float(slope),
        'total_twist_deg': float(slope*(end-start)),
        'twist_fit_rms_deg': float(fit_rms),
        'global_axis_twist_per_base_deg': float(global_slope),
        'global_axis_twist_fit_rms_deg': float(global_fit_rms),
        'bishop_twist_per_base_deg': float(bishop_slope),
        'bishop_twist_fit_rms_deg': float(bishop_fit_rms),
        'bishop_correction_applied': bool(
            design_bend['apply_bishop_correction']),
        'design_bend_classification': design_bend,
        'bend_angle_deg': float(bend_angle),
        'bend_direction_deg': float(bend_direction),
        'bend_radius_nm': radius,
        'centerline_arc_length_nm': float(arc_length),
        'bend_radius_coordinate_units': (
            radius/coordinate_scale_nm if radius is not None else None),
        'centerline_arc_length_coordinate_units': float(
            arc_length/coordinate_scale_nm),
        'coordinate_scale_nm': float(coordinate_scale_nm),
        'tangent_smoothing_half_window_bases': tangent_half_window,
        'method': method,
        'mapping': mapping_text,
        'preview_centerlines': [
            [number, [[base, planes[rank][column][0],
                       planes[rank][column][1], planes[rank][column][2]]
                      for rank, base in enumerate(bases)]]
            for column, number in enumerate(selected)],
    }


def analyze_labelled_xyz(path, helix_data, selected_helices, start, end):
    """Return bend-corrected material twist and centreline bend metrics."""
    selected = sorted(set(int(value) for value in selected_helices))
    labels, blocks = read_labelled_xyz(path)
    mapped = _map_xyz_blocks(helix_data, labels, blocks, selected)
    return _analyze_mapped(
        path, 'labelled_xyz',
        'current caDNAno double-stranded bases, indel-aware',
        mapped, helix_data, selected, start, end)


def _expected_export_records(document):
    """Rebuild the exact nucleotide ordering used by cadnano exports."""
    try:
        from .oxdnaexport import _collect, _number_records
    except (ImportError, ValueError) as error:
        raise SimulationAnalysisError(
            '无法读取当前 cadnano 导出器：%s。' % error)
    unused_records, strands, unused_assigned, unused_residual = (
        _collect(document, 2.8))
    if not strands:
        raise SimulationAnalysisError('当前设计没有可映射的 DNA 链。')
    ordered = _number_records(strands)
    return ordered, strands


def _effective_parameter(record):
    count = max(1, int(record.get('count', 1)))
    # The simulation coordinate itself follows strand direction, but the
    # caDNAno mapping key is the nominal base plus insertion sub-index for
    # both complementary strands.
    return int(record['idx']) + float(record.get('sub', 0))/count


def _minimum_image(point, reference, box):
    if not box:
        return point
    values = list(point)
    for axis in range(3):
        length = float(box[axis])
        if length > 1e-12:
            values[axis] -= round((values[axis]-reference[axis])/length)*length
    return tuple(values)


def _mapped_from_nucleotide_points(records_and_points, helix_data, selected,
                                   box=None):
    buckets = {}
    for record, point in records_and_points:
        number = int(record['strand'].virtualHelix().number())
        if number not in selected:
            continue
        parameter = _effective_parameter(record)
        buckets.setdefault((number, round(parameter, 9)), []).append(point)
    mapped = {}
    for number in selected:
        target = _double_stranded_mapping(helix_data[number])
        parameters, points = [], []
        for parameter in target:
            values = buckets.get((number, round(parameter, 9)))
            if not values:
                raise SimulationAnalysisError(
                    '模拟文件缺少 helix %d base %.6g 的坐标，无法严格映射。' %
                    (number, parameter))
            parameters.append(parameter)
            aligned = [values[0]] + [
                _minimum_image(value, values[0], box)
                for value in values[1:]]
            points.append(_mean(aligned))
        if box and points:
            unwrapped = [points[0]]
            for point in points[1:]:
                unwrapped.append(_minimum_image(point, unwrapped[-1], box))
            points = unwrapped
        mapped[number] = (parameters, points, 'H%d' % number)
    return mapped


def _read_oxdna_dat(path):
    points = []
    with open(path, 'r', encoding='utf-8', errors='replace') as source:
        lines = source.readlines()
    if len(lines) < 4 or not lines[0].lstrip().startswith('t'):
        raise SimulationAnalysisError('DAT 文件缺少标准 oxDNA 三行文件头。')
    try:
        box = tuple(float(value) for value in
                    lines[1].split('=', 1)[1].split()[:3])
    except (IndexError, ValueError):
        raise SimulationAnalysisError('DAT 第二行缺少有效的周期盒尺寸。')
    for line in lines[3:]:
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            points.append(tuple(float(value) for value in fields[:3]))
        except ValueError:
            raise SimulationAnalysisError('DAT 中存在无法读取的坐标行。')
    return points, box


def _validate_top(path, nucleotide_count, strand_count):
    root, unused = os.path.splitext(path)
    top_path = root + '.top'
    if not os.path.exists(top_path):
        raise SimulationAnalysisError(
            '请选择与同名 .top 位于同一文件夹的 oxDNA .dat 文件。')
    with open(top_path, 'r', encoding='utf-8', errors='replace') as source:
        fields = source.readline().split()
    try:
        observed = tuple(int(value) for value in fields[:2])
    except (TypeError, ValueError):
        raise SimulationAnalysisError('TOP 第一行不是有效的 nucleotide/strand 数。')
    expected = (int(nucleotide_count), int(strand_count))
    if observed != expected:
        raise SimulationAnalysisError(
            'DAT+TOP 与当前 JSON 不匹配：TOP 为 %d nt/%d strands，当前设计为 '
            '%d nt/%d strands。' % (observed + expected))
    return top_path


def _read_pdb_residues(path):
    chains, chain_points = [], []
    residue_points, current_residue, current_chain = [], None, None
    marker = False
    with open(path, 'r', encoding='utf-8', errors='replace') as source:
        for line in source:
            if 'CADNANO OXDNA-BACKMAPPED' in line:
                marker = True
            if line[:6].strip() == 'TER':
                if residue_points:
                    chain_points.append(_mean(residue_points))
                    residue_points = []
                if chain_points:
                    chains.append(chain_points)
                    chain_points = []
                current_residue, current_chain = None, None
                continue
            if line[:6].strip() not in ('ATOM', 'HETATM'):
                continue
            key = (line[21:22], line[22:27].strip())
            try:
                point = (float(line[30:38])/10.0,
                         float(line[38:46])/10.0,
                         float(line[46:54])/10.0)
            except ValueError:
                continue
            if current_chain is not None and key[0] != current_chain:
                if residue_points:
                    chain_points.append(_mean(residue_points))
                    residue_points = []
                if chain_points:
                    chains.append(chain_points)
                    chain_points = []
                current_residue = None
            if current_residue is not None and key != current_residue:
                chain_points.append(_mean(residue_points))
                residue_points = []
            current_chain, current_residue = key[0], key
            residue_points.append(point)
    if residue_points:
        chain_points.append(_mean(residue_points))
    if chain_points:
        chains.append(chain_points)
    if not chains:
        raise SimulationAnalysisError('PDB 中没有可读取的 ATOM/HETATM 坐标。')
    return chains, marker


def _read_mmcif_residues(path):
    rows, columns, in_atom_loop = [], [], False
    marker = False
    with open(path, 'r', encoding='utf-8', errors='replace') as source:
        for raw in source:
            line = raw.strip()
            if 'cadnano oxDNA-backmapped' in line:
                marker = True
            if line == 'loop_':
                columns, in_atom_loop = [], False
                continue
            if line.startswith('_atom_site.'):
                columns.append(line.split('.', 1)[1])
                in_atom_loop = True
                continue
            if in_atom_loop and line and not line.startswith(('_', '#')):
                fields = shlex.split(line)
                if len(fields) >= len(columns):
                    rows.append(dict(zip(columns, fields)))
            elif in_atom_loop and line.startswith('#'):
                in_atom_loop = False
    needed = {'label_asym_id', 'label_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'}
    if not rows or not needed <= set(rows[0]):
        raise SimulationAnalysisError('mmCIF 中没有完整的 _atom_site 坐标。')
    residues, order = {}, []
    for row in rows:
        key = (row['label_asym_id'], row['label_seq_id'])
        try:
            point = (float(row['Cartn_x'])/10.0,
                     float(row['Cartn_y'])/10.0,
                     float(row['Cartn_z'])/10.0)
        except ValueError:
            continue
        if key not in residues:
            residues[key] = []
            order.append(key)
        residues[key].append(point)
    chains, current = [], None
    for key in order:
        if key[0] != current:
            chains.append([])
            current = key[0]
        chains[-1].append(_mean(residues[key]))
    return chains, marker


def _match_structure_chains(chains, strands, marker):
    expected = [item['records'] for item in strands]
    expected_lengths = [len(item) for item in expected]
    observed_lengths = [len(item) for item in chains]
    if observed_lengths == expected_lengths and marker:
        return [(record, point) for records, points in zip(expected, chains)
                for record, point in zip(records, points)]
    # Generic structures may be accepted only when chain lengths define one
    # unique assignment. Duplicate lengths are intentionally rejected rather
    # than silently swapping indistinguishable strands/helices.
    assignments, used = [], set()
    for points in chains:
        matches = [index for index, records in enumerate(expected)
                   if index not in used and len(records) == len(points)]
        if len(matches) != 1:
            raise SimulationAnalysisError(
                '外部 PDB/mmCIF 的链无法与当前 JSON 唯一映射（链长 %d 有 %d 个'
                '候选）。请使用本软件导出的结构文件或带 helix 标签的 XYZ。' %
                (len(points), len(matches)))
        index = matches[0]
        used.add(index)
        assignments.extend(zip(expected[index], points))
    if len(used) != len(expected):
        raise SimulationAnalysisError(
            'PDB/mmCIF 链数与当前 JSON 不一致，无法完整映射。')
    return assignments


def analyze_simulation_file(path, document, helix_data, selected_helices,
                            start, end):
    """Analyze labelled XYZ, exact cadnano DAT+TOP, or mappable PDB/mmCIF."""
    extension = os.path.splitext(path)[1].lower()
    if extension == '.xyz':
        return analyze_labelled_xyz(
            path, helix_data, selected_helices, start, end)
    selected = sorted(set(int(value) for value in selected_helices))
    ordered, strands = _expected_export_records(document)
    if extension == '.dat':
        points, box = _read_oxdna_dat(path)
        _validate_top(path, len(ordered), len(strands))
        if len(points) != len(ordered):
            raise SimulationAnalysisError(
                'DAT 坐标数 %d 与当前 JSON nucleotide 数 %d 不一致。' %
                (len(points), len(ordered)))
        pairs = zip(ordered, points)
        source_type = 'oxdna_dat_top'
        mapping_text = 'exact current cadnano oxDNA export order, TOP-validated'
    elif extension in ('.pdb', '.cif', '.mmcif'):
        box = None
        if extension == '.pdb':
            chains, marker = _read_pdb_residues(path)
        else:
            chains, marker = _read_mmcif_residues(path)
        pairs = _match_structure_chains(chains, strands, marker)
        source_type = 'cadnano_pdb' if extension == '.pdb' else 'cadnano_mmcif'
        mapping_text = ('exact cadnano exporter chain/residue order' if marker
                        else 'unique external chain-length assignment')
    else:
        raise SimulationAnalysisError(
            '仅支持 .xyz、.dat + 同名 .top、.pdb、.cif/.mmcif。')
    mapped = _mapped_from_nucleotide_points(
        pairs, helix_data, set(selected), box)
    return _analyze_mapped(path, source_type, mapping_text, mapped,
                           helix_data, selected, start, end)

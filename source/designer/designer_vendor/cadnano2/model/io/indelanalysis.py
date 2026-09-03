"""Shared indel-distribution analysis and pair-aware curvature refinement."""

import base64
import collections
import csv
import gzip
import html
import json
import math
import os
import statistics


def _longitudinal(record, helix, base):
    if record == [-1, -1, -1, -1]:
        return False
    for offset in (0, 2):
        partner, partner_base = map(int, record[offset:offset + 2])
        if partner != helix or abs(partner_base - base) != 1:
            return False
    return True


def _geometry_centers(design):
    encoded = (design.get('curved_metadata', {}) or {}).get('geometry_data')
    if not encoded:
        return {}
    geometry = json.loads(gzip.decompress(base64.b64decode(encoded)))
    positions = collections.defaultdict(list)
    for key, frame in geometry.get('frames', {}).items():
        fields = key.split(':')
        if len(fields) < 2:
            continue
        strand_type, helix = fields[:2]
        if strand_type != 'scaffold':
            continue
        x, y, z = frame['pos']
        positions[int(helix)].append((math.hypot(x, y), float(z)))
    return dict((helix, (statistics.mean(p[0] for p in values),
                         statistics.mean(p[1] for p in values)))
                for helix, values in positions.items())


def _cluster(values):
    groups = []
    for value in sorted(values):
        if not groups or value - groups[-1][-1] > 1.01:
            groups.append([value])
        else:
            groups[-1].append(value)
    return [sum(group) / float(len(group)) for group in groups]


def _crossover_events(rows):
    by_pair = collections.defaultdict(lambda: collections.defaultdict(set))
    for strand_type in ('scaf', 'stap'):
        edges = set()
        for helix, row in rows.items():
            for base, entry in enumerate(row.get(strand_type, [])):
                for offset in (0, 2):
                    other, other_base = map(int, entry[offset:offset + 2])
                    if other >= 0 and other != helix:
                        edges.add(tuple(sorted(((helix, base),
                                                (other, other_base)))))
        for left, right in edges:
            pair = tuple(sorted((left[0], right[0])))
            by_pair[pair][strand_type].add(
                (left[1] + right[1]) / 2.0)
    return by_pair


def _physical_pair_rows(design, rows):
    centers = _geometry_centers(design)
    if len(centers) < 2:
        return []
    events = _crossover_events(rows)
    result = []
    for pair in sorted(events):
        if pair[0] not in centers or pair[1] not in centers:
            continue
        # An AutoCS crossover is the authoritative lattice/topology proof
        # that two helices are physical neighbours.  Do not re-filter that
        # topology by their *deformed* geometry: after a Curved/Frame target
        # geometry is embedded, a legitimate neighbour (for example H0-H1
        # at one Frame face) can be slightly farther apart than the global
        # nearest pair and was previously omitted from both optimization and
        # the final report.
        sites = sorted(set(_cluster(events[pair]['scaf']) +
                           _cluster(events[pair]['stap'])))
        if len(sites) < 2:
            continue
        inner, outer = sorted(pair, key=lambda helix: centers[helix][0])
        result.append({'pair': pair, 'inner': inner, 'outer': outer,
                       'events': sites})
    return result


def _arc_index(events, position):
    for index, start in enumerate(events):
        end = events[(index + 1) % len(events)]
        inside = (start < position < end if end > start else
                  position > start or position < end)
        if inside:
            return index
    return None


def _circular_gap_cv(positions, length):
    ordered = sorted(positions)
    if len(ordered) < 2:
        return 0.0
    gaps = [right-left for left, right in zip(ordered, ordered[1:])]
    gaps.append(float(length)-ordered[-1]+ordered[0])
    mean = statistics.mean(gaps)
    return statistics.pstdev(gaps) / mean if mean else 0.0


def _interval_bin_weights(events, length, bin_width):
    bin_count = max(1, int(math.ceil(length / float(bin_width))))
    weights = [[0.0] * bin_count for unused in events]
    for base in range(length):
        interval = _arc_index(events, base + .314159)
        if interval is not None:
            weights[interval][base // bin_width] += 1.0
    return weights


def _short_staple_deletion_protection(rows):
    nodes = set()
    for helix, row in rows.items():
        for base, entry in enumerate(row.get('stap', [])):
            if entry != [-1, -1, -1, -1]:
                nodes.add((helix, base))
    unseen = set(nodes)
    protected = collections.defaultdict(set)
    while unseen:
        seed = unseen.pop()
        component = {seed}
        stack = [seed]
        while stack:
            helix, base = stack.pop()
            entry = rows[helix]['stap'][base]
            for offset in (0, 2):
                neighbour = (int(entry[offset]), int(entry[offset + 1]))
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        actual = 0
        for helix, base in component:
            loop = int(rows[helix].get('loop', [])[base])
            skip = int(rows[helix].get('skip', [])[base])
            actual += max(0, 1 + max(0, loop) - max(0, -skip))
        if actual <= 21:
            for helix, base in component:
                protected[helix].add(base)
    return protected


def short_staple_deletion_protection(rows):
    """Return bases whose staple component cannot tolerate a deletion.

    This is the public entry point shared by Curved/Frame placement passes.
    Keeping the component walk in one implementation prevents the initial
    Frame allocation and its later pair-aware refinement from applying
    different 21-nt staple protections.
    """
    return _short_staple_deletion_protection(rows)


def _pair_analysis(pair_rows, variables, length, bin_width=42):
    by_helix = collections.defaultdict(list)
    for variable in variables:
        by_helix[variable['helix']].append(variable)
    rows = []
    for pair_row in pair_rows:
        diffs = [0] * len(pair_row['events'])
        for variable in by_helix[pair_row['inner']]:
            interval = _arc_index(pair_row['events'], variable['position'])
            if interval is not None:
                diffs[interval] -= variable['value']
        for variable in by_helix[pair_row['outer']]:
            interval = _arc_index(pair_row['events'], variable['position'])
            if interval is not None:
                diffs[interval] += variable['value']
        total = sum(diffs)
        # Only adjacent helix pairs with a designed differential strain are
        # curvature-bearing pairs.  Uniform pairs are intentionally omitted.
        if total == 0:
            continue
        inner, outer = pair_row['inner'], pair_row['outer']
        if total < 0:
            inner, outer = outer, inner
            diffs = [-value for value in diffs]
            total = -total
        floor_value, remainder = divmod(total, len(diffs))
        rows.append({
            'pair': pair_row['pair'], 'inner': inner, 'outer': outer,
            'events': list(pair_row['events']), 'differences': diffs,
            'total_difference': total,
            'mean_difference': total / float(len(diffs)),
            'ideal_floor': floor_value,
            'ideal_ceiling': floor_value + (1 if remainder else 0)})

    reverse = outside = zero = above = severity = 0
    for row in rows:
        floor_value, ceiling = row['ideal_floor'], row['ideal_ceiling']
        for difference in row['differences']:
            reverse += difference < 0
            if difference < floor_value or difference > ceiling:
                outside += 1
                severity += min(abs(difference-floor_value),
                                abs(difference-ceiling)) ** 2
            zero += difference == 0 and floor_value >= 1
            above += difference > ceiling

    bin_count = max(1, int(math.ceil(length / float(bin_width))))
    bins = [0.0] * bin_count
    weights = [0.0] * bin_count
    for base in range(length):
        sample = base + .314159
        for row in rows:
            interval = _arc_index(row['events'], sample)
            if interval is None or row['mean_difference'] <= 0:
                continue
            bin_index = base // bin_width
            bins[bin_index] += (row['differences'][interval] /
                                row['mean_difference'])
            weights[bin_index] += 1.0
    normalized = [value/weight if weight else 0.0
                  for value, weight in zip(bins, weights)]
    positive = [value for value in normalized if value > 0]
    axial_cv = (statistics.pstdev(positive) / statistics.mean(positive)
                if len(positive) > 1 else 0.0)
    interval_count = sum(len(row['differences']) for row in rows)
    return rows, {
        'physical_adjacent_curvature_pairs': len(rows),
        'crossover_intervals': interval_count,
        'reverse_curvature_intervals': reverse,
        'outside_floor_ceiling_intervals': outside,
        # This is a diagnostic subset of outside, not an extra error class.
        'zero_curvature_intervals_where_floor_is_one': zero,
        'above_ceiling_intervals': above,
        'outside_severity': severity,
        'normalized_axial_curvature_bin_width_bp': int(bin_width),
        'normalized_axial_curvature_42bp_bin_cv': axial_cv,
        'normalized_axial_curvature_42bp_bin_minimum': min(positive or [0]),
        'normalized_axial_curvature_42bp_bin_maximum': max(positive or [0]),
        'normalized_axial_curvature_42bp_bins': normalized}


def optimize_curved_pair_curvature(design, curvature_indels,
                                   maximum_passes=10):
    """Refine Curved indel coordinates without changing design topology.

    The implementation is incremental: candidate coordinates are shortlisted
    by their vector of physical helix-pair crossover intervals, and a trial
    move updates only the affected interval counters.  This preserves the
    strict priority ordering without the multi-minute full-design recompute
    used by the experimental prototype.
    """
    rows = {int(row['num']): row for row in design.get('vstrands', [])}
    length = int(design.get('num_bases') or max(
        [len(row.get('loop', [])) for row in rows.values()] or [0]))
    pair_rows = _physical_pair_rows(design, rows)
    records = list((curvature_indels or {}).get('rings', []))
    variables = []
    for record in records:
        helix = int(record['helix'])
        for item in record.get('indel_placements', []):
            variables.append({
                'helix': helix, 'value': int(item['value']),
                'position': int(item['base']),
                'partition': int(item['partition']),
                'partition_count': int(item['partition_count']),
                'ideal': float(item['ideal_base']),
                'partition_start': int(item['partition_start']),
                'partition_end': int(item['partition_end']),
                'target_domain': int(item['target_domain']),
                'domain_size': int(record.get('domain_size_bp', 7))})
    if not variables or not pair_rows:
        unused_rows, summary = _pair_analysis(
            pair_rows, variables, max(1, length))
        return summary

    deletion_protected = _short_staple_deletion_protection(rows)
    safe = {}
    for helix, row in rows.items():
        safe[helix] = [base for base in range(1, length-1)
                       if _longitudinal(row['scaf'][base], helix, base) and
                       _longitudinal(row['stap'][base], helix, base)]

    adjacent_rows = collections.defaultdict(list)
    for row_index, pair_row in enumerate(pair_rows):
        adjacent_rows[pair_row['inner']].append((row_index, -1))
        adjacent_rows[pair_row['outer']].append((row_index, 1))

    # First establish the designed curvature direction for every physical
    # pair.  Uniform-strain pairs are not bending pairs and are omitted.
    raw_differences = [[0] * len(row['events']) for row in pair_rows]
    for variable in variables:
        for row_index, sign in adjacent_rows[variable['helix']]:
            interval = _arc_index(pair_rows[row_index]['events'],
                                  variable['position'])
            if interval is not None:
                raw_differences[row_index][interval] += (
                    sign * variable['value'])
    active_indices = []
    direction = {}
    for row_index, differences in enumerate(raw_differences):
        total = sum(differences)
        if total:
            active_indices.append(row_index)
            direction[row_index] = 1 if total > 0 else -1
    if not active_indices:
        unused_rows, summary = _pair_analysis(
            pair_rows, variables, max(1, length))
        return summary

    # Remove inactive rows from the adjacency map so their crossover register
    # cannot influence either candidate signatures or the score.
    for helix in list(adjacent_rows):
        adjacent_rows[helix] = [
            (row_index, sign * direction[row_index])
            for row_index, sign in adjacent_rows[helix]
            if row_index in direction]

    for variable in variables:
        candidates = [
            base for base in safe.get(variable['helix'], [])
            if variable['partition_start'] <= base <=
            variable['partition_end'] and
            (variable['value'] > 0 or
             base not in deletion_protected[variable['helix']])]
        signatures = collections.defaultdict(list)
        for base in candidates:
            signature = tuple(
                (row_index,
                 _arc_index(pair_rows[row_index]['events'], base))
                for row_index, unused_sign in
                adjacent_rows[variable['helix']])
            signatures[signature].append(base)
        shortlisted = {variable['position']}
        for values in signatures.values():
            shortlisted.update(sorted(values, key=lambda base: (
                base // variable['domain_size'] !=
                variable['target_domain'],
                abs(base-variable['ideal']), base))[:2])
        variable['candidates'] = sorted(shortlisted)

    domain_counts = collections.Counter()
    by_helix = collections.defaultdict(list)
    for index, variable in enumerate(variables):
        domain = variable['position'] // variable['domain_size']
        domain_counts[(variable['helix'], domain)] += 1
        by_helix[variable['helix']].append(index)

    row_differences = {}
    for row_index in active_indices:
        row_differences[row_index] = [
            direction[row_index] * value
            for value in raw_differences[row_index]]
        pair_row = pair_rows[row_index]
        total = sum(row_differences[row_index])
        floor_value, remainder = divmod(
            total, len(row_differences[row_index]))
        pair_row['ideal_floor'] = floor_value
        pair_row['ideal_ceiling'] = floor_value + (1 if remainder else 0)
        pair_row['mean_difference'] = (
            total / float(len(row_differences[row_index])))
        pair_row['bin_weights'] = _interval_bin_weights(
            pair_row['events'], length, 42)

    helix_cvs = dict((helix, _circular_gap_cv(
        [variables[index]['position'] for index in indices], length))
        for helix, indices in by_helix.items())

    def evaluate():
        reverse = outside = severity = zero = above = 0
        bin_count = max(1, int(math.ceil(length / 42.0)))
        bins = [0.0] * bin_count
        weights = [0.0] * bin_count
        for row_index in active_indices:
            pair_row = pair_rows[row_index]
            floor_value = pair_row['ideal_floor']
            ceiling = pair_row['ideal_ceiling']
            mean = pair_row['mean_difference']
            for interval, difference in enumerate(
                    row_differences[row_index]):
                reverse += difference < 0
                if difference < floor_value or difference > ceiling:
                    outside += 1
                    severity += min(abs(difference-floor_value),
                                    abs(difference-ceiling)) ** 2
                zero += difference == 0 and floor_value >= 1
                above += difference > ceiling
                if mean > 0:
                    for bin_index, weight in enumerate(
                            pair_row['bin_weights'][interval]):
                        bins[bin_index] += weight * difference / mean
                        weights[bin_index] += weight
        normalized = [value/weight if weight else 0.0
                      for value, weight in zip(bins, weights)]
        positive = [value for value in normalized if value > 0]
        axial_cv = (statistics.pstdev(positive)/statistics.mean(positive)
                    if len(positive) > 1 else 0.0)
        single = statistics.mean(helix_cvs.values()) if helix_cvs else 0.0
        distance = statistics.mean(
            abs(variable['position']-variable['ideal'])
            for variable in variables)
        summary = {
            'physical_adjacent_curvature_pairs': len(active_indices),
            'crossover_intervals': sum(
                len(row_differences[index]) for index in active_indices),
            'reverse_curvature_intervals': reverse,
            'outside_floor_ceiling_intervals': outside,
            'zero_curvature_intervals_where_floor_is_one': zero,
            'above_ceiling_intervals': above,
            'outside_severity': severity,
            'normalized_axial_curvature_42bp_bin_cv': axial_cv,
            'normalized_axial_curvature_42bp_bin_minimum':
                min(positive or [0]),
            'normalized_axial_curvature_42bp_bin_maximum':
                max(positive or [0]),
            'normalized_axial_curvature_42bp_bins': normalized,
            'mean_single_helix_spacing_cv': single,
            'mean_distance_from_equal_partition_target_bp': distance}
        # Strict lexicographic priority confirmed for Curved Design:
        # reverse -> outside -> severity -> axial fluctuation -> single helix.
        score = (reverse, outside, severity, round(axial_cv, 9),
                 round(single, 9), round(distance, 6))
        return score, summary

    def move(variable_index, new_position):
        variable = variables[variable_index]
        old_position = variable['position']
        if new_position == old_position:
            return
        for row_index, sign in adjacent_rows[variable['helix']]:
            old_interval = _arc_index(
                pair_rows[row_index]['events'], old_position)
            new_interval = _arc_index(
                pair_rows[row_index]['events'], new_position)
            if old_interval != new_interval:
                row_differences[row_index][old_interval] -= (
                    sign * variable['value'])
                row_differences[row_index][new_interval] += (
                    sign * variable['value'])
        old_domain = old_position // variable['domain_size']
        new_domain = new_position // variable['domain_size']
        domain_counts[(variable['helix'], old_domain)] -= 1
        domain_counts[(variable['helix'], new_domain)] += 1
        variable['position'] = new_position
        indices = by_helix[variable['helix']]
        helix_cvs[variable['helix']] = _circular_gap_cv(
            [variables[index]['position'] for index in indices], length)

    initial_score, initial = evaluate()
    current_score, current = initial_score, initial
    accepted = 0
    for unused_pass in range(int(maximum_passes)):
        changed = 0
        for variable_index, variable in enumerate(variables):
            old = variable['position']
            old_domain = old // variable['domain_size']
            best = old
            best_score = current_score
            best_summary = current
            for candidate in variable['candidates']:
                if candidate == old:
                    continue
                candidate_domain = candidate // variable['domain_size']
                if (candidate_domain != old_domain and
                        domain_counts[(variable['helix'],
                                       candidate_domain)] >= 3):
                    continue
                move(variable_index, candidate)
                trial_score, trial_summary = evaluate()
                move(variable_index, old)
                if trial_score < best_score:
                    best, best_score = candidate, trial_score
                    best_summary = trial_summary
            if best != old:
                move(variable_index, best)
                current_score, current = best_score, best_summary
                accepted += 1
                changed += 1
        if not changed:
            break

    for row in rows.values():
        row['loop'] = [0] * len(row.get('loop', []))
        row['skip'] = [0] * len(row.get('skip', []))
    by_record = collections.defaultdict(list)
    for variable in variables:
        row = rows[variable['helix']]
        if variable['value'] > 0:
            row['loop'][variable['position']] += variable['value']
        else:
            row['skip'][variable['position']] += variable['value']
        by_record[variable['helix']].append(variable)
    for record in records:
        helix = int(record['helix'])
        selected = by_record.get(helix, [])
        record['insertions'] = sorted(
            variable['position'] for variable in selected
            if variable['value'] > 0)
        record['deletions'] = sorted(
            variable['position'] for variable in selected
            if variable['value'] < 0)
        record['indel_placements'] = [{
            'base': variable['position'], 'value': variable['value'],
            'partition': variable['partition'],
            'partition_count': variable['partition_count'],
            'ideal_base': variable['ideal'],
            'partition_start': variable['partition_start'],
            'partition_end': variable['partition_end'],
            'target_domain': variable['target_domain'],
            'domain': variable['position']//variable['domain_size']}
            for variable in selected]
    audit = {'method': 'pair-aware equal-partition coordinate refinement',
             'accepted_moves': accepted,
             'initial_metrics': initial, 'final_metrics': current,
             'preserved': ['scaffold', 'staples', 'nicks', 'crossovers',
                           'signed indel total on every helix']}
    design.setdefault('curved_metadata', {})[
        'pair_aware_indel_optimization'] = audit
    curvature_indels['pair_aware_indel_optimization'] = audit
    return current


def _linear_interval(events, position):
    for index in range(len(events)-1):
        if events[index] < position < events[index+1]:
            return index
    return None


def _bend_pair_rows(region, helix_data):
    selected = set(int(value) for value in region.get('helices', ()))
    centers = dict((helix, tuple(helix_data[helix]['coord']))
                   for helix in selected if helix in helix_data)
    if len(centers) < 2:
        return []
    distances = [math.dist(centers[left], centers[right])
                 for left in centers for right in centers if right > left]
    if not distances:
        return []
    nearest = min(distances)
    start, end = int(region['start']), int(region['end'])
    result = []
    for left in sorted(centers):
        for right in sorted(centers):
            if right <= left or math.dist(
                    centers[left], centers[right]) > nearest*1.08:
                continue
            sites = []
            for helix, other in ((left, right), (right, left)):
                sites.extend(
                    int(base) for base, partner in
                    helix_data[helix].get('crossovers', ())
                    if int(partner) == other and start < int(base) < end)
            sites = _cluster(sites)
            if not sites:
                continue
            # Synthetic selection boundaries make every new edit part of one
            # axial bin.  Interior entries remain the actual scaffold/staple
            # crossover locations and are never moved by this analysis.
            events = [float(start)-.5] + sites + [float(end)+.5]
            result.append({'pair': (left, right), 'events': events})
    return result


def bend_plan_pair_curvature_data(region, helix_data, plan):
    """Analyze one Add Bending plan over physical adjacent helix pairs."""
    pair_rows = _bend_pair_rows(region, helix_data)
    edits = [edit for edit in plan.get('edits', ())
             if edit.get('operation') != 'remove_existing']
    rows = []
    for pair_row in pair_rows:
        left, right = pair_row['pair']
        differences = [0] * (len(pair_row['events'])-1)
        for edit in edits:
            helix = int(edit['helix'])
            if helix not in (left, right):
                continue
            interval = _linear_interval(
                pair_row['events'], int(edit['idx']))
            if interval is not None:
                differences[interval] += (
                    int(edit['length']) * (1 if helix == right else -1))
        total = sum(differences)
        if not total:
            continue
        inner, outer = left, right
        if total < 0:
            inner, outer = right, left
            differences = [-value for value in differences]
            total = -total
        floor_value, remainder = divmod(total, len(differences))
        rows.append({
            'pair': pair_row['pair'], 'inner': inner, 'outer': outer,
            'events': list(pair_row['events']),
            'differences': differences, 'total_difference': total,
            'mean_difference': total/float(len(differences)),
            'ideal_floor': floor_value,
            'ideal_ceiling': floor_value + (1 if remainder else 0)})
    summary = _pair_rows_summary(
        rows, int(region['end'])-int(region['start'])+1,
        21 if str(plan.get('lattice', '')).lower() == 'honeycomb' else 32,
        origin=float(region['start'])-.5)
    return rows, summary


def _pair_rows_summary(rows, length, bin_width, origin=0.0):
    reverse = outside = zero = above = severity = 0
    for row in rows:
        floor_value, ceiling = row['ideal_floor'], row['ideal_ceiling']
        for difference in row['differences']:
            reverse += difference < 0
            if difference < floor_value or difference > ceiling:
                outside += 1
                severity += min(abs(difference-floor_value),
                                abs(difference-ceiling)) ** 2
            zero += difference == 0 and floor_value >= 1
            above += difference > ceiling
    bin_count = max(1, int(math.ceil(length/float(bin_width))))
    bins = [0.0] * bin_count
    weights = [0.0] * bin_count
    for row in rows:
        mean = float(row.get('mean_difference', 0.0))
        if mean <= 0:
            continue
        for local_base in range(length):
            position = origin + local_base + .314159
            interval = _linear_interval(row['events'], position)
            if interval is None:
                continue
            bin_index = local_base // bin_width
            bins[bin_index] += row['differences'][interval]/mean
            weights[bin_index] += 1.0
    normalized = [value/weight if weight else 0.0
                  for value, weight in zip(bins, weights)]
    positive = [value for value in normalized if value > 0]
    axial_cv = (statistics.pstdev(positive)/statistics.mean(positive)
                if len(positive) > 1 else 0.0)
    return {
        'physical_adjacent_curvature_pairs': len(rows),
        'crossover_intervals': sum(len(row['differences']) for row in rows),
        'reverse_curvature_intervals': reverse,
        'outside_floor_ceiling_intervals': outside,
        'zero_curvature_intervals_where_floor_is_one': zero,
        'above_ceiling_intervals': above,
        'outside_severity': severity,
        'normalized_axial_curvature_bin_width_bp': int(bin_width),
        'normalized_axial_curvature_42bp_bin_cv': axial_cv,
        'normalized_axial_curvature_42bp_bin_minimum': min(positive or [0]),
        'normalized_axial_curvature_42bp_bin_maximum': max(positive or [0]),
        'normalized_axial_curvature_42bp_bins': normalized}


def optimize_bend_plan_curvature(region, helix_data, plan,
                                 maximum_passes=3):
    """Refine Add Bending edit coordinates using the Curved priorities."""
    edits = [edit for edit in plan.get('edits', ())
             if edit.get('operation') != 'remove_existing']
    if not edits:
        rows, summary = bend_plan_pair_curvature_data(
            region, helix_data, plan)
        plan['pair_curvature_rows'] = rows
        plan['pair_curvature_summary'] = summary
        return summary
    pair_rows = _bend_pair_rows(region, helix_data)
    if not pair_rows:
        rows, summary = bend_plan_pair_curvature_data(
            region, helix_data, plan)
        plan['pair_curvature_rows'] = rows
        plan['pair_curvature_summary'] = summary
        return summary

    adjacent = collections.defaultdict(list)
    for row_index, row in enumerate(pair_rows):
        for sign, helix in ((-1, row['pair'][0]), (1, row['pair'][1])):
            adjacent[helix].append((row_index, sign))
    raw = [[0] * (len(row['events'])-1) for row in pair_rows]
    for edit in edits:
        for row_index, sign in adjacent[int(edit['helix'])]:
            interval = _linear_interval(
                pair_rows[row_index]['events'], int(edit['idx']))
            if interval is not None:
                raw[row_index][interval] += sign*int(edit['length'])
    directions = {}
    for row_index, values in enumerate(raw):
        total = sum(values)
        if total:
            directions[row_index] = 1 if total > 0 else -1
    for helix in list(adjacent):
        adjacent[helix] = [
            (row_index, sign*directions[row_index])
            for row_index, sign in adjacent[helix]
            if row_index in directions]
    if not directions:
        rows, summary = bend_plan_pair_curvature_data(
            region, helix_data, plan)
        plan['pair_curvature_rows'] = rows
        plan['pair_curvature_summary'] = summary
        return summary

    row_differences = {}
    for row_index, direction in directions.items():
        values = [direction*value for value in raw[row_index]]
        row_differences[row_index] = values
        total = sum(values)
        floor_value, remainder = divmod(total, len(values))
        pair_rows[row_index]['ideal_floor'] = floor_value
        pair_rows[row_index]['ideal_ceiling'] = (
            floor_value + (1 if remainder else 0))
        pair_rows[row_index]['mean_difference'] = total/float(len(values))

    domain_load = collections.Counter()
    for helix, data in helix_data.items():
        for base, value in data.get('insertions', {}).items():
            domain_size = int(plan.get('domain_size_bp', 8))
            domain_load[(int(helix), int(base)//domain_size)] += int(value)
    for edit in edits:
        domain_size = int(edit.get(
            'domain_size', plan.get('domain_size_bp', 8)))
        domain_load[(int(edit['helix']), int(edit['idx'])//domain_size)] += \
            int(edit['length'])

    by_helix = collections.defaultdict(list)
    planned_occupancy = collections.Counter(
        (int(edit['helix']), int(edit['idx'])) for edit in edits)
    start, end = int(region['start']), int(region['end'])
    for edit_index, edit in enumerate(edits):
        helix = int(edit['helix'])
        domain_size = int(edit.get(
            'domain_size', plan.get('domain_size_bp', 8)))
        partition_start = int(edit.get('partition_start', start))
        partition_end = int(edit.get('partition_end', end))
        forbidden = set(helix_data[helix].get('forbidden', ()))
        existing = set(int(value) for value in
                       helix_data[helix].get('insertions', {}))
        protected = (set(helix_data[helix].get('deletion_protected', ()))
                     if int(edit['length']) < 0 else set())
        candidates = [
            base for base in range(partition_start, partition_end+1)
            if base not in forbidden and base not in existing and
            base not in protected]
        # Candidate signatures preserve the pair interval assignment while
        # keeping only the sites nearest the original equal-partition target.
        signatures = collections.defaultdict(list)
        for base in candidates:
            signature = tuple(
                (row_index,
                 _linear_interval(pair_rows[row_index]['events'], base))
                for row_index, unused_sign in adjacent[helix])
            signatures[signature].append(base)
        shortlisted = {int(edit['idx'])}
        ideal = float(edit.get('ideal_base', edit['idx']))
        for values in signatures.values():
            shortlisted.update(sorted(
                values, key=lambda base: (abs(base-ideal), base))[:2])
        edit['_pair_candidates'] = sorted(shortlisted)
        edit['_pair_ideal'] = ideal
        edit['_pair_domain_size'] = domain_size
        by_helix[helix].append(edit_index)

    helix_cvs = dict((helix, _circular_gap_cv(
        [int(edits[index]['idx']) for index in indices], end-start+1))
        for helix, indices in by_helix.items())
    length = end-start+1
    bin_width = (21 if str(plan.get('lattice', '')).lower() == 'honeycomb'
                 else 32)
    bin_count = max(1, int(math.ceil(length/float(bin_width))))
    interval_bin_weights = {}
    bin_weights = [0.0] * bin_count
    bin_sums = [0.0] * bin_count
    for row_index in sorted(directions):
        source = pair_rows[row_index]
        weights = [[0.0] * bin_count
                   for unused in source['events'][:-1]]
        for local_base in range(length):
            position = float(start)-.5+local_base+.314159
            interval = _linear_interval(source['events'], position)
            if interval is not None:
                weights[interval][local_base//bin_width] += 1.0
        interval_bin_weights[row_index] = weights
        mean = float(source['mean_difference'])
        for interval, difference in enumerate(row_differences[row_index]):
            for bin_index, weight in enumerate(weights[interval]):
                bin_sums[bin_index] += difference/mean*weight
                bin_weights[bin_index] += weight
    distance_sum = sum(
        abs(int(edit['idx'])-float(edit['_pair_ideal'])) for edit in edits)

    def build_rows():
        rows = []
        for row_index in sorted(directions):
            source = pair_rows[row_index]
            values = list(row_differences[row_index])
            rows.append({
                'pair': source['pair'],
                'inner': (source['pair'][0] if directions[row_index] > 0
                          else source['pair'][1]),
                'outer': (source['pair'][1] if directions[row_index] > 0
                          else source['pair'][0]),
                'events': list(source['events']),
                'differences': values,
                'total_difference': sum(values),
                'mean_difference': source['mean_difference'],
                'ideal_floor': source['ideal_floor'],
                'ideal_ceiling': source['ideal_ceiling']})
        return rows

    def metrics_and_score():
        reverse = outside = zero = above = severity = intervals = 0
        for row_index in sorted(directions):
            source = pair_rows[row_index]
            floor_value = int(source['ideal_floor'])
            ceiling = int(source['ideal_ceiling'])
            for difference in row_differences[row_index]:
                intervals += 1
                reverse += difference < 0
                if difference < floor_value or difference > ceiling:
                    outside += 1
                    severity += min(abs(difference-floor_value),
                                    abs(difference-ceiling)) ** 2
                zero += difference == 0 and floor_value >= 1
                above += difference > ceiling
        normalized = [value/weight if weight else 0.0
                      for value, weight in zip(bin_sums, bin_weights)]
        positive = [value for value in normalized if value > 0]
        axial_cv = (statistics.pstdev(positive)/statistics.mean(positive)
                    if len(positive) > 1 else 0.0)
        single = (statistics.mean(helix_cvs.values())
                  if helix_cvs else 0.0)
        distance = distance_sum/float(max(1, len(edits)))
        summary = {
            'physical_adjacent_curvature_pairs': len(directions),
            'crossover_intervals': intervals,
            'reverse_curvature_intervals': reverse,
            'outside_floor_ceiling_intervals': outside,
            'zero_curvature_intervals_where_floor_is_one': zero,
            'above_ceiling_intervals': above,
            'outside_severity': severity,
            'normalized_axial_curvature_bin_width_bp': int(bin_width),
            'normalized_axial_curvature_42bp_bin_cv': axial_cv,
            'normalized_axial_curvature_42bp_bin_minimum':
                min(positive or [0]),
            'normalized_axial_curvature_42bp_bin_maximum':
                max(positive or [0]),
            'normalized_axial_curvature_42bp_bins': normalized,
            'mean_single_helix_spacing_cv': single,
            'mean_distance_from_equal_partition_target_bp': distance}
        score = (
            summary['reverse_curvature_intervals'],
            summary['outside_floor_ceiling_intervals'],
            summary['outside_severity'],
            round(summary['normalized_axial_curvature_42bp_bin_cv'], 9),
            round(single, 9), round(distance, 6))
        return summary, score

    def move(edit_index, position):
        nonlocal distance_sum
        edit = edits[edit_index]
        old = int(edit['idx'])
        if old == position:
            return
        helix, value = int(edit['helix']), int(edit['length'])
        planned_occupancy[(helix, old)] -= 1
        planned_occupancy[(helix, int(position))] += 1
        for row_index, sign in adjacent[helix]:
            old_interval = _linear_interval(
                pair_rows[row_index]['events'], old)
            new_interval = _linear_interval(
                pair_rows[row_index]['events'], position)
            if old_interval != new_interval:
                if old_interval is not None:
                    row_differences[row_index][old_interval] -= sign*value
                if new_interval is not None:
                    row_differences[row_index][new_interval] += sign*value
                mean = float(pair_rows[row_index]['mean_difference'])
                weights = interval_bin_weights[row_index]
                for bin_index in range(bin_count):
                    old_weight = (weights[old_interval][bin_index]
                                  if old_interval is not None else 0.0)
                    new_weight = (weights[new_interval][bin_index]
                                  if new_interval is not None else 0.0)
                    bin_sums[bin_index] += (
                        sign*value/mean*(new_weight-old_weight))
        domain_size = int(edit['_pair_domain_size'])
        domain_load[(helix, old//domain_size)] -= value
        domain_load[(helix, position//domain_size)] += value
        distance_sum += (
            abs(int(position)-float(edit['_pair_ideal'])) -
            abs(old-float(edit['_pair_ideal'])))
        edit['idx'] = int(position)
        indices = by_helix[helix]
        helix_cvs[helix] = _circular_gap_cv(
            [int(edits[index]['idx']) for index in indices], end-start+1)

    initial_summary, current_score = metrics_and_score()
    current_summary = initial_summary
    accepted = 0
    for unused_pass in range(int(maximum_passes)):
        changed = 0
        for edit_index, edit in enumerate(edits):
            old = int(edit['idx'])
            best, best_score = old, current_score
            best_summary = current_summary
            helix, value = int(edit['helix']), int(edit['length'])
            domain_size = int(edit['_pair_domain_size'])
            for candidate in edit['_pair_candidates']:
                if candidate == old:
                    continue
                if planned_occupancy[(helix, int(candidate))]:
                    continue
                old_domain, new_domain = old//domain_size, candidate//domain_size
                if old_domain != new_domain:
                    new_load = domain_load[(helix, new_domain)] + value
                    if abs(new_load) > 3:
                        continue
                move(edit_index, candidate)
                trial_summary, trial_score = metrics_and_score()
                move(edit_index, old)
                if trial_score < best_score:
                    best, best_score = candidate, trial_score
                    best_summary = trial_summary
            if best != old:
                move(edit_index, best)
                current_score = best_score
                current_summary = best_summary
                accepted += 1
                changed += 1
        if not changed:
            break
    for edit in edits:
        edit['domain'] = int(edit['idx'])//int(edit['_pair_domain_size'])
        for private_key in [key for key in edit if key.startswith('_pair_')]:
            edit.pop(private_key, None)
    current_rows = build_rows()
    # Recompute once through the public analyzer as a final consistency
    # check; candidate trials use the equivalent incremental metric above.
    verified_summary = _pair_rows_summary(
        current_rows, length, bin_width, origin=float(start)-.5)
    verified_summary['mean_single_helix_spacing_cv'] = current_summary[
        'mean_single_helix_spacing_cv']
    verified_summary['mean_distance_from_equal_partition_target_bp'] = \
        current_summary['mean_distance_from_equal_partition_target_bp']
    current_summary = verified_summary
    audit = {
        'method': 'pair-aware equal-partition coordinate refinement',
        'accepted_moves': accepted,
        'initial_metrics': initial_summary,
        'final_metrics': current_summary,
        'preserved': ['existing indels', 'crossovers', 'nicks',
                      'signed planned indel total on every helix']}
    plan['pair_curvature_rows'] = current_rows
    plan['pair_curvature_summary'] = current_summary
    plan['pair_aware_indel_optimization'] = audit
    return current_summary


def curved_pair_curvature_data(design):
    rows = {int(row['num']): row for row in design.get('vstrands', [])}
    length = int(design.get('num_bases') or max(
        [len(row.get('loop', [])) for row in rows.values()] or [0]))
    variables = []
    for helix, row in rows.items():
        for base, value in enumerate(row.get('loop', [])):
            variables.extend({'helix': helix, 'position': base, 'value': 1}
                             for unused in range(max(0, int(value))))
        for base, value in enumerate(row.get('skip', [])):
            variables.extend({'helix': helix, 'position': base, 'value': -1}
                             for unused in range(max(0, -int(value))))
    return _pair_analysis(_physical_pair_rows(design, rows), variables,
                          max(1, length))


def _frame_windows(design):
    """Return the vertex-local bend windows stored by Frame Design."""
    metadata = dict(design.get('curved_metadata', {}) or {})
    plan = dict(metadata.get('frame_plan', {}) or {})
    centres = list(plan.get('vertex_native_centres', ()))
    lengths = list(plan.get('bend_length_bp', ()))
    return [
        {'vertex': index, 'start': float(centre)-float(length)/2.0,
         'end': float(centre)+float(length)/2.0,
         'length': float(length)}
        for index, (centre, length) in enumerate(zip(centres, lengths))]


def _frame_variables(design):
    variables = []
    for record in (design.get('curvature_indels', {}) or {}).get(
            'rings', ()):
        helix = int(record['helix'])
        domain_size = int(record.get('domain_size_bp') or
                          (design.get('curvature_indels', {}) or {}).get(
                              'domain_size_bp', 7))
        for vertex_record in record.get('frame_vertices', ()):
            vertex = int(vertex_record['vertex'])
            value = 1 if int(vertex_record.get('signed_indel', 0)) > 0 \
                else -1
            sites = list(map(int, vertex_record.get('sites', ())))
            count = len(sites)
            start = float(vertex_record['window_start'])
            end = float(vertex_record['window_end'])
            for order, site in enumerate(sites):
                left = start + (end-start)*order/float(max(1, count))
                right = start + (end-start)*(order+1)/float(max(1, count))
                variables.append({
                    'helix': helix, 'vertex': vertex, 'value': value,
                    'position': site, 'order': order,
                    'partition_count': count, 'partition_start': left,
                    'partition_end': right,
                    'ideal': 0.5*(left+right),
                    'domain_size': domain_size})
    return variables


def _frame_interval_rows(design, rows, variables):
    """Build only the physical crossover intervals inside Frame corners.

    Straight-edge intervals are deliberately absent: zero curvature there is
    the design target and must never be reported as a below-range defect.
    """
    windows = _frame_windows(design)
    length = int(design.get('num_bases') or max(
        [len(row.get('loop', [])) for row in rows.values()] or [0]))
    physical = _physical_pair_rows(design, rows)
    by_helix_vertex = collections.defaultdict(list)
    for variable in variables:
        by_helix_vertex[(int(variable['helix']),
                         int(variable['vertex']))].append(variable)
    result = []
    for source in physical:
        for window in windows:
            vertex = int(window['vertex'])
            # Select every crossover interval touched by a base in the bend
            # window.  This includes the two boundary intervals and ensures
            # that every permitted Frame indel belongs to a reported cell.
            global_indices = []
            first = max(0, int(math.floor(window['start'])))
            last = min(length-1, int(math.ceil(window['end'])))
            for base in range(first, last+1):
                index = _arc_index(source['events'], base+.314159)
                if index is not None and index not in global_indices:
                    global_indices.append(index)
            global_indices.sort(key=lambda index: (
                (source['events'][index]-window['start']) % length))
            if not global_indices:
                continue
            index_map = dict((global_index, local_index)
                             for local_index, global_index in
                             enumerate(global_indices))
            differences = [0] * len(global_indices)
            for helix, sign in ((source['inner'], -1),
                                (source['outer'], 1)):
                for variable in by_helix_vertex.get((helix, vertex), ()):
                    global_index = _arc_index(
                        source['events'], float(variable['position']))
                    if global_index in index_map:
                        differences[index_map[global_index]] += \
                            sign*int(variable['value'])
            total = sum(differences)
            if not total:
                continue
            inner, outer = source['inner'], source['outer']
            if total < 0:
                inner, outer = outer, inner
                differences = [-value for value in differences]
                total = -total
            floor_value, remainder = divmod(total, len(differences))
            starts = [source['events'][index] for index in global_indices]
            ends = [source['events'][(index+1) % len(source['events'])]
                    for index in global_indices]
            result.append({
                'pair': source['pair'], 'inner': inner, 'outer': outer,
                'vertex': vertex, 'label': 'H%d–H%d · V%d' % (
                    source['pair'][0], source['pair'][1], vertex+1),
                'events': starts, 'interval_starts': starts,
                'interval_ends': ends,
                'global_intervals': global_indices,
                'source_events': list(source['events']),
                'differences': differences, 'total_difference': total,
                'mean_difference': total/float(len(differences)),
                'ideal_floor': floor_value,
                'ideal_ceiling': floor_value+(1 if remainder else 0)})
    return result, max(1, length)


def _frame_rows_summary(rows, bin_width):
    reverse = outside = zero = above = severity = 0
    by_vertex = collections.defaultdict(list)
    for row in rows:
        floor_value, ceiling = row['ideal_floor'], row['ideal_ceiling']
        by_vertex[int(row.get('vertex', 0))].append(row)
        for difference in row['differences']:
            reverse += difference < 0
            if difference < floor_value or difference > ceiling:
                outside += 1
                severity += min(abs(difference-floor_value),
                                abs(difference-ceiling)) ** 2
            zero += difference == 0 and floor_value >= 1
            above += difference > ceiling
    # Compare like-for-like axial fractions across the bend windows.  This
    # excludes straight edges and makes vertices of different length
    # comparable without pretending that they are one continuous bend.
    profiles = []
    for vertex in sorted(by_vertex):
        vertex_rows = by_vertex[vertex]
        cell_count = max([len(row['differences']) for row in vertex_rows] or
                         [0])
        for cell in range(cell_count):
            values = []
            for row in vertex_rows:
                if cell >= len(row['differences']) or \
                        float(row['mean_difference']) <= 0:
                    continue
                values.append(row['differences'][cell] /
                              float(row['mean_difference']))
            if values:
                profiles.append(statistics.mean(values))
    positive = [value for value in profiles if value > 0]
    axial_cv = (statistics.pstdev(positive)/statistics.mean(positive)
                if len(positive) > 1 else 0.0)
    return {
        'scope': 'frame-bend-windows-only',
        'physical_adjacent_curvature_pairs': len(set(
            tuple(row['pair']) for row in rows)),
        'reported_pair_vertex_rows': len(rows),
        'crossover_intervals': sum(len(row['differences']) for row in rows),
        'reverse_curvature_intervals': reverse,
        'outside_floor_ceiling_intervals': outside,
        'zero_curvature_intervals_where_floor_is_one': zero,
        'above_ceiling_intervals': above,
        'outside_severity': severity,
        'normalized_axial_curvature_bin_width_bp': int(bin_width),
        'normalized_axial_curvature_42bp_bin_cv': axial_cv,
        'normalized_axial_curvature_42bp_bin_minimum': min(positive or [0]),
        'normalized_axial_curvature_42bp_bin_maximum': max(positive or [0]),
        'normalized_axial_curvature_42bp_bins': positive}


def frame_pair_curvature_data(design):
    rows = {int(row['num']): row for row in design.get('vstrands', [])}
    variables = _frame_variables(design)
    pair_rows, unused_length = _frame_interval_rows(
        design, rows, variables)
    lattice = str(design.get('lattice') or
                  (design.get('curved_metadata', {}) or {}).get(
                      'lattice', 'honeycomb')).lower()
    return pair_rows, _frame_rows_summary(
        pair_rows, 21 if lattice == 'honeycomb' else 32)


def generated_single_helix_distribution_data(design, frame_only=False):
    """Single-helix equal-spacing diagnostics for generated indels."""
    variables = _frame_variables(design) if frame_only else []
    if not frame_only:
        for record in (design.get('curvature_indels', {}) or {}).get(
                'rings', ()):
            helix = int(record['helix'])
            for base in record.get('insertions', ()):
                variables.append({'helix': helix, 'vertex': -1,
                                  'position': int(base), 'value': 1})
            for base in record.get('deletions', ()):
                variables.append({'helix': helix, 'vertex': -1,
                                  'position': int(base), 'value': -1})
    grouped = collections.defaultdict(list)
    for variable in variables:
        grouped[(int(variable['helix']),
                 int(variable.get('vertex', -1)))].append(
                     int(variable['position']))
    # A Frame report is a complete helix-by-vertex audit of bend sections,
    # including neutral helices with zero indels.  Straight sections never
    # become groups and therefore cannot enter this report or its statistics.
    if frame_only:
        vertices = [int(item['vertex']) for item in _frame_windows(design)]
        for record in (design.get('curvature_indels', {}) or {}).get(
                'rings', ()):
            helix = int(record['helix'])
            for vertex in vertices:
                grouped.setdefault((helix, vertex), [])
    length = int(design.get('num_bases') or max(
        [len(row.get('loop', [])) for row in design.get('vstrands', [])]
        or [1]))
    windows = dict((int(item['vertex']), item)
                   for item in _frame_windows(design))
    result = []
    for (helix, vertex), positions in sorted(grouped.items()):
        ordered = sorted(positions)
        gaps = [right-left for left, right in zip(ordered, ordered[1:])]
        if vertex < 0 and len(ordered) > 1:
            gaps.append(length-ordered[-1]+ordered[0])
        mean = statistics.mean(gaps) if gaps else 0.0
        cv = (statistics.pstdev(gaps)/mean
              if len(gaps) > 1 and mean else 0.0)
        result.append({
            'helix': helix, 'vertex': vertex, 'count': len(ordered),
            'positions': ordered, 'mean_spacing_bp': mean,
            'spacing_cv': cv, 'minimum_spacing_bp': min(gaps or [0]),
            'maximum_spacing_bp': max(gaps or [0]),
            'window_start': windows.get(vertex, {}).get('start', 0.0),
            'window_end': windows.get(vertex, {}).get('end', float(length)),
            'scope': ('bend-window' if vertex >= 0 else 'whole-loop')})
    return result


def write_generated_single_helix_distribution_csv(rows, path):
    """Write Frame/Curved generated-indel spacing diagnostics."""
    directory = os.path.dirname(os.path.abspath(path))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, 'w', newline='', encoding='utf-8') as output:
        writer = csv.writer(output)
        writer.writerow([
            'helix', 'vertex', 'scope', 'indel_count', 'positions',
            'window_start', 'window_end', 'mean_spacing_bp', 'spacing_cv',
            'minimum_spacing_bp', 'maximum_spacing_bp'])
        for row in rows:
            vertex = int(row.get('vertex', -1))
            writer.writerow([
                row['helix'], vertex+1 if vertex >= 0 else '',
                row.get('scope', ''), row['count'],
                ';'.join(str(value) for value in row.get('positions', ())),
                row.get('window_start', ''), row.get('window_end', ''),
                row.get('mean_spacing_bp', 0.0),
                row.get('spacing_cv', 0.0),
                row.get('minimum_spacing_bp', 0),
                row.get('maximum_spacing_bp', 0)])


def optimize_frame_pair_curvature(design, maximum_passes=20):
    """Refine Frame corner indels without touching AutoCS topology.

    Every move remains inside the same vertex bend window.  The strict score
    order is reverse curvature, outside theoretical floor/ceiling, severity,
    axial fluctuation, single-helix spacing CV, then distance from the
    equal-partition target.  Signed indel totals and the user-selected
    +/-1, +/-2 or +/-3 per-domain limit are invariant.  Iteration stops early
    at the first complete pass with no
    improvement; the larger bound only prevents the former three-pass cap
    from abandoning still-improving Frame solutions.
    """
    rows = {int(row['num']): row for row in design.get('vstrands', [])}
    variables = _frame_variables(design)
    state_rows, unused_length = _frame_interval_rows(
        design, rows, variables)
    if not variables or not state_rows:
        return _frame_rows_summary(state_rows, 42)
    lattice = str(design.get('lattice') or
                  (design.get('curved_metadata', {}) or {}).get(
                      'lattice', 'honeycomb')).lower()
    bin_width = 21 if lattice == 'honeycomb' else 32
    maximum_per_domain = max(1, min(3, int(
        (design.get('curvature_indels', {}) or {}).get(
            'maximum_indel_per_domain_allowed', 3))))
    deletion_protected = _short_staple_deletion_protection(rows)
    safe = {}
    for helix, row in rows.items():
        size = min(len(row.get('scaf', ())), len(row.get('stap', ())))
        safe[helix] = [
            base for base in range(1, max(1, size-1))
            if _longitudinal(row['scaf'][base], helix, base) and
            _longitudinal(row['stap'][base], helix, base)]

    adjacent = collections.defaultdict(list)
    global_to_local = []
    for row_index, row in enumerate(state_rows):
        mapping = dict((global_index, local_index)
                       for local_index, global_index in enumerate(
                           row['global_intervals']))
        global_to_local.append(mapping)
        adjacent[(int(row['inner']), int(row['vertex']))].append(
            (row_index, -1))
        adjacent[(int(row['outer']), int(row['vertex']))].append(
            (row_index, 1))

    domain_counts = collections.Counter()
    site_counts = collections.Counter()
    by_group = collections.defaultdict(list)
    for variable_index, variable in enumerate(variables):
        helix = int(variable['helix'])
        position = int(variable['position'])
        domain = position//int(variable['domain_size'])
        domain_counts[(helix, domain)] += 1
        site_counts[(helix, position, int(variable['value']))] += 1
        by_group[(helix, int(variable['vertex']))].append(variable_index)

    # Candidate signatures represent the physical crossover cells affected
    # by a coordinate.  Keep the closest safe coordinate for each signature;
    # this permits pair correction without sacrificing equal partition more
    # than necessary.
    for variable in variables:
        helix = int(variable['helix'])
        vertex = int(variable['vertex'])
        start = min(float(variable['partition_start']),
                    float(variable['partition_end']))
        end = max(float(variable['partition_start']),
                  float(variable['partition_end']))
        # First honour the assigned equal partition.  If it cannot change the
        # pair-cell signature, allow the full bend window as a fallback; the
        # equal-spacing CV and ideal-distance terms remain active penalties.
        window_record = next(
            item for item in _frame_windows(design)
            if int(item['vertex']) == vertex)
        partition_candidates = [
            base for base in safe.get(helix, ()) if start <= base < end and
            (int(variable['value']) > 0 or
             base not in deletion_protected[helix])]
        window_candidates = [
            base for base in safe.get(helix, ())
            if float(window_record['start']) <= base <=
            float(window_record['end']) and
            (int(variable['value']) > 0 or
             base not in deletion_protected[helix])]

        def signature(base):
            values = []
            for row_index, unused_sign in adjacent.get((helix, vertex), ()):
                source_events = state_rows[row_index]['source_events']
                global_index = _arc_index(source_events, base)
                values.append((row_index,
                               global_to_local[row_index].get(global_index)))
            return tuple(values)

        current_signature = signature(int(variable['position']))
        partition_signatures = set(signature(base)
                                   for base in partition_candidates)
        pool = (partition_candidates if len(partition_signatures) > 1
                else window_candidates)
        signatures = collections.defaultdict(list)
        for base in pool:
            signatures[signature(base)].append(base)
        candidates = {int(variable['position'])}
        for values in signatures.values():
            candidates.update(sorted(values, key=lambda base: (
                abs(base-float(variable['ideal'])), base))[:2])
        # Always retain the current equal-partition cell even when a safe-site
        # peculiarity makes its signature unique.
        if current_signature not in signatures:
            candidates.add(int(variable['position']))
        variable['_frame_candidates'] = sorted(candidates)

    def group_cv(key):
        positions = sorted(int(variables[index]['position'])
                           for index in by_group[key])
        gaps = [right-left for left, right in zip(positions, positions[1:])]
        mean = statistics.mean(gaps) if gaps else 0.0
        return (statistics.pstdev(gaps)/mean
                if len(gaps) > 1 and mean else 0.0)

    group_cvs = dict((key, group_cv(key)) for key in by_group)

    def summary_and_score():
        summary = _frame_rows_summary(state_rows, bin_width)
        single_cv = statistics.mean(group_cvs.values()) if group_cvs else 0
        distance = statistics.mean(
            abs(float(variable['position'])-float(variable['ideal']))
            for variable in variables)
        summary['mean_single_helix_spacing_cv'] = single_cv
        summary['mean_distance_from_equal_partition_target_bp'] = distance
        score = (
            int(summary['reverse_curvature_intervals']),
            int(summary['outside_floor_ceiling_intervals']),
            int(summary['outside_severity']),
            round(float(summary['normalized_axial_curvature_42bp_bin_cv']),
                  9),
            round(float(single_cv), 9), round(float(distance), 6))
        return summary, score

    def move(variable_index, new_position):
        variable = variables[variable_index]
        old_position = int(variable['position'])
        new_position = int(new_position)
        if old_position == new_position:
            return
        helix = int(variable['helix'])
        vertex = int(variable['vertex'])
        value = int(variable['value'])
        for row_index, sign in adjacent.get((helix, vertex), ()):
            source_events = state_rows[row_index]['source_events']
            old_global = _arc_index(source_events, old_position)
            new_global = _arc_index(source_events, new_position)
            old_local = global_to_local[row_index].get(old_global)
            new_local = global_to_local[row_index].get(new_global)
            if old_local == new_local:
                continue
            if old_local is not None:
                state_rows[row_index]['differences'][old_local] -= sign*value
            if new_local is not None:
                state_rows[row_index]['differences'][new_local] += sign*value
        domain_size = int(variable['domain_size'])
        domain_counts[(helix, old_position//domain_size)] -= 1
        domain_counts[(helix, new_position//domain_size)] += 1
        site_counts[(helix, old_position, value)] -= 1
        site_counts[(helix, new_position, value)] += 1
        variable['position'] = new_position
        group_cvs[(helix, vertex)] = group_cv((helix, vertex))

    initial_summary, current_score = summary_and_score()
    current_summary = initial_summary
    accepted = 0
    for unused_pass in range(max(1, int(maximum_passes))):
        changed = 0
        for variable_index, variable in enumerate(variables):
            old = int(variable['position'])
            best = old
            best_summary = current_summary
            best_score = current_score
            helix = int(variable['helix'])
            value = int(variable['value'])
            domain_size = int(variable['domain_size'])
            for candidate in variable.get('_frame_candidates', (old,)):
                candidate = int(candidate)
                if candidate == old:
                    continue
                if value < 0 and site_counts[(helix, candidate, value)] > 0:
                    continue
                old_domain, new_domain = (old//domain_size,
                                          candidate//domain_size)
                if (old_domain != new_domain and
                        domain_counts[(helix, new_domain)] >=
                        maximum_per_domain):
                    continue
                move(variable_index, candidate)
                trial_summary, trial_score = summary_and_score()
                move(variable_index, old)
                if trial_score < best_score:
                    best, best_score = candidate, trial_score
                    best_summary = trial_summary
            if best != old:
                move(variable_index, best)
                current_score = best_score
                current_summary = best_summary
                accepted += 1
                changed += 1
        if not changed:
            break

    # Rewrite only the Frame-managed indel arrays and metadata.  Scaffold and
    # staple links are deliberately untouched.
    curvature = dict(design.get('curvature_indels', {}) or {})
    records = list(curvature.get('rings', ()))
    for record in records:
        helix = int(record['helix'])
        row = rows[helix]
        for base in record.get('insertions', ()):
            row['loop'][int(base)] = 0
        for base in record.get('deletions', ()):
            row['skip'][int(base)] = 0
    positioned = collections.defaultdict(list)
    for variable in variables:
        positioned[(int(variable['helix']),
                    int(variable['vertex']))].append(variable)
    for record in records:
        helix = int(record['helix'])
        row = rows[helix]
        inserted, deleted = [], []
        domain_load = collections.Counter()
        for vertex_record in record.get('frame_vertices', ()):
            vertex = int(vertex_record['vertex'])
            values = sorted(positioned[(helix, vertex)],
                            key=lambda item: int(item['order']))
            sites = [int(item['position']) for item in values]
            vertex_record['sites'] = sites
            vertex_record['placements'] = [{
                'base': int(item['position']),
                'ideal_base': float(item['ideal']),
                'partition_start': float(item['partition_start']),
                'partition_end': float(item['partition_end'])}
                for item in values]
            for item in values:
                base, value = int(item['position']), int(item['value'])
                domain_load[base//int(item['domain_size'])] += 1
                if value > 0:
                    row['loop'][base] += 1
                    inserted.append(base)
                else:
                    row['skip'][base] = -1
                    deleted.append(base)
        record['insertions'] = sorted(inserted)
        record['deletions'] = sorted(deleted)
        record['domain_indel_quotas'] = [
            domain_load[index] for index in range(int(math.ceil(
                int(record.get('nominal_bases', len(row['loop']))) /
                float(int(record.get('domain_size_bp', 7))))))]
        record['maximum_indel_in_one_domain'] = max(
            domain_load.values() or [0])
    final_rows, unused_length = _frame_interval_rows(
        design, rows, _frame_variables(design))
    final_summary = _frame_rows_summary(final_rows, bin_width)
    final_summary['mean_single_helix_spacing_cv'] = current_summary.get(
        'mean_single_helix_spacing_cv', 0.0)
    curvature['pair_aware_indel_optimization'] = {
        'method': 'Frame vertex-local pair-aware coordinate refinement',
        'accepted_moves': accepted, 'initial_metrics': initial_summary,
        'final_metrics': final_summary,
        'preserved': ['AutoCS topology', 'vertex bend windows',
                      'signed per-helix/per-vertex indel totals',
                      'maximum %d indels per native domain' %
                      maximum_per_domain]}
    design['curvature_indels'] = curvature
    return final_summary


def write_pair_curvature_csv(rows, path):
    directory = os.path.dirname(os.path.abspath(path))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, 'w', newline='', encoding='utf-8') as output:
        writer = csv.writer(output)
        writer.writerow(['pair', 'vertex', 'inner_helix', 'outer_helix',
                         'interval_index', 'start_base', 'end_base',
                         'outer_minus_inner', 'ideal_floor', 'ideal_ceiling',
                         'reverse', 'outside_floor_ceiling'])
        for row in rows:
            for index, difference in enumerate(row['differences']):
                starts = row.get('interval_starts', row['events'])
                ends = row.get('interval_ends')
                start = starts[index]
                end = (ends[index] if ends is not None else
                       row['events'][(index+1) % len(row['events'])])
                writer.writerow([
                    '%d-%d' % row['pair'],
                    (int(row['vertex'])+1 if 'vertex' in row else ''),
                    row['inner'], row['outer'], index,
                    start, end, difference, row['ideal_floor'],
                    row['ideal_ceiling'], difference < 0,
                    difference < row['ideal_floor'] or
                    difference > row['ideal_ceiling']])


def write_pair_curvature_svg(rows, summary, path):
    """Write the final adjacent-helix-pair curvature map as editable SVG."""
    width = 1400
    margin = 42
    label_width = 100
    row_height = 30
    header_height = 112
    height = header_height + max(1, len(rows)) * row_height + 42
    plot_width = width - margin * 2 - label_width
    colors = {
        'valid': '#2A9D8F', 'reverse': '#D73027',
        'below': '#F4A261', 'above': '#7B2CBF'}

    def esc(value):
        return html.escape(str(value), quote=True)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
        'viewBox="0 0 %d %d">' % (width, height, width, height),
        '<rect width="100%" height="100%" fill="#FAFBFD"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#24272C}'
        '.title{font-size:20px;font-weight:700}.summary{font-size:13px}'
        '.pair{font-size:12px;font-weight:700}.cell{font-size:9px;'
        'font-weight:700;fill:#fff;text-anchor:middle}</style>',
        '<text class="title" x="%d" y="32">Adjacent-helix-pair curvature '
        'distribution</text>' % margin,
        '<text class="summary" x="%d" y="58">Pairs: %s · crossover '
        'intervals: %s · reverse: %s · outside theoretical interval: %s · '
        '%d-bp axial CV: %.2f%%</text>' % (
            margin,
            esc(summary.get('physical_adjacent_curvature_pairs', 0)),
            esc(summary.get('crossover_intervals', 0)),
            esc(summary.get('reverse_curvature_intervals', 0)),
            esc(summary.get('outside_floor_ceiling_intervals', 0)),
            int(summary.get(
                'normalized_axial_curvature_bin_width_bp', 42)),
            100.0 * float(summary.get(
                'normalized_axial_curvature_42bp_bin_cv', 0.0))),
        '<text class="summary" x="%d" y="82">Green: within theoretical '
        'floor/ceiling · red: reverse · orange: below · purple: above</text>'
        % margin]
    for row_index, row in enumerate(rows):
        y = header_height + row_index * row_height
        pair_label = row.get('label', 'H%d–H%d' % tuple(row['pair']))
        pair_label += ' · target %d–%d' % (
            int(row.get('ideal_floor', 0)),
            int(row.get('ideal_ceiling', row.get('ideal_floor', 0))))
        lines.append('<text class="pair" x="%d" y="%d">%s</text>' % (
            margin, y + 19, esc(pair_label)))
        values = list(row.get('differences', []))
        cell_width = plot_width / float(max(1, len(values)))
        floor_value = int(row.get('ideal_floor', 0))
        ceiling = int(row.get('ideal_ceiling', floor_value))
        for index, value in enumerate(values):
            state = ('reverse' if value < 0 else
                     'below' if value < floor_value else
                     'above' if value > ceiling else 'valid')
            x = margin + label_width + index * cell_width
            lines.append(
                '<rect x="%.3f" y="%d" width="%.3f" height="23" '
                'rx="2" fill="%s" stroke="#FAFBFD" stroke-width="0.8"/>'
                % (x, y, max(1.0, cell_width), colors[state]))
            if cell_width >= 12:
                lines.append(
                    '<text class="cell" x="%.3f" y="%d">%s</text>' % (
                        x + cell_width / 2.0, y + 16, esc(value)))
    lines.append('</svg>')
    directory = os.path.dirname(os.path.abspath(path))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, 'w', encoding='utf-8', newline='\n') as output:
        output.write('\n'.join(lines))


def single_helix_distribution_data(plans):
    result = []
    for plan_index, plan in enumerate(plans or []):
        if plan.get('kind') not in ('add_twist', 'remove_twist', 'bend'):
            continue
        by_helix = collections.defaultdict(list)
        for edit in plan.get('edits', []):
            if edit.get('operation') == 'remove_existing':
                continue
            by_helix[int(edit['helix'])].append(int(edit['idx']))
        for helix, positions in sorted(by_helix.items()):
            ordered = sorted(positions)
            gaps = [right-left for left, right in zip(ordered, ordered[1:])]
            mean = statistics.mean(gaps) if gaps else 0.0
            cv = (statistics.pstdev(gaps)/mean
                  if len(gaps) > 1 and mean else 0.0)
            result.append({
                'plan': plan_index + 1, 'kind': plan.get('kind'),
                'helix': helix, 'count': len(ordered),
                'positions': ordered, 'mean_spacing_bp': mean,
                'spacing_cv': cv,
                'minimum_spacing_bp': min(gaps or [0]),
                'maximum_spacing_bp': max(gaps or [0])})
    return result


def write_single_helix_distribution_csv(rows, path):
    directory = os.path.dirname(os.path.abspath(path))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, 'w', newline='', encoding='utf-8') as output:
        writer = csv.writer(output)
        writer.writerow([
            'plan', 'kind', 'helix', 'indel_count', 'positions',
            'mean_spacing_bp', 'spacing_cv', 'minimum_spacing_bp',
            'maximum_spacing_bp'])
        for row in rows:
            writer.writerow([
                row['plan'], row['kind'], row['helix'], row['count'],
                ';'.join(str(value) for value in row['positions']),
                row['mean_spacing_bp'], row['spacing_cv'],
                row['minimum_spacing_bp'], row['maximum_spacing_bp']])

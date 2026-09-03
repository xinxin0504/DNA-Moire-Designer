#!/usr/bin/env python3
"""Isolated caDNAno model worker used by the Moire Designer GUI."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
import importlib.util
import itertools
import json
import math
import os
import sys
from pathlib import Path


HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parents[2]
CADNANO_CANDIDATES = (
    PROJECT_ROOT / "work" / "cadnano2-modified",       # source checkout
    HERE.parents[1] / "designer_vendor" / "cadnano2", # packaged engine
)
CADNANO_ROOT = next(
    (candidate for candidate in CADNANO_CANDIDATES
     if (candidate / "__init__.py").is_file()), None)
if CADNANO_ROOT is None:
    raise RuntimeError("找不到用于结构生成的cadnano模型源码。")
sys.path.insert(0, str(CADNANO_ROOT))
package_spec = importlib.util.spec_from_file_location(
    "cadnano2", CADNANO_ROOT / "__init__.py",
    submodule_search_locations=[str(CADNANO_ROOT)])
package = importlib.util.module_from_spec(package_spec)
sys.modules["cadnano2"] = package
package_spec.loader.exec_module(package)

import cadnano2.util as _cadnano_util

# The GUI virtualenv contains PyQt6, but the structure worker must remain
# headless and architecture-independent.  Force cadnano's Dummy binding before
# importing any model module; otherwise Qt can be loaded merely because it is
# installed in the parent GUI environment.
_cadnano_util.qtFrameworkList = ["Dummy"]

import cadnano2.cadnano as cadnano


app = cadnano.initAppWithoutGui([])
app.prefs.squareRows = 40
app.prefs.squareCols = 40
app.prefs.squareSteps = 12


class DummySignal:
    def emit(self, *unused_args):
        pass


app.documentWasCreatedSignal = DummySignal()

import dummyqt.QtCore as _dummy_qt_core
import dummyqt.QtGui as _dummy_qt_gui


class _DummyWeight:
    Bold = None


_dummy_qt_gui.QFont.Weight = _DummyWeight
_dummy_qt_gui.QColor.name = lambda self: "#000000"
# Recent cadnano crossover code queries the undo stack even in headless model
# workers.  DummyQt predates those read-only methods; provide them locally so
# structure generation never falls through to a GUI Qt dependency.
if not hasattr(_dummy_qt_gui.QUndoStack, "canUndo"):
    _dummy_qt_gui.QUndoStack.canUndo = lambda self: False
if not hasattr(_dummy_qt_gui.QUndoStack, "canRedo"):
    _dummy_qt_gui.QUndoStack.canRedo = lambda self: False


class _DummyQObject:
    def __init__(self, *unused_args, **unused_kwargs):
        self._dummy_parent = unused_kwargs.get("parent")

    def setParent(self, parent):
        self._dummy_parent = parent

    def parent(self):
        return self._dummy_parent

    def deleteLater(self):
        pass


_dummy_qt_core.QObject = _DummyQObject


class _DummyBoundSignal:
    def __init__(self):
        self.targets = []

    def connect(self, target):
        if target not in self.targets:
            self.targets.append(target)

    def disconnect(self, target):
        if target in self.targets:
            self.targets.remove(target)

    def emit(self, *args):
        for target in list(self.targets):
            target(*args)


class _DummySignalDescriptor:
    def __init__(self, *unused_args):
        self.instances = {}

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        return self.instances.setdefault(id(instance), _DummyBoundSignal())


_dummy_qt_core.pyqtSignal = _DummySignalDescriptor

from cadnano2.model.document import Document

sys.path.insert(0, str(HERE.parents[1]))
from moire_design_core.calculations import (
    MAX_SEED_DELETIONS_PER_DOMAIN,
    SQUARE_DOMAIN_BP,
    minimum_seed_deletion_per_helix,
)
from moire_design_core.structure import (
    CAPTURE_COLUMN_COLORS,
    CAPTURE_DIRECT_POSITIONS,
    CAPTURE_EXTENSION_NT,
    CAPTURE_FACE_DEFINITIONS,
    CAPTURE_PAIR_COLORS,
    CAPTURE_REFERENCE_COLUMN_BY_COLOR,
    CAPTURE_EXPORT_PHASE_CYCLE,
    CAPTURE_PHASE_CYCLE,
    CAPTURE_PHASE_MAPPINGS,
    SCAFFOLD_CAPACITY_P8064,
    SCAFFOLD_CAPACITY_ORTHOGONAL,
    SCAFFOLD_MAX_COUNT,
    SEED_CAPTURE_REFERENCE,
    SQUARE_SCAF_HIGH,
    SQUARE_SCAF_LOW,
    SST_ARRAY_LENGTH,
    build_capture_ready_sst_payload,
    build_shifted_sst_payload,
    capture_column_color,
    capture_extension_nt,
    capture_pair_index,
    capture_export_site_assignments,
    fixed_seed_overlap_layout,
    scaffold_capacity_plan,
    capture_site_assignments,
    payload_to_internal_numbering,
    payload_to_sst_first_numbering,
    structure_layout,
    _copy_periodic_segment,
    _periodic_source_index,
)


def _reference():
    # This is deliberately the frozen two-layer Square routing baseline.
    # Kagome's historical three-layer JSON is consulted only by the
    # lattice-specific capture projection and is never returned here.
    return json.loads(SEED_CAPTURE_REFERENCE.read_text(encoding="utf-8"))


def _seed_canvas_shift(layout):
    """Return the one allowed geometric change to the frozen Seed.

    ``centered_square_sst_geometry`` can move the complete design by a whole
    32-bp routing period when a long left SST layer would otherwise cross
    base 0.  The Seed topology is still immutable: every base, partner base,
    nick, crossover, colour and indel coordinate is translated by the same
    amount.  No local Seed crop/growth/repair is performed here.
    """
    shift = int(layout.get("coordinate_shift_bp", 0) or 0)
    if shift < 0 or shift % 32:
        raise RuntimeError(
            "固定Seed只能按非负32-bp完整周期整体平移；当前为%d bp。" %
            shift)
    return shift


def _seed_capture_support_ranges(layout):
    """Return the actual Seed Z1/Z3 contact intervals for this spacing.

    The frozen Seed routing is translated as a whole, but its Z2 interval is
    locked to the requested layer spacing.  Capture support is consequently
    the first and third Seed partitions, not the nominal 128-bp template
    ranges and not the complete SST duplex ranges.
    """
    geometry = layout.get(
        "square_centered_geometry", layout.get("centered_geometry", {}))
    partitions = (layout.get("seed_partition_ranges") or
                  geometry.get("seed_partition_ranges"))
    if partitions and len(partitions) == 3:
        return [copy.deepcopy(partitions[0]), copy.deepcopy(partitions[2])]
    overlap_ranges = (layout.get("overlap_ranges") or
                      geometry.get("optimized_overlap_ranges"))
    if overlap_ranges and len(overlap_ranges) == 2:
        return copy.deepcopy(overlap_ranges)
    return copy.deepcopy(layout.get(
        "capture_support_ranges", layout.get("seed_layer_ranges", [
            [48 + _seed_canvas_shift(layout),
             175 + _seed_canvas_shift(layout)],
            [208 + _seed_canvas_shift(layout),
             335 + _seed_canvas_shift(layout)],
        ])))


def _translated_frozen_reference(layout):
    """Return the frozen 2L reference translated on the current canvas."""
    source = copy.deepcopy(_reference())
    shift = _seed_canvas_shift(layout)
    source_length = max(
        int(source.get("num_bases", 0) or 0),
        max(len(row.get("scaf", ())) for row in source["vstrands"]))
    target_length = max(
        source_length + shift,
        int(layout.get("array_length", source_length + shift) or 0))

    for row in source["vstrands"]:
        for field in ("scaf", "stap"):
            translated = [[-1, -1, -1, -1]
                          for unused in range(target_length)]
            for source_index, record in enumerate(row.get(field, [])):
                target_index = source_index + shift
                if target_index >= target_length:
                    break
                translated_record = list(record)
                for offset in (0, 2):
                    if int(translated_record[offset]) >= 0:
                        translated_record[offset + 1] = int(
                            translated_record[offset + 1]) + shift
                translated[target_index] = translated_record
            row[field] = translated
        for field in ("loop", "skip"):
            translated = [0] * target_length
            for source_index, value in enumerate(row.get(field, [])):
                target_index = source_index + shift
                if target_index >= target_length:
                    break
                translated[target_index] = int(value)
            row[field] = translated
        row["stap_colors"] = [
            [int(index) + shift, int(color)]
            for index, color in row.get("stap_colors", [])
            if int(index) + shift < target_length]
    source["num_bases"] = target_length
    return source










def _scaffold_components(rows):
    nodes = {
        (number, index)
        for number, row in rows.items()
        for index, record in enumerate(row["scaf"])
        if record != [-1, -1, -1, -1]}
    adjacency = {node: [] for node in nodes}
    for number, index in nodes:
        record = rows[number]["scaf"][index]
        for offset in (0, 2):
            other = (int(record[offset]), int(record[offset + 1]))
            if other in nodes:
                adjacency[(number, index)].append(other)
    components = []
    visited = set()
    for first in sorted(nodes):
        if first in visited:
            continue
        component = set()
        stack = [first]
        visited.add(first)
        while stack:
            node = stack.pop()
            component.add(node)
            for other in adjacency[node]:
                if other not in visited:
                    visited.add(other)
                    stack.append(other)
        components.append((component, adjacency))
    return components


def _component_helix_ranges(components):
    """Describe the immutable scaffold bands from their actual nodes.

    This intentionally derives the two Path-view bands from the accepted
    template.  No historical edge-balancing table is allowed to redefine or
    repair the Seed.
    """
    result = []
    for component, unused_adjacency in components:
        ranges = {}
        for helix in range(48):
            bases = sorted(base for number, base in component
                           if number == helix)
            if bases:
                ranges[str(helix)] = [bases[0], bases[-1]]
        result.append(ranges)
    return result


def _seed_quota_score(helix_numbers, coordinates, quotas):
    """Return no-bending metrics for one deterministic Seed quota vector."""
    helix_numbers = list(map(int, helix_numbers))
    values = [float(quotas[number]) for number in helix_numbers]
    centre_x = sum(float(coordinates[number][0])
                   for number in helix_numbers)/float(len(helix_numbers))
    centre_y = sum(float(coordinates[number][1])
                   for number in helix_numbers)/float(len(helix_numbers))
    centred = dict((number, (
        float(coordinates[number][0])-centre_x,
        float(coordinates[number][1])-centre_y))
        for number in helix_numbers)
    radius2 = sum(x*x+y*y for x, y in centred.values()) / \
        float(len(helix_numbers))
    radius = math.sqrt(max(radius2, 1e-12))
    average = sum(values)/float(len(values))
    deviations = dict((number, float(quotas[number])-average)
                      for number in helix_numbers)
    x_moment = sum(deviations[number]*centred[number][0]
                   for number in helix_numbers)
    y_moment = sum(deviations[number]*centred[number][1]
                   for number in helix_numbers)
    deviation_total = max(1.0, sum(abs(value)
                                   for value in deviations.values()))
    first_moment = math.hypot(x_moment, y_moment) / (
        radius*deviation_total)
    variance = sum(value*value for value in deviations.values()) / \
        float(len(helix_numbers))
    xx_minus_yy = sum(
        deviations[number]*(centred[number][0]**2-centred[number][1]**2)
        for number in helix_numbers)
    two_xy = sum(
        deviations[number]*2.0*centred[number][0]*centred[number][1]
        for number in helix_numbers)
    anisotropy = math.hypot(xx_minus_yy, two_xy) / (
        max(radius2, 1e-12)*deviation_total)
    return first_moment, variance, anisotropy


def _moment_balanced_seed_targets(magnitude, helix_numbers, coordinates,
                                  capacities=None):
    """Integerize a fractional Seed mean without adding a bend moment.

    The allocation is deterministic.  Every complete common-mode layer is
    assigned to all eligible helices.  For a fractional final layer, the
    ceiling recipients are selected from the real frozen-Seed coordinates by
    minimizing the cross-section first moment, then quota variance and
    second-moment anisotropy.  Capacity-limited helices remain hard limits.
    """
    helix_numbers = sorted(set(map(int, helix_numbers)))
    magnitude = max(0, int(magnitude))
    capacities = (dict((number, magnitude) for number in helix_numbers)
                  if capacities is None else
                  dict((int(number), max(0, int(value)))
                       for number, value in capacities.items()))
    if magnitude > sum(capacities.get(number, 0)
                       for number in helix_numbers):
        raise RuntimeError(
            "Z2安全位点容量不足，无法完成截面平衡整数化。")
    quotas = dict((number, 0) for number in helix_numbers)

    centre_x = sum(float(coordinates[number][0])
                   for number in helix_numbers)/float(len(helix_numbers))
    centre_y = sum(float(coordinates[number][1])
                   for number in helix_numbers)/float(len(helix_numbers))
    angles = dict((number, math.atan2(
        float(coordinates[number][1])-centre_y,
        float(coordinates[number][0])-centre_x))
        for number in helix_numbers)

    def subset_score(subset):
        future = dict(quotas)
        for number in subset:
            future[number] += 1
        return _seed_quota_score(helix_numbers, coordinates, future)

    def improve(initial, eligible):
        selected = set(initial)
        best_score = subset_score(selected)
        while True:
            best_swap = None
            best_key = (best_score, tuple(sorted(selected)))
            for removed in sorted(selected):
                for added in eligible:
                    if added in selected:
                        continue
                    candidate = (selected-{removed}) | {added}
                    key = (subset_score(candidate), tuple(sorted(candidate)))
                    if key < best_key:
                        best_key = key
                        best_swap = (removed, added)
            if best_swap is None:
                return selected, best_score
            selected.remove(best_swap[0])
            selected.add(best_swap[1])
            best_score = best_key[0]

    remaining = magnitude
    while remaining:
        eligible = [number for number in helix_numbers
                    if quotas[number] < capacities.get(number, 0)]
        if not eligible:
            raise RuntimeError(
                "Z2安全位点容量不足，无法完成截面平衡整数化。")
        take = min(remaining, len(eligible))
        if take == len(eligible):
            selected = set(eligible)
        else:
            candidate_sets = []
            exact_search = math.comb(len(eligible), take) <= 5000
            if exact_search:
                candidate_sets.extend(
                    set(values) for values in itertools.combinations(
                        eligible, take))
            else:
                angular = sorted(eligible, key=lambda number: (
                    angles[number], number))
                for phase in range(min(len(angular), 64)):
                    chosen = set()
                    for rank in range(take):
                        target = (phase+(rank+.5)*len(angular)/float(take)) % \
                            len(angular)
                        for offset in range(len(angular)):
                            candidate = angular[int(round(target+offset)) %
                                                len(angular)]
                            if candidate not in chosen:
                                chosen.add(candidate)
                                break
                    candidate_sets.append(chosen)
                candidate_sets.append(set(sorted(eligible, key=lambda number: (
                    quotas[number], angles[number], number))[:take]))
            best = None
            for candidate in candidate_sets:
                if exact_search:
                    improved, score = set(candidate), subset_score(candidate)
                else:
                    improved, score = improve(candidate, eligible)
                key = (score, tuple(sorted(improved)))
                if best is None or key < best[0]:
                    best = (key, improved)
            selected = set(best[1])
        for number in selected:
            quotas[number] += 1
            remaining -= 1

    score = _seed_quota_score(helix_numbers, coordinates, quotas)
    return quotas, {
        "method": "deterministic-first-moment-balanced-integer-rounding",
        "cross_section_first_moment_residual": score[0],
        "cross_section_quota_variance": score[1],
        "cross_section_quota_anisotropy": score[2],
    }


def _apply_frozen_seed_z2_indels(payload, layout):
    """Place the requested integer indels without changing Seed topology.

    The nominal Z2 coordinate band comes from the shared SST/Seed geometry.
    A usable site must be duplex-occupied and strictly longitudinal in both
    polymers, so it cannot be a nick or a scaffold/staple crossover.  The
    requested mean is converted to one integer total across 48 Seed helices
    and spread across helices and Z2 coordinates.  Both insertions and
    deletions are limited to at most three edits in each nominal 8-bp domain
    on each helix.  Insertions are additionally constrained by the 7557-nt
    capacity of each scaffold.
    """
    requested = float(layout.get("mean_indel_per_helix", 0.0) or 0.0)
    nominal_z2_bp = int(layout.get("z2_bp", 32) or 0)
    minimum_deletion = minimum_seed_deletion_per_helix(nominal_z2_bp)
    if requested < minimum_deletion - 1e-9:
        raise RuntimeError(
            "Mean deletion %.1f/helix exceeds the %.1f/helix limit for "
            "%d-bp spacing. Each 8-bp domain permits at most %d evenly "
            "distributed deletions." % (
                requested, minimum_deletion, nominal_z2_bp,
                MAX_SEED_DELETIONS_PER_DOMAIN))
    z2_range = layout.get("seed_z2_range") or layout.get("spacing_range")
    if not z2_range:
        layer_ranges = layout.get("layer_ranges") or []
        if len(layer_ranges) >= 2 and all(
                isinstance(value, (list, tuple)) and len(value) == 2
                for value in layer_ranges[:2]):
            z2_range = [
                int(layer_ranges[0][1]) + 1,
                int(layer_ranges[1][0]) - 1,
            ]
            layout["spacing_range"] = list(z2_range)
            layout["seed_z2_range"] = list(z2_range)
    if not z2_range or len(z2_range) != 2:
        raise RuntimeError("The nominal Seed Z2 coordinate range is missing.")
    z2_low, z2_high = map(int, z2_range)
    if nominal_z2_bp != max(0, z2_high-z2_low+1):
        raise RuntimeError(
            "The nominal Seed Z2 range does not match the selected spacing.")
    # New SST generators always expose ``spacing_range`` while older saved
    # projects may expose only ``seed_z2_range``.  Normalize both aliases
    # before the indel audit is copied into the scaffold-review payload.
    layout["spacing_range"] = [z2_low, z2_high]
    layout["seed_z2_range"] = [z2_low, z2_high]
    total = int(round(requested * 48.0))
    rows = {int(row["num"]): row for row in payload["vstrands"]
            if 0 <= int(row["num"]) < 48}
    components = _scaffold_components(rows)
    labels = {node: component_index
              for component_index, (component, unused_adjacency)
              in enumerate(components) for node in component}

    native_lengths = [len(component) for component, unused in components]
    if total == 0:
        layout.update({
            "mean_indel_per_helix_requested": requested,
            "mean_indel_per_helix_actual": 0.0,
            "actual_z2_spacing_bp": float(nominal_z2_bp),
            "seed_z2_indel_range": [z2_low, z2_high],
            "seed_z2_indel_placements": [],
            "seed_z2_domain_bp": SQUARE_DOMAIN_BP,
            "maximum_seed_indels_per_domain":
                MAX_SEED_DELETIONS_PER_DOMAIN,
            "maximum_seed_deletions_per_domain":
                MAX_SEED_DELETIONS_PER_DOMAIN,
            "minimum_seed_deletion_per_helix": minimum_deletion,
            "seed_scaffold_lengths_after_indel": native_lengths,
            "seed_z2_indel_distribution": {
                "method": (
                    "equal-partition centres with scaffold-capacity repair"),
                "cross_section_allocation_method":
                    "deterministic-first-moment-balanced-integer-rounding",
                "cross_section_first_moment_residual": 0.0,
                "cross_section_quota_variance": 0.0,
                "cross_section_quota_anisotropy": 0.0,
                "per_helix_counts": {
                    str(helix): 0 for helix in range(48)},
                "maximum_distance_from_bin_center": 0.0,
                "scaffold_capacity_nt": 7557,
                "scaffold_lengths_nt": native_lengths,
            },
            "seed_indel_policy": (
                "Z2 only; duplex longitudinal bases only; avoids every "
                "nick and scaffold/staple crossover; evenly distributed; "
                "at most 3 insertions or deletions per nominal 8-bp domain "
                "per helix; "
                "each scaffold <=7557 nt"),
        })
        return native_lengths

    # Scaffold-review files deliberately hide staples.  Candidate safety is
    # nevertheless defined by the immutable template staple topology, not by
    # the presentation stage of the output file.
    topology_rows = {
        int(row["num"]): row
        for row in _translated_frozen_reference(layout)["vstrands"]
        if 0 <= int(row["num"]) < 48}

    def longitudinal(record, helix, base):
        if record == [-1, -1, -1, -1]:
            return False
        for offset in (0, 2):
            partner, partner_base = map(int, record[offset:offset + 2])
            if partner < 0 or partner != helix or \
                    abs(partner_base - base) != 1:
                return False
        return True

    candidates = defaultdict(list)
    for helix in range(48):
        row = rows[helix]
        for base in range(z2_low, z2_high + 1):
            node = (helix, base)
            if node not in labels:
                continue
            if not longitudinal(row["scaf"][base], helix, base) or \
                    not longitudinal(
                        topology_rows[helix]["stap"][base], helix, base):
                continue
            if int(row["loop"][base]) or int(row["skip"][base]):
                continue
            candidates[helix].append((base, labels[node]))
    if any(not candidates[helix] for helix in range(48)):
        missing = [helix for helix in range(48) if not candidates[helix]]
        raise RuntimeError("固定Seed Z2缺少安全indel位点：%s。" % missing)

    capacities = [7557 - length for length in native_lengths]
    if total > sum(capacities):
        raise RuntimeError(
            "Z2 insertion需要%d nt，但两条固定scaffold仅剩%d nt容量；"
            "每条scaffold必须≤7557 nt。" % (total, sum(capacities)))
    if total < 0 and -total > sum(len(values)
                                  for values in candidates.values()):
        raise RuntimeError("Z2安全deletion位点不足：需要%d，只有%d。" % (
            -total, sum(len(values) for values in candidates.values())))

    def ideal_centres(required):
        """Return one target near the centre of every equal Z2 partition."""
        if required <= 0:
            return []
        width = float(nominal_z2_bp) / float(required)
        result = []
        for rank in range(required):
            centre = float(z2_low) + (rank + 0.5) * width - 0.5
            centre = max(float(z2_low), min(float(z2_high), centre))
            partition_start = z2_low + int(math.floor(
                rank*nominal_z2_bp/float(required)))
            partition_end = z2_low + int(math.floor(
                (rank+1)*nominal_z2_bp/float(required))) - 1
            result.append({
                "target": centre,
                "partition": rank,
                "partition_start": partition_start,
                "partition_end": partition_end,
                "domain": min(
                    max(0, int((centre-z2_low) // SQUARE_DOMAIN_BP)),
                    max(0, int(math.ceil(
                        nominal_z2_bp / float(SQUARE_DOMAIN_BP))) - 1)),
            })
        return result

    def ordered_options(values, ideals, tracked_component=None):
        """Minimum-cost ordered site selections for one helix.

        The returned dictionary is keyed by the number of sites placed on
        ``tracked_component``.  Keeping every feasible key allows the later
        global pass to enforce scaffold capacity before it optimizes spatial
        uniformity.
        """
        required = len(ideals)
        if not required:
            # Keep the zero-selection cost structurally identical to every
            # nonzero dynamic-programming state.  A positive structure-wide
            # insertion total below 48 gives some Seed helices one insertion
            # and the remainder zero; the global capacity optimizer always
            # combines three-component lexicographic cost tuples.
            return {0: ((0, 0, 0.0), [])}
        values = sorted(values)
        if len(values) < required:
            return {}
        # Cost is lexicographic: leaving the equal partition is most costly,
        # then leaving the target absolute domain, then squared axial error.
        # This makes partition/domain intersection a real placement rule
        # rather than merely a soft nearest-centre tendency.
        states = {(0, 0): ((0, 0, 0.0), [])}
        for base, component_index in values:
            updated = dict(states)
            for (selected, component_count), (cost, chosen) in states.items():
                if selected >= required:
                    continue
                next_count = component_count + int(
                    tracked_component is not None and
                    component_index == tracked_component)
                # Squared deviation strongly avoids sacrificing one bin while
                # several other bins happen to be exact.  Capacity remains a
                # hard constraint in the global pass below.
                ideal = ideals[selected]
                domain = int((base-z2_low) // SQUARE_DOMAIN_BP)
                # The 8-bp domain capacity is a symmetric physical limit:
                # insertions and deletions may each occupy at most three
                # legal sites in one native domain on one Seed helix.
                if sum(
                        int((chosen_base-z2_low) // SQUARE_DOMAIN_BP) ==
                        domain
                        for chosen_base, unused_component in chosen) >= \
                        MAX_SEED_DELETIONS_PER_DOMAIN:
                    continue
                outside_partition = int(not (
                    ideal["partition_start"] <= base <=
                    ideal["partition_end"]))
                outside_domain = int(domain != ideal["domain"])
                deviation = float(base) - ideal["target"]
                next_cost = (
                    cost[0] + outside_partition,
                    cost[1] + outside_domain,
                    cost[2] + deviation * deviation)
                key = (selected + 1, next_count)
                previous = updated.get(key)
                candidate_path = chosen + [(base, component_index)]
                if previous is None or next_cost < previous[0] or (
                        next_cost == previous[0] and
                        candidate_path < previous[1]):
                    updated[key] = (next_cost, candidate_path)
            states = updated
        return {
            component_count: value
            for (selected, component_count), value in states.items()
            if selected == required
        }

    def choose_domain_limited_sites(values, ideals, domain_count):
        """Choose deletion sites with bin-first/domain-second priorities."""
        remaining = sorted(values)
        room = dict((domain, MAX_SEED_DELETIONS_PER_DOMAIN)
                    for domain in range(domain_count))
        chosen = {}
        unresolved = []

        def options_for(ideal, inside_partition):
            result = []
            for base, component_index in remaining:
                domain = int((base-z2_low) // SQUARE_DOMAIN_BP)
                if not 0 <= domain < domain_count or room[domain] <= 0:
                    continue
                inside = (ideal["partition_start"] <= base <=
                          ideal["partition_end"])
                if inside == inside_partition:
                    result.append((base, component_index, domain))
            return result

        def reserve(ideal, options):
            selected = min(options, key=lambda item: (
                item[2] != ideal["domain"],
                abs(float(item[0])-ideal["target"]),
                abs(item[2]-ideal["domain"]), item[0]))
            base, component_index, domain = selected
            remaining.remove((base, component_index))
            room[domain] -= 1
            chosen[ideal["partition"]] = selected

        for ideal in ideals:
            options = options_for(ideal, True)
            if options:
                reserve(ideal, options)
            else:
                unresolved.append(ideal)
        for ideal in unresolved:
            options = options_for(ideal, False)
            if not options:
                return None
            reserve(ideal, options)
        return [chosen[rank] for rank in range(len(ideals))]

    placements = []
    domain_count = nominal_z2_bp // SQUARE_DOMAIN_BP
    maximum_per_helix_from_domains = (
        domain_count * MAX_SEED_DELETIONS_PER_DOMAIN)
    coordinates = dict((helix, (
        float(rows[helix].get("col", 0)),
        float(rows[helix].get("row", 0)))) for helix in range(48))
    helix_capacities = dict((helix, min(
        maximum_per_helix_from_domains, len(candidates[helix])))
        for helix in range(48))
    targets, cross_section_balance = _moment_balanced_seed_targets(
        abs(total), range(48), coordinates, helix_capacities)
    if any(required > maximum_per_helix_from_domains
           for required in targets.values()):
        operation = "insertions" if total > 0 else "deletions"
        raise RuntimeError(
            "The requested %s cannot be distributed across %d nominal "
            "8-bp Z2 domains with at most %d edits per domain per helix."
            % (operation, domain_count,
               MAX_SEED_DELETIONS_PER_DOMAIN))
    if total < 0:
        if domain_count <= 0:
            raise RuntimeError(
                "Deletions cannot be placed when nominal spacing is 0 bp.")

        for helix in range(48):
            required = targets[helix]
            ideals = ideal_centres(required)
            selected_values = choose_domain_limited_sites(
                candidates[helix], ideals, domain_count)
            if selected_values is None:
                raise RuntimeError(
                    "Seed helix %d cannot distribute deletions across equal "
                    "Z2 partitions while retaining the +/-3 per-domain "
                    "limit and all safe-site rules." % helix)
            for ideal, (base, component_index, domain) in zip(
                    ideals, selected_values):
                rows[helix]["skip"][base] = -1
                placements.append({
                    "helix": helix, "base": base, "value": -1,
                    "scaffold_component": component_index,
                    "z2_domain": domain + 1,
                    "z2_bin": ideal["partition"] + 1,
                    "z2_bin_count": required,
                    "ideal_base": round(ideal["target"], 3),
                    "distance_from_bin_center": round(
                        abs(float(base)-ideal["target"]), 3),
                })
    else:
        if any(capacity < 0 for capacity in capacities):
            raise RuntimeError(
                "固定Seed scaffold已超过7557 nt，不能继续加入insertion。")
        if len(capacities) not in (1, 2):
            raise RuntimeError(
                "固定Seed必须含1或2条scaffold，当前检测到%d条。" %
                len(capacities))

        per_helix_options = {}
        per_helix_ideals = {}
        for helix in range(48):
            required = targets[helix]
            ideals = ideal_centres(required)
            per_helix_ideals[helix] = ideals
            options = ordered_options(
                candidates[helix], ideals,
                tracked_component=0 if len(capacities) == 2 else None)
            if not options:
                raise RuntimeError(
                    "Seed helix %d lacks enough legal Z2 insertion sites."
                    % helix)
            per_helix_options[helix] = options

        # Globally choose one option per helix.  The state tracks component-0
        # use; component-1 use follows from the fixed total.  Capacity is
        # filtered before spatial cost, making <=7557 nt a hard constraint.
        global_states = {0: ((0, 0, 0.0), [])}
        for helix in range(48):
            next_states = {}
            for used_component_zero, (global_cost, choices) \
                    in global_states.items():
                for component_zero_count, (local_cost, selected_values) \
                        in per_helix_options[helix].items():
                    next_zero = used_component_zero + component_zero_count
                    if len(capacities) == 2 and \
                            next_zero > capacities[0]:
                        continue
                    key = next_zero
                    candidate = (
                        tuple(global_cost[index] + local_cost[index]
                              for index in range(3)),
                        choices + [(helix, selected_values)])
                    previous = next_states.get(key)
                    if previous is None or candidate[0] < previous[0]:
                        next_states[key] = candidate
            global_states = next_states
            if not global_states:
                raise RuntimeError(
                    "The requested insertions cannot satisfy scaffold "
                    "capacity while remaining on legal Seed Z2 sites.")

        feasible = []
        for used_component_zero, state in global_states.items():
            if len(capacities) == 1:
                if total <= capacities[0]:
                    feasible.append((used_component_zero, state))
            else:
                used_component_one = total - used_component_zero
                if used_component_one <= capacities[1]:
                    feasible.append((used_component_zero, state))
        if not feasible:
            raise RuntimeError(
                "The requested insertions cannot satisfy the <=7557-nt "
                "capacity of every scaffold.")
        unused_count, (unused_cost, chosen_by_helix) = min(
            feasible, key=lambda item: (item[1][0], item[0]))

        for helix, selected_values in chosen_by_helix:
            ideals = per_helix_ideals[helix]
            for rank, ((base, component_index), ideal) in enumerate(
                    zip(selected_values, ideals)):
                rows[helix]["loop"][base] += 1
                placements.append({
                    "helix": helix, "base": base, "value": 1,
                    "scaffold_component": component_index,
                    "z2_domain": (
                        (base-z2_low) // SQUARE_DOMAIN_BP) + 1,
                    "z2_bin": rank + 1,
                    "z2_bin_count": len(ideals),
                    "ideal_base": round(ideal["target"], 3),
                    "distance_from_bin_center": round(
                        abs(float(base)-ideal["target"]), 3),
                })

    actual_lengths = list(native_lengths)
    for item in placements:
        actual_lengths[item["scaffold_component"]] += item["value"]
    if max(actual_lengths) > 7557:
        raise RuntimeError("Z2 indel后scaffold超过7557 nt：%s。" %
                           actual_lengths)
    mean_actual = sum(item["value"] for item in placements) / 48.0
    layout.update({
        "mean_indel_per_helix_requested": requested,
        "mean_indel_per_helix_actual": mean_actual,
        "actual_z2_spacing_bp": float(nominal_z2_bp) + mean_actual,
        "seed_z2_indel_range": [z2_low, z2_high],
        "seed_z2_indel_placements": placements,
        "seed_z2_domain_bp": SQUARE_DOMAIN_BP,
        "maximum_seed_indels_per_domain":
            MAX_SEED_DELETIONS_PER_DOMAIN,
        "maximum_seed_deletions_per_domain":
            MAX_SEED_DELETIONS_PER_DOMAIN,
        "minimum_seed_deletion_per_helix": minimum_deletion,
        "seed_scaffold_lengths_after_indel": actual_lengths,
        "seed_z2_indel_distribution": {
            "method": (
                "equal-partition/native-domain intersection first; nearest "
                "safe-site and scaffold-capacity repair; no forced stagger"),
            "per_helix_counts": {
                str(helix): targets[helix] for helix in range(48)},
            "maximum_distance_from_bin_center": max(
                (item.get("distance_from_bin_center", 0.0)
                 for item in placements), default=0.0),
            "scaffold_capacity_nt": 7557,
            "scaffold_lengths_nt": actual_lengths,
            "cross_section_allocation_method":
                cross_section_balance["method"],
            "cross_section_first_moment_residual":
                cross_section_balance[
                    "cross_section_first_moment_residual"],
            "cross_section_quota_variance":
                cross_section_balance["cross_section_quota_variance"],
            "cross_section_quota_anisotropy":
                cross_section_balance[
                    "cross_section_quota_anisotropy"],
        },
        "seed_indel_policy": (
            "Z2 only; duplex longitudinal bases only; avoids every nick and "
            "scaffold/staple crossover; one site near every equal-partition "
            "centre; insertions and deletions have at most 3 per nominal "
            "8-bp domain per helix; scaffold capacity is enforced before "
            "spatial uniformity; "
            "each scaffold <=7557 nt"),
    })
    return actual_lengths
























def _frozen_two_layer_seed_scaffold_payload(layout):
    """Load the accepted two-layer routing without synthesising new paths.

    The physical Seed is immutable and independent of the current SST length
    or spacing.  Scaffold review intentionally shows only its fixed scaffold;
    finalization restores the accepted staple/capture template verbatim.
    """
    source = _translated_frozen_reference(layout)
    canvas_shift = _seed_canvas_shift(layout)
    array_length = max(
        int(source.get("num_bases", 0) or 0),
        max(len(row.get("scaf", ())) for row in source["vstrands"]),
        int(layout.get("array_length", SST_ARRAY_LENGTH)))
    rows = []
    for source_row in source["vstrands"]:
        number = int(source_row["num"])
        if not 0 <= number < 48:
            continue
        row = dict(source_row)
        _resize_row(row, array_length)
        row["scaf"] = [list(record) for record in row["scaf"]]
        row["stap"] = [[-1, -1, -1, -1] for unused in range(array_length)]
        row["stap_colors"] = []
        row["loop"] = [0] * array_length
        row["skip"] = [0] * array_length
        rows.append(row)
    by_number = {int(row["num"]): row for row in rows}
    components = _scaffold_components(by_number)
    lengths = sorted(len(component) for component, unused in components)
    if lengths != [7300, 7336]:
        raise RuntimeError(
            "冻结2L Seed routing已被改动：scaffold lengths=%r。" % lengths)
    ranges = {}
    for number, row in by_number.items():
        occupied = [index for index, record in enumerate(row["scaf"])
                    if record != [-1, -1, -1, -1]]
        ranges[number] = (min(occupied), max(occupied))
    low_stagger = max(low for low, unused in ranges.values()) - min(
        low for low, unused in ranges.values())
    high_stagger = max(high for unused, high in ranges.values()) - min(
        high for unused, high in ranges.values())
    if low_stagger > 11 or high_stagger > 21:
        raise RuntimeError(
            "冻结2L Seed边缘范围异常：low=%d high=%d。" %
            (low_stagger, high_stagger))
    total_nt = sum(lengths)
    capacity = scaffold_capacity_plan(total_nt)
    if capacity["count"] != 2 or max(lengths) > 7557:
        raise RuntimeError("冻结2L Seed不再满足两条scaffold容量。")
    # The reference physical ranges are the accepted result.  Report them in
    # the layout so all later capture and validation code uses the same truth.
    # The frozen Seed is defined in the coordinate frame of the accepted
    # 128/32/128 Square 2L reference.  A long SST may move the whole Seed by
    # one or more complete 32-bp periods; that global translation does not
    # alter any routing relationship inside the reference.
    first_end = 48 + canvas_shift + 128 - 1
    second_start = first_end + 1 + 32
    capture_helices = tuple(range(0, 8)) + tuple(range(24, 32))
    physical_z1 = max(first_end - ranges[number][0] + 1
                      for number in capture_helices)
    physical_z3 = max(ranges[number][1] - second_start + 1
                      for number in capture_helices)
    square_geometry = layout.get("square_centered_geometry", {})
    capture_support_ranges = _seed_capture_support_ranges(layout)
    overlap = fixed_seed_overlap_layout(
        layout.get("layer_ranges", [
            [48 + canvas_shift, 175 + canvas_shift],
            [208 + canvas_shift, 335 + canvas_shift]]),
        (layout.get("theoretical_capture_positions_by_layer")
         if layout.get("lattice_type") in
         ("kagome", "square_kagome") else None),
        lattice_type=layout.get("lattice_type", "square"),
        seed_layer_ranges=capture_support_ranges,
        seed_capture_positions_by_layer=(
            layout.get("seed_capture_positions_by_layer")
            if layout.get("lattice_type") in
            ("kagome", "square_kagome") else
            square_geometry.get("seed_capture_positions_by_layer")))
    actual_z1, actual_z3 = map(int, layout.get(
        "seed_sst_overlap_bp", overlap["seed_sst_overlap_bp"]))
    component_ranges = _component_helix_ranges(components)
    layout.update({
        "seed_z1_actual_bp": actual_z1,
        "seed_z3_actual_bp": actual_z3,
        "seed_z1_edge_growth_bp": 0,
        "seed_z3_edge_growth_bp": 0,
        "seed_edge_stagger_limit_used_bp": max(low_stagger, high_stagger),
        "seed_edge_selected_scaffold_count": 2,
        "capture_helices_define_maximum_actual_length": True,
        "seed_routing_reference": "Square_Seed_2L_newtemplate.json",
        "seed_routing_is_frozen_reference": True,
        "seed_template_physical_support_bp": [physical_z1, physical_z3],
    })
    payload = {
        "name": "square_moire_seed_scaffold.json",
        "vstrands": [by_number[number] for number in range(48)],
        "num_bases": array_length,
        "lattice": "square",
        "moire_edge_metadata": {
            "scope": "Moiré Designer only; cadnano AutoCS unchanged",
            "partition": (
                "frozen two-layer left/right capacity-safe Path-view "
                "bands; all 48 helices"),
            "policy": (
                "immutable Square_Seed_2L_newtemplate.json routing; no dynamic "
                "crossover synthesis"),
            "reference": "Square_Seed_2L_newtemplate.json",
            "coordinate_shift_bp": canvas_shift,
            "total_seed_scaffold_nt": total_nt,
            "balanced_scaffold_lengths": lengths,
            "maximum_edge_stagger_bp": max(low_stagger, high_stagger),
            "edge_stagger_limit_used_bp": max(low_stagger, high_stagger),
            "selected_scaffold_count": 2,
            "capture_helices_define_maximum_actual_length": True,
            "seed_z1_actual_bp": actual_z1,
            "seed_z3_actual_bp": actual_z3,
            "seed_template_physical_support_bp": [physical_z1, physical_z3],
            "seed_z1_edge_growth_bp": layout["seed_z1_edge_growth_bp"],
            "seed_z3_edge_growth_bp": layout["seed_z3_edge_growth_bp"],
            "seed_routing_is_frozen_reference": True,
            "ranges": {str(number): list(value)
                       for number, value in ranges.items()},
            # Audit-only description of the accepted components.  It is
            # derived from the template and is never used to regenerate it.
            "band_ranges": component_ranges,
            "dynamic_seed_adjustments": False,
        },
    }
    return payload, lengths












def _reference_seed_scaffold_payload(layout=None):
    """Return the one immutable two-layer Seed scaffold template.

    SST length and spacing affect only SST placement and the reported
    overlap.  They never select a Seed crop/growth/seam/AutoCS branch.
    """
    layout = layout or structure_layout()
    layout["seed_z1_requested_bp"] = 128
    layout["seed_z3_requested_bp"] = 128
    layout["seed_length_adjustment_enabled"] = False
    layout["seed_geometry_policy"] = "immutable_2L_reference"
    return _frozen_two_layer_seed_scaffold_payload(layout)


def _seed_routing_layout(sst_layout):
    """Return the native 8x8 Seed routing inputs for either SST lattice.

    Kagome SST owns a different active-helix mask and capture catalogue, but
    the Seed is still the validated Square S8-R4x4C object.  Rebuild only the
    Seed edge/routing fields through the normal Square layout constructor;
    lattice-specific capture assignments remain in the untouched SST payload
    and are resolved later by ``_combine_seed_sst``.
    """
    routing = copy.deepcopy(sst_layout)
    routing["seed_z1_requested_bp"] = 128
    routing["seed_z3_requested_bp"] = 128
    routing["seed_z1_actual_bp"] = 128
    routing["seed_z3_actual_bp"] = 128
    routing["seed_z1_edge_growth_bp"] = 0
    routing["seed_z3_edge_growth_bp"] = 0
    routing["seed_length_adjustment_enabled"] = False
    routing["seed_geometry_policy"] = "immutable_2L_reference"
    routing["seed_cross_section_preset"] = "s8_r4x4"
    return routing


def _resize_row(row, length):
    for field in ("scaf", "stap"):
        records = list(row.get(field, []))[:length]
        records.extend([[-1, -1, -1, -1] for unused in
                        range(length - len(records))])
        row[field] = records
    for field in ("loop", "skip"):
        values = list(row.get(field, []))[:length]
        values.extend([0] * (length - len(values)))
        row[field] = values
    row["stap_colors"] = [
        item for item in row.get("stap_colors", [])
        if int(item[0]) < length]


def _combine_seed_sst(seed_payload, stage, sst_payload=None):
    sst = (sst_payload if sst_payload is not None else
           build_shifted_sst_payload(
               "square_moire_sst_2L.json", reserve_capture_gaps=True))
    if not sst.get("moire_structure_metadata", {}).get(
            "capture_gaps_reserved", False):
        sst = build_capture_ready_sst_payload(
            sst, "%s_capture_ready.json" %
            Path(sst.get("name", "square_moire_sst_2L.json")).stem)
    sst = payload_to_internal_numbering(sst)
    inherited_metadata = sst.get("moire_structure_metadata", {})
    configured_layout = sst.get("moire_structure_metadata", {}).get(
        "variable_length_layout", {})
    # Rebuild derived Capture assignments with the current policy.  Saved
    # scaffold-review files may contain the former A0/B0/A1/B1 physical cycle;
    # reusing that cached dictionary would silently preserve the defect when
    # the user regenerates Staple/Capture from an older accepted scaffold.
    inherited_lattice = inherited_metadata.get("lattice_type")
    if inherited_lattice in ("kagome", "square_kagome"):
        centered_geometry = configured_layout.get(
            "square_centered_geometry",
            configured_layout.get("centered_geometry", {}))
        coordinate_shift = int(configured_layout.get(
            "coordinate_shift_bp",
            centered_geometry.get("coordinate_shift_bp", 0)) or 0)
        seed_layer_ranges = copy.deepcopy(configured_layout.get(
            "seed_layer_ranges",
            centered_geometry.get("seed_layer_ranges", [
                [48 + coordinate_shift, 175 + coordinate_shift],
                [208 + coordinate_shift, 335 + coordinate_shift]])))
        capture_support_ranges = _seed_capture_support_ranges(
            configured_layout)
        layout = copy.deepcopy(configured_layout)
        layout.update({
            "lattice_type": inherited_lattice,
            "coordinate_shift_bp": coordinate_shift,
            "capture_extension_nt": int(configured_layout.get(
                "capture_extension_nt", 32)),
            "seed_z1_requested_bp": 128,
            "seed_z3_requested_bp": 128,
            "seed_z1_actual_bp": 128,
            "seed_z3_actual_bp": 128,
            "spacing_range": [
                int(configured_layout["layer_ranges"][0][1]) + 1,
                int(configured_layout["layer_ranges"][1][0]) - 1],
            "seed_layer_ranges": seed_layer_ranges,
            "capture_support_ranges": copy.deepcopy(
                capture_support_ranges),
            "seed_capture_positions_by_layer": copy.deepcopy(
                configured_layout.get(
                    "seed_capture_positions_by_layer", [[], []])),
            "capture_positions": [
                value for layer in configured_layout.get(
                    "theoretical_capture_positions_by_layer", [])
                for value in layer],
            "kagome_capture_anchor_assignments": copy.deepcopy(
                configured_layout.get(
                    "kagome_capture_anchor_assignments",
                    inherited_metadata.get(
                        "kagome_capture_anchor_assignments_sst_only", []))),
            "auxiliary_sst_routing": copy.deepcopy(
                configured_layout.get(
                    "auxiliary_sst_routing",
                    inherited_metadata.get("auxiliary_sst_routing", {}))),
        })
        overlap = fixed_seed_overlap_layout(
            configured_layout["layer_ranges"],
            configured_layout.get(
                "theoretical_capture_positions_by_layer"),
            lattice_type=inherited_lattice,
            seed_layer_ranges=capture_support_ranges,
            seed_capture_positions_by_layer=configured_layout.get(
                "seed_capture_positions_by_layer"))
        layout.update(overlap)
        layout["layer_ranges"] = copy.deepcopy(
            configured_layout["layer_ranges"])
        # ``fixed_seed_overlap_layout`` reports the mask supplied to it as
        # ``seed_layer_ranges``.  Keep that actual contact mask in its own
        # field, while preserving the nominal frozen-Seed ranges used for
        # translating and validating the routing template.
        layout["capture_support_ranges"] = copy.deepcopy(
            overlap["overlap_ranges"])
        layout["seed_layer_ranges"] = copy.deepcopy(seed_layer_ranges)
        # Capture assignment uses every physically legal Kagome column in
        # the current SST duplex.  The reported Z1/Z3 values are a different
        # quantity: the two actual overlaps from the fixed 288-bp Seed
        # partition.  Preserve the generator's overlap-optimized placement
        # metadata so reopening or regenerating a design cannot replace it
        # with the full SST layer lengths.
        layout["overlap_ranges"] = copy.deepcopy(
            configured_layout.get(
                "overlap_ranges", overlap["overlap_ranges"]))
        layout["seed_sst_overlap_bp"] = list(
            map(int, configured_layout.get(
                "seed_sst_overlap_bp", overlap["seed_sst_overlap_bp"])))
        layout["capture_positions"] = [
            value for values in layout["capture_positions_by_layer"]
            for value in values]
        layout["capture_column_count"] = len(layout["capture_positions"])
        layout["capture_columns_by_layer"] = [
            len(layer)
            for layer in layout.get("capture_positions_by_layer", [])]
        layout["capture_sites_per_column_by_layer"] = copy.deepcopy(
            overlap["capture_sites_per_column_by_layer"])
        layout["capture_sites_by_layer"] = list(
            overlap["capture_sites_by_layer"])
        layout["capture_pair_equivalents_by_layer"] = list(
            overlap["capture_pair_equivalents_by_layer"])
        layout["capture_pair_equivalents"] = overlap[
            "capture_pair_equivalents"]
        layout["pair_count_by_layer"] = [
            int(math.ceil(len(layer) / 2.0))
            for layer in layout.get("capture_positions_by_layer", [])]
        layout["pair_count"] = sum(layout["pair_count_by_layer"])
        layout["capture_site_assignments"] = capture_site_assignments(layout)
        layout["capture_export_site_assignments"] = \
            capture_export_site_assignments(layout)
        layout["expected_capture_bridges"] = sum(
            len(item["bridges"])
            for item in layout["capture_site_assignments"])
        layout["expected_capture_export_sequences"] = sum(
            len(item["bridges"])
            for item in layout["capture_export_site_assignments"])
    else:
        layout = structure_layout(
            configured_layout.get("z1_bp", 128),
            configured_layout.get("z2_bp", 32),
            configured_layout.get("z3_bp", 128),
            configured_layout.get("seed_z1_requested_bp"),
            configured_layout.get("seed_z3_requested_bp"),
            configured_layout.get("capture_extension_nt",
                                  CAPTURE_EXTENSION_NT))
        layout["mean_indel_per_helix"] = float(
            configured_layout.get("mean_indel_per_helix", 0.0) or 0.0)
        layout["auxiliary_sst_routing"] = copy.deepcopy(
            configured_layout.get(
                "auxiliary_sst_routing",
                inherited_metadata.get("auxiliary_sst_routing", {})))
    for key in (
            "mean_indel_per_helix_requested",
            "mean_indel_per_helix_actual",
            "actual_z2_spacing_bp",
            "seed_z2_indel_range",
            "seed_z2_indel_placements",
            "seed_z2_domain_bp",
            "maximum_seed_indels_per_domain",
            "maximum_seed_deletions_per_domain",
            "minimum_seed_deletion_per_helix",
            "seed_scaffold_lengths_after_indel",
            "seed_z2_indel_distribution",
            "seed_indel_policy"):
        if key in configured_layout:
            layout[key] = copy.deepcopy(configured_layout[key])
    edge_metadata = seed_payload.get("moire_edge_metadata", {})
    # Seed topology is immutable.  Legacy edge optimizer metadata is never
    # imported back into a newly generated design.
    layout.update({
        "seed_z1_requested_bp": 128,
        "seed_z3_requested_bp": 128,
        "seed_z1_actual_bp": int(edge_metadata.get(
            "seed_z1_actual_bp", layout.get("seed_z1_actual_bp", 128))),
        "seed_z3_actual_bp": int(edge_metadata.get(
            "seed_z3_actual_bp", layout.get("seed_z3_actual_bp", 128))),
        "seed_z1_edge_growth_bp": 0,
        "seed_z3_edge_growth_bp": 0,
        "seed_length_adjustment_enabled": False,
        "seed_geometry_policy": "immutable_2L_reference",
        "seed_routing_is_frozen_reference": True,
        "seed_routing_reference": "Square_Seed_2L_newtemplate.json",
    })
    seed_preset = configured_layout.get(
        "seed_cross_section_preset",
        inherited_metadata.get("seed_cross_section_preset", "s8_r4x4"))
    layout["seed_cross_section_preset"] = seed_preset
    layout["lattice_type"] = inherited_metadata.get(
        "lattice_type", configured_layout.get("lattice_type", "square"))
    length = int(layout.get("array_length", sst.get("num_bases",
                                                    SST_ARRAY_LENGTH)))
    rows = []
    for row in sorted(seed_payload["vstrands"], key=lambda item: item["num"]):
        # Finalization imports the combined review file, so its encoder also
        # contains the 16 SST helices.  They are replaced below by the exact
        # shifted/capture-aware SST template and must not be duplicated.
        if int(row["num"]) >= 48:
            continue
        row = dict(row)
        _resize_row(row, length)
        if stage == "scaffold_review":
            row["stap"] = [[-1, -1, -1, -1] for unused in range(length)]
            row["stap_colors"] = []
        rows.append(row)
    for row in sst["vstrands"]:
        if int(row["num"]) < 48:
            continue
        row = dict(row)
        _resize_row(row, length)
        rows.append(row)
    result = {
        "name": seed_payload.get("name", "square_moire_seed.json"),
        "vstrands": sorted(rows, key=lambda item: item["num"]),
        "num_bases": length,
        "lattice": "square",
        "scaffold_colors": seed_payload.get("scaffold_colors", []),
        "moire_structure_metadata": {
            "stage": stage,
            "seed_cross_section": "S8-R4x4C (48 helices)",
            "seed_cross_section_preset": seed_preset,
            "sst_cross_section": (
                "Square layer 1 + Kagome layer 2"
                if inherited_metadata.get("lattice_type") ==
                "square_kagome" else
                ("Kagome 12/16 active helices" if
                 inherited_metadata.get("lattice_type") == "kagome" else
                 "Square 4x4 (16 helices)")),
            "lattice_type": inherited_metadata.get(
                "lattice_type", "square"),
            "kagome_capture_anchor_assignments_sst_only": copy.deepcopy(
                inherited_metadata.get(
                    "kagome_capture_anchor_assignments_sst_only", [])),
            "auxiliary_sst_routing": copy.deepcopy(
                inherited_metadata.get(
                    "auxiliary_sst_routing",
                    layout.get("auxiliary_sst_routing", {}))),
            "sst_layers": 2,
            "sst_scaffold_ranges": inherited_metadata.get(
                "sst_scaffold_ranges", layout.get(
                    "scaffold_ranges", layout.get(
                        "sst_scaffold_ranges", layout["layer_ranges"]))),
            "sst_staple_ranges": inherited_metadata.get(
                "sst_staple_ranges", layout["staple_ranges"]),
            "sst_complete_source_unchanged": True,
            "sst_combined_topology": "capture-ready derived copy",
            "output_sst_export": (
                "generated later as a temporary same-coordinate snapshot"),
            "sst_source_base_start": 16,
            "sst_applied_shift_bp": 32,
            "sst_generated_base_start": layout["layer_ranges"][0][0],
            "seed_scaffold_partition": "left/right Path-view bands",
            "seed_edge_routing": seed_payload.get(
                "moire_edge_metadata",
                inherited_metadata.get("seed_edge_routing", {})),
            "seed_adjustment_policy": (
                "immutable accepted 2L template; no trim, growth, edge seam, "
                "AutoCS rerun or staple repair"),
            "capture_anchor_positions": layout["capture_positions"],
            "capture_positions_by_layer": layout[
                "capture_positions_by_layer"],
            "capture_phase_cycle": list(CAPTURE_PHASE_CYCLE),
            "capture_export_phase_cycle": list(CAPTURE_EXPORT_PHASE_CYCLE),
            "capture_face_definitions": [
                {
                    key: (list(value) if isinstance(value, tuple) else value)
                    for key, value in face.items()
                }
                for face in CAPTURE_FACE_DEFINITIONS],
            "capture_site_assignments_internal":
                layout["capture_site_assignments"],
            "capture_export_site_assignments_internal":
                layout["capture_export_site_assignments"],
            "capture_connection_policy": (
                "structure connects only origin face helices 16-19 and "
                "47-44; translated face helices are added during sequence "
                "export"),
            "seed_support_ranges": layout["seed_layer_ranges"],
            "seed_sst_overlap_ranges": layout["overlap_ranges"],
            "spacing_range": layout["spacing_range"],
            "seed_z2_range": layout["spacing_range"],
            "spacing_seed_z2_coincident": True,
            "variable_length_layout": layout,
            "seed_scaffold_count": (
                len(seed_payload.get("moire_edge_metadata", {}).get(
                    "band_ranges", [])) or
                int(inherited_metadata.get("seed_scaffold_count", 2))),
            "sequence_assignment": "pending step 3",
        },
    }
    return result


def _reference_staples_for_layout(layout, capture_translation="origin"):
    source_rows = {int(row["num"]): row for row in _reference()["vstrands"]}
    source_ranges = ((40, 183), (200, 343))
    length = int(layout["array_length"])
    exact_reference = (
        (layout["z1_bp"], layout["z2_bp"], layout["z3_bp"]) ==
        (128, 32, 128) and layout.get("coordinate_shift_bp", 0) == 0)
    result = {}
    for number, source in source_rows.items():
        row = copy.deepcopy(source)
        _resize_row(row, length)
        if not exact_reference:
            records = [[-1, -1, -1, -1] for unused in range(length)]
            for source_range, target_range in zip(
                    source_ranges, layout["staple_ranges"]):
                _copy_periodic_segment(
                    source["stap"], source_range, records, target_range,
                    minimum_internal_helix=0, interior_start=8,
                    repeat_bp=128)
            row["stap"] = records
            source_colors = {int(index): int(color)
                             for index, color in
                             source.get("stap_colors", [])}
            colors = {}
            for source_range, target_range in zip(
                    source_ranges, layout["staple_ranges"]):
                target_length = target_range[1] - target_range[0] + 1
                for offset, target_index in enumerate(
                        range(target_range[0], target_range[1] + 1)):
                    source_index = _periodic_source_index(
                        source_range[0], source_range[1], offset,
                        target_length, interior_start=8, repeat_bp=128)
                    if source_index in source_colors:
                        colors[target_index] = source_colors[source_index]
            row["stap_colors"] = [list(item)
                                  for item in sorted(colors.items())]
        result[number] = row

    # Seed rows are not a variable-length SST routing program.  They must be
    # the accepted 2L Seed template, translated only by a complete 32-bp
    # canvas period.  The periodic construction above remains relevant only
    # to the SST rows.
    shifted_seed_rows = {
        int(row["num"]): row
        for row in _translated_frozen_reference(layout)["vstrands"]
        if 0 <= int(row["num"]) < 48}
    for number, row in shifted_seed_rows.items():
        row = copy.deepcopy(row)
        _resize_row(row, length)
        result[number] = row

    if not exact_reference:
        # A longer SST appends another complete 128-base staple-routing
        # program.  The reference program already contains its own nick
        # pattern; without an explicit boundary nick, the last staple of one
        # program is joined to the first staple of the next and can contain
        # two capture targets (typically a 96-nt product).  Open only that
        # longitudinal boundary.  Crossovers and every capture edge remain
        # untouched.
        for target_low, target_high in layout["staple_ranges"]:
            boundary = int(target_low) + 8 + 128
            while boundary <= int(target_high) - 8:
                for number, row in result.items():
                    left = boundary - 1
                    right = boundary
                    if not (0 <= left < len(row["stap"]) and
                            0 <= right < len(row["stap"])):
                        continue
                    for offset in (0, 2):
                        if row["stap"][left][offset:offset + 2] == [
                                number, right]:
                            row["stap"][left][offset:offset + 2] = [-1, -1]
                    for offset in (0, 2):
                        if row["stap"][right][offset:offset + 2] == [
                                number, left]:
                            row["stap"][right][offset:offset + 2] = [-1, -1]
                boundary += 128

    assignments = capture_site_assignments(layout, capture_translation)
    capture_positions = {item["position"] for item in assignments}
    # Remove copied reference bridges first.  Every realized column is then
    # rebuilt from exactly one A0/B0/A1/B1 mapping, preventing duplicate or
    # branched Seed–SST capture connections.
    for number, row in result.items():
        for position in capture_positions:
            record = row["stap"][position]
            for offset in (0, 2):
                partner = int(record[offset])
                if partner >= 0 and ((number < 48 <= partner) or
                                     (partner < 48 <= number)):
                    record[offset:offset + 2] = [-1, -1]

    # A1/B1 are the 4x4-unit translations of the validated A0/B0 template.
    if layout.get("lattice_type") in ("kagome", "square_kagome"):
        for assignment in assignments:
            position = int(assignment["position"])
            for bridge in assignment["bridges"]:
                output_helix = int(bridge["sst_helix"])
                seed_helix = int(bridge["seed_helix"])
                output_offset = int(bridge["sst_slot"])
                seed_record = result[seed_helix]["stap"][position]
                seed_offsets = [
                    offset for offset in (0, 2)
                    if int(seed_record[offset]) < 0]
                if not seed_offsets:
                    raise RuntimeError(
                        "Kagome capture位点的Seed端没有可用nick："
                        "h%d:%d。" % (seed_helix, position))
                seed_offset = seed_offsets[0]
                result[output_helix]["stap"][position][
                    output_offset:output_offset + 2] = [seed_helix, position]
                result[seed_helix]["stap"][position][
                    seed_offset:seed_offset + 2] = [output_helix, position]
        return result

    template_mappings = {
        phase: dict(CAPTURE_PHASE_MAPPINGS[phase + "0"])
        for phase in ("A", "B")}
    canonical = {"A": 72, "B": 56}
    for assignment in assignments:
        position = int(assignment["position"])
        phase = assignment["phase"]
        source_position = canonical[phase]
        for bridge in assignment["bridges"]:
            output_helix = int(bridge["sst_helix"])
            logical_output_helix = int(bridge.get(
                "logical_sst_helix", output_helix))
            seed_helix = int(bridge["seed_helix"])
            template_seed = int(
                template_mappings[phase][logical_output_helix])
            output_source = source_rows[logical_output_helix]["stap"][
                source_position]
            seed_source = source_rows[template_seed]["stap"][
                source_position]
            output_offset = next(
                offset for offset in (0, 2)
                if int(output_source[offset]) == template_seed)
            seed_offset = next(
                offset for offset in (0, 2)
                if int(seed_source[offset]) == logical_output_helix)
            result[output_helix]["stap"][position][
                output_offset:output_offset + 2] = [seed_helix, position]
            result[seed_helix]["stap"][position][
                seed_offset:seed_offset + 2] = [output_helix, position]
    return result


def _crop_reference_staples(payload, capture_translation="origin"):
    """Install the validated two-layer staple/capture topology.

    Scaffold routing is the accepted/generated routing; staple routing and
    capture colors are the first two layers of the supplied 8x8 Seed
    reference.  Cropping at 40..343 removes the third SST layer cleanly.
    """
    layout = payload.get("moire_structure_metadata", {}).get(
        "variable_length_layout", structure_layout())
    if layout.get("lattice_type") in ("kagome", "square_kagome"):
        rows = {int(row["num"]): row for row in payload["vstrands"]}
        # Export translation is a replacement of the physical K0 bridges,
        # not an addition to them.  Open every existing Seed--SST edge first
        # so the Kagome SST slot is available for its K1 counterpart.  This
        # intentionally leaves the old Seed core unextended in the derived
        # export-only copy; the accepted physical JSON is never modified.
        existing_edges = []
        for number in range(48):
            for index, record in enumerate(rows[number]["stap"]):
                for offset in (0, 2):
                    if int(record[offset]) >= 48:
                        existing_edges.append((number, index, offset))
        for number, index, offset in existing_edges:
            if int(rows[number]["stap"][index][offset]) >= 48:
                _disconnect_staple_slot(rows, number, index, offset)
        installed = 0
        assignments = capture_site_assignments(layout, capture_translation)
        kagome_assignments = [
            item for item in assignments
            if item.get("phase") == "K"]
        square_assignments = [
            item for item in assignments
            if item.get("phase") != "K"]
        for assignment in kagome_assignments:
            position = int(assignment["position"])
            for bridge in assignment["bridges"]:
                output_helix = int(bridge["sst_helix"])
                seed_helix = int(bridge["seed_helix"])
                output_offset = int(bridge["sst_slot"])
                output_record = rows[output_helix]["stap"][position]
                if int(output_record[output_offset]) >= 0:
                    raise RuntimeError(
                        "Kagome SST capture端点未打开：h%d:%d。" %
                        (output_helix, position))
                seed_record = rows[seed_helix]["stap"][position]
                # The translated Kagome core already ends at this capture
                # base in the frozen Seed.  Use its genuinely open record
                # slot first, preserving the longitudinal core segment.  An
                # older branch always cut that segment and produced a
                # one-base singleton which caDNAno could not join to the SST
                # during sequence export.
                seed_offsets = [
                    offset for offset in (0, 2)
                    if int(seed_record[offset]) < 0]
                longitudinal_offsets = [
                    offset for offset in (0, 2)
                    if int(seed_record[offset]) == seed_helix and
                    abs(int(seed_record[offset + 1]) - position) == 1]
                if not seed_offsets:
                    seed_offsets = longitudinal_offsets
                if not seed_offsets:
                    raise RuntimeError(
                        "Kagome capture位点的Seed端没有合法开放槽或纵向"
                        "连接：h%d:%d。" % (seed_helix, position))
                seed_offset = seed_offsets[0]
                if seed_offset in longitudinal_offsets:
                    _disconnect_staple_slot(
                        rows, seed_helix, position, seed_offset)
                rows[seed_helix]["stap"][position][
                    seed_offset:seed_offset + 2] = [output_helix, position]
                output_record[output_offset:output_offset + 2] = [
                    seed_helix, position]
                installed += 1
        if square_assignments:
            # Mixed Square--Kagome designs use the exact Square template
            # slots for layer 1.  Reuse the same canonical A/B lookup as the
            # normal Square branch, but install only the already opened
            # layer-1 bridge endpoints; layer 2 remains purely Kagome.
            source_rows = {int(row["num"]): row
                           for row in _reference()["vstrands"]}
            template_mappings = {
                phase: dict(CAPTURE_PHASE_MAPPINGS[phase + "0"])
                for phase in ("A", "B")}
            canonical = {"A": 72, "B": 56}
            for assignment in square_assignments:
                position = int(assignment["position"])
                phase = assignment["phase"]
                source_position = canonical[phase]
                for bridge in assignment["bridges"]:
                    output_helix = int(bridge["sst_helix"])
                    logical_output = int(bridge.get(
                        "logical_sst_helix", output_helix))
                    seed_helix = int(bridge["seed_helix"])
                    template_seed = int(
                        template_mappings[phase][logical_output])
                    output_source = source_rows[logical_output]["stap"][
                        source_position]
                    seed_source = source_rows[template_seed]["stap"][
                        source_position]
                    output_offset = next(
                        offset for offset in (0, 2)
                        if int(output_source[offset]) == template_seed)
                    seed_offset = next(
                        offset for offset in (0, 2)
                        if int(seed_source[offset]) == logical_output)
                    output_record = rows[output_helix]["stap"][position]
                    if int(output_record[output_offset]) >= 0:
                        _disconnect_staple_slot(
                            rows, output_helix, position, output_offset)
                    seed_record = rows[seed_helix]["stap"][position]
                    if int(seed_record[seed_offset]) >= 0:
                        _disconnect_staple_slot(
                            rows, seed_helix, position, seed_offset)
                    output_record[output_offset:output_offset + 2] = [
                        seed_helix, position]
                    seed_record[seed_offset:seed_offset + 2] = [
                        output_helix, position]
                    installed += 1
        return installed
    reference_rows = _reference_staples_for_layout(
        layout, capture_translation)
    rows = {int(row["num"]): row for row in payload["vstrands"]}
    length = int(layout["array_length"])
    protected_indices = {
        index for low, high in layout["staple_ranges"]
        for index in range(int(low), int(high) + 1)}
    physical_assignments = capture_site_assignments(
        layout, capture_translation)
    physical_targets = {
        (int(bridge["seed_helix"]), int(assignment["position"]))
        for assignment in physical_assignments
        for bridge in assignment["bridges"]}
    physical_sst_targets = {
        (int(bridge["sst_helix"]), int(assignment["position"]))
        for assignment in physical_assignments
        for bridge in assignment["bridges"]}
    for number, row in rows.items():
        logical_number = number - 16 if 64 <= number < 80 else number
        source = reference_rows[logical_number]
        if number < 48:
            # AutoStaple/Autobreak already produced normal staples over every
            # Seed scaffold base.  Replace only the capture-bearing windows;
            # retaining the generated records outside those windows prevents
            # the former 8/16-nt edge fragments and covers Z2/outer routing.
            records = [list(record) for record in row["stap"]]
            # AutoStaple may have connected the temporary generated SST to
            # Seed. Remove every such provisional link first. For a variable
            # design, restore only the actual capture slot below.  The same
            # rule now also applies to the canonical 128-bp design: copying
            # both complete reference windows here used to overwrite the
            # valid cadnano AutoCS-staple topology whenever Scaffold routing
            # changed, causing missing or mismatched staple crossovers.
            for index, record in enumerate(records):
                for offset in (0, 2):
                    if int(record[offset]) >= 48:
                        record[offset:offset + 2] = [-1, -1]
            for helix, index in sorted(physical_targets):
                if helix != number:
                    continue
                reference_record = source["stap"][index]
                for offset in (0, 2):
                    if int(reference_record[offset]) >= 48:
                        # A capture base must be a genuine longitudinal
                        # Seed-staple end.  It may never share a nucleotide
                        # with an AutoCS staple crossover, even if the other
                        # record slot is longitudinal and could be opened.
                        if any(
                                int(records[index][slot]) >= 0 and
                                int(records[index][slot]) != number
                                for slot in (0, 2)):
                            raise RuntimeError(
                                "capture候选h%d[%d]落在AutoCS staple "
                                "crossover碱基上；必须调整capture相位/位点，"
                                "禁止在该base连接。" % (number, index))
                        current_partner = records[index][offset]
                        if int(current_partner) >= 0:
                            if not (
                                    int(current_partner) == number and
                                    abs(int(records[index][offset + 1]) -
                                        index) == 1):
                                raise RuntimeError(
                                    "capture位点h%d[%d]不是合法纵向端点。" %
                                    (number, index))
                            _disconnect_staple_slot(
                                rows, number, index, offset)
                        records[index][offset:offset + 2] = list(
                            reference_record[offset:offset + 2])
            row["stap"] = records
            colors = {
                int(index): int(color)
                for index, color in row.get("stap_colors", [])
                if int(index) not in protected_indices}
            colors.update({
                int(index): int(color)
                for index, color in source.get("stap_colors", [])
                if (number, int(index)) in physical_targets})
            row["stap_colors"] = [list(item)
                                  for item in sorted(colors.items())]
            continue
        # Preserve expert SST edits.  Only restore the Seed-facing half of
        # each validated capture crossover and its pair color.
        for low, high in layout["staple_ranges"]:
            for index in range(low, high + 1):
                # The periodic reference contains theoretical capture
                # crossovers throughout the SST program.  Only targets in
                # the actual Seed/SST overlap are physical connections.
                # Copying every theoretical Seed-facing slot opened SST U
                # crossovers that were neither connected nor export sites.
                if (number, index) not in physical_sst_targets:
                    continue
                reference_record = source["stap"][index]
                current = list(row["stap"][index])
                for offset in (0, 2):
                    if 0 <= int(reference_record[offset]) < 48:
                        current[offset:offset + 2] = list(
                            reference_record[offset:offset + 2])
                row["stap"][index] = current
        capture_colors = {
            int(index): int(color)
            for index, color in source.get("stap_colors", [])
            if any(low <= int(index) <= high
                   for low, high in layout["staple_ranges"])}
        existing_colors = {
            int(index): int(color) for index, color in row.get("stap_colors", [])
            if any(low <= int(index) <= high
                   for low, high in layout["staple_ranges"])}
        existing_colors.update(capture_colors)
        row["stap_colors"] = [list(item) for item in
                              sorted(existing_colors.items())]

    # The periodic reference contains a capture crossover at every candidate
    # SST column.  Remove only cross-interface links outside the realized
    # Seed/SST overlap; internal Seed and SST staple routing is untouched.
    allowed_positions = set(map(int, layout["capture_positions"]))
    removals = []
    for number in range(48):
        row = rows[number]
        for index, record in enumerate(row["stap"]):
            for offset in (0, 2):
                partner, partner_base = map(int, record[offset:offset + 2])
                if partner >= 48 and partner_base not in allowed_positions:
                    removals.append(
                        (number, index, offset, partner, partner_base))
    for number, index, offset, partner, partner_base in removals:
        rows[number]["stap"][index][offset:offset + 2] = [-1, -1]
        partner_record = rows[partner]["stap"][partner_base]
        for partner_offset in (0, 2):
            if partner_record[partner_offset:partner_offset + 2] == \
                    [number, index]:
                partner_record[partner_offset:partner_offset + 2] = [-1, -1]
                break

    _remove_nonreciprocal_staple_links(rows)
    # The Seed side is immutable template data.  In particular, do not call
    # the former join/fill/merge post-processors here: although intended as
    # edge repair, they changed ordinary Seed staple nicks and crossovers for
    # variable SST lengths.  Only cross-interface capture edges above may be
    # installed/removed.
    return len(physical_targets)


def build_translated_capture_export_payload(payload):
    """Build an in-memory A1/B1 alternative used only for sequence export.

    The accepted design remains A0/B0 and therefore connects only Seed
    helices 16–19 and 47–44 in public SST-first numbering.  This derived copy
    swaps those bridges for their 4x4-unit translations, allowing cadnano's
    sequence engine to calculate the other four helices on each face without
    adding those bridges to the accepted structure JSON.
    """
    source_metadata = payload.get("moire_structure_metadata", {})
    source_was_sst_first = source_metadata.get(
        "helix_numbering") == "sst_first"
    result = payload_to_internal_numbering(payload)
    layout = result.get("moire_structure_metadata", {}).get(
        "variable_length_layout", structure_layout())
    _crop_reference_staples(result, capture_translation="translated")
    metadata = result.setdefault("moire_structure_metadata", {})
    metadata.update({
        "stage": "capture_export_translation",
        "capture_export_only": True,
        "capture_translation": "translated",
        "capture_site_assignments_internal":
            capture_site_assignments(layout, "translated"),
        "capture_export_site_assignments_internal":
            capture_export_site_assignments(layout),
        "capture_connection_policy": (
            "A1/B1 translated alternative generated in memory for sequence "
            "export; not part of accepted physical structure"),
        "capture_core_policy": (
            "immutable template core; length is not evaluated"),
    })
    return (payload_to_sst_first_numbering(result)
            if source_was_sst_first else result)


def _remove_nonreciprocal_staple_links(rows):
    """Open stale serialized links after mixing automatic/reference records."""
    for number, row in rows.items():
        for index, record in enumerate(row["stap"]):
            for offset in (0, 2):
                partner, partner_index = map(
                    int, record[offset:offset + 2])
                reciprocal = (
                    partner in rows and
                    0 <= partner_index < len(rows[partner]["stap"]) and
                    any(rows[partner]["stap"][partner_index][
                            partner_offset:partner_offset + 2] ==
                        [number, index]
                        for partner_offset in (0, 2)))
                if partner >= 0 and not reciprocal:
                    record[offset:offset + 2] = [-1, -1]


def _staple_components_from_rows(rows):
    nodes = {
        (number, index)
        for number, row in rows.items()
        for index, record in enumerate(row["stap"])
        if list(record) != [-1, -1, -1, -1]}
    adjacency = {node: set() for node in nodes}
    for number, index in nodes:
        record = rows[number]["stap"][index]
        for offset in (0, 2):
            other = (int(record[offset]), int(record[offset + 1]))
            if other in nodes:
                adjacency[(number, index)].add(other)
    components = []
    labels = {}
    visited = set()
    for first in sorted(nodes):
        if first in visited:
            continue
        component = set([first])
        stack = [first]
        visited.add(first)
        while stack:
            node = stack.pop()
            for other in adjacency[node]:
                if other not in visited:
                    visited.add(other)
                    component.add(other)
                    stack.append(other)
        index = len(components)
        for node in component:
            labels[node] = index
        components.append(component)
    return components, labels, adjacency


def _disconnect_staple_slot(rows, number, index, offset):
    partner, partner_index = map(
        int, rows[number]["stap"][index][offset:offset + 2])
    rows[number]["stap"][index][offset:offset + 2] = [-1, -1]
    if partner not in rows or not (
            0 <= partner_index < len(rows[partner]["stap"])):
        return
    partner_record = rows[partner]["stap"][partner_index]
    for partner_offset in (0, 2):
        if partner_record[partner_offset:partner_offset + 2] == \
                [number, index]:
            partner_record[partner_offset:partner_offset + 2] = [-1, -1]


def _staple_edge_is_reciprocal(rows, edge):
    """Return whether both records of a protected staple edge still exist."""
    first, second = tuple(edge)
    if first[0] not in rows or second[0] not in rows:
        return False
    if not (0 <= first[1] < len(rows[first[0]]["stap"]) and
            0 <= second[1] < len(rows[second[0]]["stap"])):
        return False
    first_record = rows[first[0]]["stap"][first[1]]
    second_record = rows[second[0]]["stap"][second[1]]
    return (any(first_record[offset:offset + 2] == [second[0], second[1]]
                for offset in (0, 2)) and
            any(second_record[offset:offset + 2] == [first[0], first[1]]
                for offset in (0, 2)))


def _cross_helix_staple_edges(rows, seed_only=False):
    """Return every reciprocal staple crossover as an undirected edge."""
    edges = set()
    for number, row in rows.items():
        for index, record in enumerate(row["stap"]):
            for offset in (0, 2):
                partner, partner_index = map(
                    int, record[offset:offset + 2])
                if partner < 0 or partner == number:
                    continue
                if seed_only and (number >= 48 or partner >= 48):
                    continue
                edge = tuple(sorted(((number, index),
                                     (partner, partner_index))))
                if _staple_edge_is_reciprocal(rows, edge):
                    edges.add(edge)
    return edges


def _restore_seed_template_colors(rows, layout):
    """Restore immutable Capture colours and make ordinary staples black.

    The accepted template contains light-blue/pink presentation colours on
    ordinary support staples.  They are not Capture identities and must not
    appear in a generated design.  Conversely, the 128 Capture-core colours
    and sixteen gray potential-Z2 cores are immutable identities.  Resolve
    colour by the complete template staple component so this also works when
    the template's marker happens to sit on the SST half of a bridge.
    """
    ordinary_template_colors = {0x000000, 0x60C9F6, 0xF49AE5}
    reference_rows = {
        int(row["num"]): row
        for row in _translated_frozen_reference(layout)["vstrands"]}
    reference_components, unused_labels, unused_adjacency = \
        _staple_components_from_rows(reference_rows)
    reference_markers = {
        (int(number), int(index)): int(color)
        for number, row in reference_rows.items()
        for index, color in row.get("stap_colors", [])}
    template_color_by_seed_node = {}
    for component in reference_components:
        seed_nodes = {node for node in component if node[0] < 48}
        if not seed_nodes:
            continue
        colors = {
            reference_markers[node] for node in component
            if node in reference_markers}
        if len(colors) != 1:
            raise RuntimeError(
                "冻结Seed范本的每条staple必须且只能有一个颜色标记；"
                "当前组件有%d个。" % len(colors))
        template_color = next(iter(colors))
        final_color = (0x000000 if template_color in
                       ordinary_template_colors else template_color)
        for node in seed_nodes:
            template_color_by_seed_node[node] = final_color

    components, unused_labels, unused_adjacency = \
        _staple_components_from_rows(rows)
    marker_colors = {
        (int(number), int(index)): int(color)
        for number, row in rows.items()
        for index, color in row.get("stap_colors", [])}
    for component in components:
        seed_nodes = {node for node in component if node[0] < 48}
        if not seed_nodes:
            # SST-output staples that are not Capture extensions are ordinary
            # staples too.  Keep them black in the design JSON.
            final_color = 0x000000
        else:
            mapped_colors = {
                template_color_by_seed_node[node]
                for node in seed_nodes if node in template_color_by_seed_node}
            if len(mapped_colors) != 1:
                raise RuntimeError(
                    "生成的Seed staple无法唯一映射到冻结范本颜色；"
                    "检测到%d种候选颜色。" % len(mapped_colors))
            final_color = next(iter(mapped_colors))
        existing = sorted(node for node in component if node in marker_colors)
        marker = existing[0] if existing else min(component)
        for node in component:
            marker_colors.pop(node, None)
        marker_colors[marker] = final_color

    for number, row in rows.items():
        row["stap_colors"] = [
            [index, color]
            for (helix, index), color in sorted(marker_colors.items())
            if helix == number]


def _normalize_capture_component_colors(rows, layout):
    """Give every immutable Capture column one independent colour.

    Pair grouping remains available to the routing rules, but it is not a
    display identity.  Recover each Capture core's real-space column from
    the frozen 2L template, then recolour every generated component mapped
    to that core.  This covers physical bridges, translated candidates,
    Kagome-hole candidates and both unconnected Z2 reserve columns.  Chain
    topology and endpoints are never changed.
    """
    reference_rows = {
        int(row["num"]): row
        for row in _translated_frozen_reference(layout)["vstrands"]}
    reference_components, unused_labels, unused_adjacency = \
        _staple_components_from_rows(reference_rows)
    reference_markers = {
        (int(number), int(index)): int(color)
        for number, row in reference_rows.items()
        for index, color in row.get("stap_colors", [])}
    ordinary_reference_colors = {0x000000, 0x60C9F6, 0xF49AE5}
    shift = _seed_canvas_shift(layout)
    template_column_by_seed_node = {}
    for component in reference_components:
        seed_nodes = {node for node in component if node[0] < 48}
        if not seed_nodes:
            continue
        colors = {
            reference_markers[node] for node in component
            if node in reference_markers}
        if len(colors) != 1:
            raise RuntimeError(
                "冻结Seed范本的每条staple必须且只能有一个颜色标记；"
                "当前组件有%d个。" % len(colors))
        reference_color = next(iter(colors))
        if reference_color in ordinary_reference_colors:
            continue
        if reference_color == 0x999999:
            logical_candidates = (184, 200)
        else:
            logical_column = CAPTURE_REFERENCE_COLUMN_BY_COLOR.get(
                reference_color)
            if logical_column is None:
                raise RuntimeError(
                    "冻结Seed范本包含无法映射到capture列的颜色#%06x。" %
                    reference_color)
            logical_candidates = (logical_column,)
        matched_columns = []
        for logical_column in logical_candidates:
            column = int(logical_column) + shift
            unit_index = (int(logical_column) -
                          CAPTURE_DIRECT_POSITIONS[0]) // 16
            phase = "B" if unit_index % 2 == 0 else "A"
            capture_helices = {
                int(seed_helix)
                for cycle in (phase + "0", phase + "1")
                for unused_sst, seed_helix in
                CAPTURE_PHASE_MAPPINGS[cycle]}
            if any((helix, column) in component
                   for helix in capture_helices):
                matched_columns.append(column)
        if len(matched_columns) != 1:
            raise RuntimeError(
                "冻结Seed capture核心无法唯一映射到真实列：%s。" %
                matched_columns)
        column = matched_columns[0]
        for node in seed_nodes:
            template_column_by_seed_node[node] = column

    components, unused_labels, unused_adjacency = \
        _staple_components_from_rows(rows)
    marker_colors = {
        (int(number), int(index)): int(color)
        for number, row in rows.items()
        for index, color in row.get("stap_colors", [])}
    for component in components:
        seed_nodes = {node for node in component if node[0] < 48}
        columns = {
            template_column_by_seed_node[node]
            for node in seed_nodes if node in template_column_by_seed_node}
        if not columns:
            continue
        if len(columns) != 1:
            raise RuntimeError(
                "一条Seed capture核心跨越了多个真实列：%s。" %
                sorted(columns))
        column = next(iter(columns))
        color = capture_column_color(column, layout)
        if color is None:
            raise RuntimeError(
                "Seed capture列base %d没有独立显示颜色。" % column)
        for node in component:
            marker_colors.pop(node, None)
        marker_colors[min(component)] = color
    for number, row in rows.items():
        row["stap_colors"] = [
            [index, color]
            for (helix, index), color in sorted(marker_colors.items())
            if helix == number]






def _capture_pair_index(base_index, layout):
    """Return the zero-based cooperative pair for a capture oligo."""
    return capture_pair_index(base_index, layout)




def scaffold(output, input_json=None):
    sst = (json.loads(Path(input_json).read_text(encoding="utf-8"))
           if input_json else None)
    layout = ((sst or {}).get("moire_structure_metadata", {}).get(
        "variable_length_layout") or structure_layout())
    layout["seed_z1_requested_bp"] = 128
    layout["seed_z3_requested_bp"] = 128
    layout["seed_z1_actual_bp"] = 128
    layout["seed_z3_actual_bp"] = 128
    layout["seed_z1_edge_growth_bp"] = 0
    layout["seed_z3_edge_growth_bp"] = 0
    layout["seed_length_adjustment_enabled"] = False
    layout["seed_geometry_policy"] = "immutable_2L_reference"
    square_geometry = layout.get("square_centered_geometry", {})
    # Kagome and mixed SST files historically stored the centred Seed shift
    # only inside ``square_centered_geometry``.  Recover it before creating
    # the frozen Seed so all three lattice combinations share one coordinate
    # frame through SST, scaffold and finalization.
    layout["coordinate_shift_bp"] = int(layout.get(
        "coordinate_shift_bp",
        square_geometry.get("coordinate_shift_bp", 0)) or 0)
    canvas_shift = _seed_canvas_shift(layout)
    seed_ranges = _seed_capture_support_ranges(layout)
    overlap = fixed_seed_overlap_layout(
        layout.get("layer_ranges", [
            [48 + canvas_shift, 175 + canvas_shift],
            [208 + canvas_shift, 335 + canvas_shift]]),
        (layout.get("theoretical_capture_positions_by_layer")
         if layout.get("lattice_type") in
         ("kagome", "square_kagome") else None),
        lattice_type=layout.get("lattice_type", "square"),
        seed_layer_ranges=seed_ranges,
        seed_capture_positions_by_layer=(
            layout.get("seed_capture_positions_by_layer")
            if layout.get("lattice_type") in
            ("kagome", "square_kagome") else
            square_geometry.get("seed_capture_positions_by_layer")))
    reported_ranges = square_geometry.get(
        "optimized_overlap_ranges", layout.get(
            "overlap_ranges", overlap["overlap_ranges"]))
    reported_overlap = square_geometry.get(
        "optimized_seed_overlap_bp", layout.get(
            "seed_sst_overlap_bp", overlap["seed_sst_overlap_bp"]))
    layout["overlap_ranges"] = copy.deepcopy(reported_ranges)
    layout["seed_sst_overlap_bp"] = list(map(int, reported_overlap))
    layout["capture_positions_by_layer"] = overlap[
        "capture_positions_by_layer"]
    layout["capture_positions"] = overlap["capture_positions"]
    layout["capture_columns_by_layer"] = overlap[
        "capture_columns_by_layer"]
    layout["capture_sites_per_column_by_layer"] = overlap[
        "capture_sites_per_column_by_layer"]
    layout["capture_sites_by_layer"] = overlap[
        "capture_sites_by_layer"]
    layout["capture_pair_equivalents_by_layer"] = overlap[
        "capture_pair_equivalents_by_layer"]
    layout["capture_pair_equivalents"] = overlap[
        "capture_pair_equivalents"]
    layout["pair_count_by_layer"] = overlap["pair_count_by_layer"]
    layout["pair_count"] = overlap["pair_count"]
    if not overlap["capture_pairs_valid"]:
        raise ValueError(
            "每层至少需要2列capture pair；当前为%s。请减小SST spacing。" %
            "/".join(map(str, overlap["pair_count_by_layer"])))
    positions_by_layer = layout.get(
        "capture_positions_by_layer",
        (sst or {}).get("moire_structure_metadata", {}).get(
            "kagome_theoretical_capture_positions_by_layer", [[], []]))
    layout.setdefault("capture_positions_by_layer", positions_by_layer)
    layout.setdefault("capture_positions", [
        value for layer in positions_by_layer for value in layer])
    layout.setdefault("capture_column_count",
                      len(layout.get("capture_positions", [])))
    routing_layout = _seed_routing_layout(layout)
    seed, loops = _reference_seed_scaffold_payload(routing_layout)
    indel_lengths = _apply_frozen_seed_z2_indels(seed, routing_layout)
    for key in (
            "mean_indel_per_helix",
            "mean_indel_per_helix_requested",
            "mean_indel_per_helix_actual",
            "actual_z2_spacing_bp",
            "seed_z2_range",
            "seed_z2_indel_range",
            "seed_z2_indel_placements",
            "seed_z2_domain_bp",
            "maximum_seed_indels_per_domain",
            "maximum_seed_deletions_per_domain",
            "minimum_seed_deletion_per_helix",
            "seed_scaffold_lengths_after_indel",
            "seed_z2_indel_distribution",
            "seed_indel_policy"):
        if key in routing_layout:
            layout[key] = copy.deepcopy(routing_layout[key])
    seed["name"] = Path(output).name
    payload = _combine_seed_sst(seed, "scaffold_review", sst)
    # ``_seed_routing_layout`` returns an isolated Square routing view for a
    # Kagome SST.  The combined payload contains the authoritative values
    # copied back from the frozen Seed; report that layout rather than the
    # pre-routing Kagome input dictionary.
    layout = payload["moire_structure_metadata"]["variable_length_layout"]
    for key in (
            "mean_indel_per_helix",
            "mean_indel_per_helix_requested",
            "mean_indel_per_helix_actual",
            "actual_z2_spacing_bp",
            "seed_z2_range",
            "seed_z2_indel_range",
            "seed_z2_indel_placements",
            "seed_z2_domain_bp",
            "maximum_seed_indels_per_domain",
            "maximum_seed_deletions_per_domain",
            "minimum_seed_deletion_per_helix",
            "seed_scaffold_lengths_after_indel",
            "seed_z2_indel_distribution",
            "seed_indel_policy"):
        if key in routing_layout:
            layout[key] = copy.deepcopy(routing_layout[key])
    payload = payload_to_sst_first_numbering(payload)
    Path(output).write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                            encoding="utf-8")
    return {
        "output": str(Path(output).resolve()),
        "stage": "scaffold_review",
        "seed_scaffold_components": len(loops),
        "seed_scaffold_single_nicks": len(loops),
        "seed_scaffold_lengths": sorted(indel_lengths),
        "seed_scaffold_partition": seed["moire_edge_metadata"]["partition"],
        "sst_first_base": layout["layer_ranges"][0][0],
        "sst_ranges": layout["layer_ranges"],
        "seed_requested_lengths": [
            layout["seed_z1_requested_bp"],
            layout["seed_z3_requested_bp"]],
        "seed_routing_lengths": [
            layout["seed_z1_actual_bp"], layout["seed_z3_actual_bp"]],
        "seed_support_ranges": layout["seed_layer_ranges"],
        "overlap_ranges": layout["overlap_ranges"],
        "capture_column_count": layout["capture_column_count"],
        "capture_columns_by_layer": [
            len(item) for item in layout["capture_positions_by_layer"]],
        "maximum_edge_stagger_bp": seed["moire_edge_metadata"].get(
            "maximum_edge_stagger_bp", 8),
        "variable_length_layout": layout,
    }


def capacity(input_json=None):
    """Report capacity from the same frozen Seed routing used for output."""
    sst = (json.loads(Path(input_json).read_text(encoding="utf-8"))
           if input_json else None)
    layout = ((sst or {}).get("moire_structure_metadata", {}).get(
        "variable_length_layout") or structure_layout())
    layout["seed_z1_requested_bp"] = 128
    layout["seed_z3_requested_bp"] = 128
    routing_layout = _seed_routing_layout(layout)
    seed, band_lengths = _reference_seed_scaffold_payload(routing_layout)
    edge = seed["moire_edge_metadata"]
    total_nt = int(edge["total_seed_scaffold_nt"])
    plan = scaffold_capacity_plan(total_nt)
    insertion_headroom = sum(
        max(0, int(plan["per_scaffold_capacity_nt"]) - int(length))
        for length in band_lengths)
    return {
        "stage": "capacity_precheck",
        "seed_scaffold_total_nt": total_nt,
        "seed_scaffold_count": plan["count"],
        "per_scaffold_capacity_nt": plan["per_scaffold_capacity_nt"],
        "total_capacity_nt": plan["total_capacity_nt"],
        "planned_balanced_lengths": list(band_lengths),
        "seed_insertion_headroom_nt": insertion_headroom,
        "maximum_mean_insertion_from_scaffold_capacity":
            insertion_headroom / 48.0,
        "maximum_edge_stagger_bp": int(edge["maximum_edge_stagger_bp"]),
        "seed_routing_reference": edge["reference"],
    }


def _remove_orphan_legacy_color_anchors(payload):
    """Remove colour markers that caDNAno cannot associate with a strand.

    Short-spacing SST routing uses auxiliary virtual helices whose topology
    contains only the detoured interval.  Final colour normalization may
    nevertheless copy markers for empty positions on those rows.  caDNAno's
    legacy decoder dereferences the strand for every marker without checking
    for ``None``, so a single empty anchor makes an otherwise valid JSON
    impossible to open.  Filtering these markers is serialization-only: no
    scaffold, staple, crossover, nick, sequence or valid colour is changed.
    """
    removed = {"stap_colors": 0, "scaf_colors": 0}
    empty = [-1, -1, -1, -1]
    for row in payload.get("vstrands", []):
        for topology_name, colors_name in (
                ("stap", "stap_colors"), ("scaf", "scaf_colors")):
            topology = row.get(topology_name, [])
            valid = []
            for item in row.get(colors_name, []):
                try:
                    index = int(item[0])
                    int(item[1])
                except (IndexError, TypeError, ValueError):
                    removed[colors_name] += 1
                    continue
                if index < 0 or index >= len(topology) or \
                        topology[index] == empty:
                    removed[colors_name] += 1
                    continue
                valid.append(list(item))
            if colors_name in row or valid:
                row[colors_name] = valid
    return removed




def _finalize_frozen_square_reference(source, layout, output):
    """Finalize canonical 128/32/128 from the accepted two-layer template.

    The reference is already the reviewed result of the normal caDNAno
    AutoCS/Autobreak workflow.  This fast path never reruns those commands:
    it freezes scaffold, ordinary staple and capture-core topology.  No Seed
    edge trimming, length repair, seam synthesis or AutoCS rerun occurs.
    """
    exact_reference = (
        (int(layout.get("z1_bp", 128)), int(layout.get("z2_bp", 32)),
         int(layout.get("z3_bp", 128))) == (128, 32, 128) and
        int(layout.get("coordinate_shift_bp", 0)) == 0)
    # Freeze only the Seed half of the accepted reference.  The SST half is
    # the already validated, length/spacing-specific source produced by the
    # lattice-specific SST builder and must never be replaced by the
    # canonical 128/32/128 SST rows during finalization.
    frozen = payload_to_internal_numbering(
        _translated_frozen_reference(layout))
    payload = (copy.deepcopy(frozen) if exact_reference
               else copy.deepcopy(source))
    source_rows = {int(row["num"]): row for row in source["vstrands"]}
    reference_rows = {int(row["num"]): row for row in frozen["vstrands"]}
    for number in range(48):
        source_scaf = source_rows[number]["scaf"]
        reference_scaf = reference_rows[number]["scaf"]
        length = max(len(source_scaf), len(reference_scaf))
        empty = [-1, -1, -1, -1]
        normalized_source = list(source_scaf) + [empty] * (
            length - len(source_scaf))
        normalized_reference = list(reference_scaf) + [empty] * (
            length - len(reference_scaf))
        if normalized_source != normalized_reference:
            raise RuntimeError(
                "已接受Seed scaffold与冻结2L范本不一致；"
                "请先重新生成或在专家模式确认，不得静默覆盖。")

    payload_rows = {int(row["num"]): row for row in payload["vstrands"]}
    array_length = int(layout.get("array_length", payload.get(
        "num_bases", SST_ARRAY_LENGTH)))
    # caDNAno requires every virtual helix array to have exactly the declared
    # document length.  The frozen template itself is 384 bases; pad both
    # Seed and SST rows when the project canvas is larger without changing
    # any occupied record.
    for row in payload_rows.values():
        _resize_row(row, array_length)
    for number in range(48):
        row = copy.deepcopy(reference_rows[number])
        _resize_row(row, array_length)
        payload_rows[number].clear()
        payload_rows[number].update(row)
    payload["vstrands"] = [
        payload_rows[number] for number in sorted(payload_rows)]
    payload["num_bases"] = array_length

    if not exact_reference:
        # Install only the template-defined capture endpoints on top of the
        # current SST topology.  This operation does not inspect, shorten,
        # split, rebreak or reject a capture core by length.
        _crop_reference_staples(payload, capture_translation="origin")
    # Design colours are a separate concern from sequence-text rich colours:
    # all ordinary staples are black, while Capture/Z2 cores retain the exact
    # accepted-template identity colour until the current Capture columns are
    # normalized below.
    _restore_seed_template_colors(payload_rows, layout)
    _normalize_capture_component_colors(payload_rows, layout)

    # Reproduce the same deterministic Z2-only indels shown in scaffold
    # review.  No scaffold/staple edge or color is changed.
    actual_scaffold_lengths = _apply_frozen_seed_z2_indels(payload, layout)

    metadata = copy.deepcopy(source.get("moire_structure_metadata", {}))
    layout = copy.deepcopy(layout)
    layout.update({
        "seed_routing_is_frozen_reference": True,
        "seed_routing_reference": "Square_Seed_2L_newtemplate.json",
        "seed_edge_stagger_limit_used_bp": 14,
        "seed_edge_selected_scaffold_count": 2,
        "capture_helices_define_maximum_actual_length": False,
        "seed_z1_requested_bp": 128,
        "seed_z3_requested_bp": 128,
        "seed_z1_actual_bp": int(layout.get("seed_z1_actual_bp", 128)),
        "seed_z3_actual_bp": int(layout.get("seed_z3_actual_bp", 128)),
        "seed_z1_edge_growth_bp": 0,
        "seed_z3_edge_growth_bp": 0,
        "seed_length_adjustment_enabled": False,
        "seed_geometry_policy": "immutable_2L_reference",
        "seed_staple_support_ranges": copy.deepcopy(
            layout.get("seed_layer_ranges", [[48, 175], [208, 335]])),
        # The two support regions plus Z2 form one continuous physical Seed
        # duplex band.  Scaffold outside this range is deliberately retained
        # from the immutable routing reference as scaffold-only support.
        "seed_staple_physical_range": [
            48 + _seed_canvas_shift(layout),
            335 + _seed_canvas_shift(layout)],
        # The immutable reference staggers several helix ends by one 8-nt
        # domain.  Bases 64..319 are therefore the common all-helix duplex
        # region whose scaffold coverage must be complete; 48..63 and
        # 320..335 remain edge-repair windows governed by the frozen
        # reference and capture-nick exceptions.
        "seed_staple_required_coverage_range": [
            64 + _seed_canvas_shift(layout),
            319 + _seed_canvas_shift(layout)],
    })
    metadata.update({
        "stage": "structure_complete",
        "lattice_type": "square",
        "seed_cross_section": "S8-R4x4C (48 helices)",
        "seed_cross_section_preset": "s8_r4x4",
        "sst_layers": 2,
        "sst_scaffold_ranges": copy.deepcopy(layout.get(
            "scaffold_ranges", layout.get(
                "sst_scaffold_ranges", layout["layer_ranges"]))),
        "sst_staple_ranges": copy.deepcopy(layout["staple_ranges"]),
        "variable_length_layout": layout,
        "seed_scaffold_count": 2,
        "seed_scaffold_lengths_after_indel": actual_scaffold_lengths,
        "accepted_scaffold_source": str(Path(output).resolve()),
        "staple_generation": (
            "frozen accepted 2L topology; no Seed trim/growth/seam, no "
            "AutoCS/Autobreak rerun, no capture-core length processing"),
        "seed_staple_coverage": "immutable accepted 2L template",
        "capture_pair_color_policy": (
            "all ordinary staples black; Capture/Z2 cores retain exact "
            "Square_Seed_2L_newtemplate.json colours"),
        "capture_pair_colors": [
            CAPTURE_PAIR_COLORS[index % len(CAPTURE_PAIR_COLORS)]
            for index in range(layout["pair_count"])],
        "capture_column_colors": list(CAPTURE_COLUMN_COLORS),
    })
    payload["name"] = Path(output).name
    payload["moire_structure_metadata"] = metadata
    payload["moire_structure_metadata"].update({
        "capture_anchor_positions": list(layout["capture_positions"]),
        "capture_positions_by_layer": copy.deepcopy(
            layout["capture_positions_by_layer"]),
        "capture_site_assignments_internal": copy.deepcopy(
            layout["capture_site_assignments"]),
        "capture_export_site_assignments_internal": copy.deepcopy(
            layout["capture_export_site_assignments"]),
        "autobreak": {
            "supported": True,
            "already_applied": True,
            "source": "Square_Seed_2L_newtemplate.json",
            "rerun": False,
        },
        "capture_core_policy": (
            "immutable template core; length is not evaluated"),
    })
    payload = payload_to_sst_first_numbering(payload)
    removed_color_anchors = _remove_orphan_legacy_color_anchors(payload)
    payload["moire_structure_metadata"][
        "legacy_color_anchor_sanitization"] = removed_color_anchors
    Path(output).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return {
        "output": str(Path(output).resolve()),
        "stage": "structure_complete",
        "sst_first_base": 48,
        "capture_bridge_positions": list(layout["capture_positions"]),
        "legacy_color_anchor_sanitization": removed_color_anchors,
        "variable_length_layout": layout,
        "autobreak": payload["moire_structure_metadata"]["autobreak"],
    }




def _frozen_kagome_capture_reference(seed_source=None):
    """Build the canonical two-layer Kagome final topology.

    Seed scaffold and ordinary Seed staples always come from the accepted
    two-layer Square reference.  The historical three-layer Kagome file is
    consulted only for (a) the Kagome SST rows and (b) Seed--SST capture
    links in its first two layers.  This separation is intentional: the 3L
    file is a capture catalogue, never a Seed routing template.
    """
    resource = HERE.with_name("resources") / \
        "Kagome_Seed_Ka-seed-pore_3L.json"
    if not resource.is_file():
        raise FileNotFoundError(resource)
    payload = payload_to_internal_numbering(copy.deepcopy(_reference()))
    catalogue = payload_to_internal_numbering(
        json.loads(resource.read_text(encoding="utf-8")))
    rows = {int(row["num"]): row for row in payload["vstrands"]}
    catalogue_rows = {
        int(row["num"]): row for row in catalogue["vstrands"]}

    # The Kagome catalogue owns only the SST/capture topology.  Seed scaffold
    # routing comes from the already reviewed 2L reference transformation
    # (including its +11-bp edge compression); never silently restore the
    # longer physical edges carried by the historical catalogue file.
    if seed_source is not None:
        accepted_rows = {
            int(row["num"]): row for row in seed_source["vstrands"]}
        for number in range(48):
            rows[number]["scaf"] = copy.deepcopy(
                accepted_rows[number]["scaf"])
            rows[number]["loop"] = copy.deepcopy(
                accepted_rows[number].get("loop", rows[number]["loop"]))
            rows[number]["skip"] = copy.deepcopy(
                accepted_rows[number].get("skip", rows[number]["skip"]))

    # Replace only the first-two-layer SST topology.  Everything at or beyond
    # base 344 belongs to the third layer in the catalogue and is discarded.
    for number in range(48, 64):
        row = rows[number]
        source_row = catalogue_rows[number]
        for field in ("scaf", "stap"):
            for index in range(min(344, len(row[field]))):
                row[field][index] = copy.deepcopy(source_row[field][index])
            for index in range(344, len(row[field])):
                row[field][index] = [-1, -1, -1, -1]
        row["stap_colors"] = [
            [int(index), int(color)]
            for index, color in source_row.get("stap_colors", [])
            if int(index) < 344 and
            row["stap"][int(index)] != [-1, -1, -1, -1]]

    # Remove Square Seed--SST capture links, leaving ordinary accepted 2L
    # Seed staples untouched, then install only the Kagome capture catalogue.
    for number in range(48):
        for index, record in enumerate(rows[number]["stap"]):
            for offset in (0, 2):
                if int(record[offset]) >= 48:
                    _disconnect_staple_slot(
                        rows, number, index, offset)
    for number in range(48):
        source_row = catalogue_rows[number]
        for index in range(min(344, len(source_row["stap"]))):
            record = source_row["stap"][index]
            for offset in (0, 2):
                partner, partner_index = map(
                    int, record[offset:offset + 2])
                if partner < 48:
                    continue
                rows[number]["stap"][index][offset:offset + 2] = [
                    partner, partner_index]
                partner_record = catalogue_rows[partner]["stap"][
                    partner_index]
                reciprocal_offset = (
                    0 if partner_record[:2] == [number, index] else 2)
                rows[partner]["stap"][partner_index][
                    reciprocal_offset:reciprocal_offset + 2] = [
                        number, index]
    _remove_nonreciprocal_staple_links(rows)

    # The catalogue's capture cores are authoritative.  Do not shorten,
    # split, rebreak or otherwise normalize them according to a length rule.
    return payload, []


def _finalize_frozen_kagome_reference(source, layout, output):
    """Finalize Kagome from its frozen two-layer capture catalogue.

    Seed routing, ordinary staples and Kagome capture cores are accepted
    template data.  No Seed edge or capture-length optimizer is permitted.
    """
    # Preserve the length/spacing-specific Kagome SST, including any real
    # h64--79 sequence detour.  Only the immutable 2L Seed rows are imported
    # from the accepted reference; the historical 3L Kagome file remains a
    # phase/capture catalogue and must never replace a generated SST.
    payload = copy.deepcopy(source)
    frozen = payload_to_internal_numbering(
        _translated_frozen_reference(layout))
    source_rows = {int(row["num"]): row for row in source["vstrands"]}
    frozen_rows = {int(row["num"]): row for row in frozen["vstrands"]}
    for number in range(48):
        source_scaf = source_rows[number]["scaf"]
        frozen_scaf = frozen_rows[number]["scaf"]
        length = max(len(source_scaf), len(frozen_scaf))
        empty = [-1, -1, -1, -1]
        if (list(source_scaf) + [empty] * (length - len(source_scaf))) != \
                (list(frozen_scaf) + [empty] * (length - len(frozen_scaf))):
            raise RuntimeError(
                "已接受Seed scaffold与冻结2L范本不一致；"
                "禁止为Kagome SST静默覆盖Seed。")
    payload_rows = {int(row["num"]): row for row in payload["vstrands"]}
    for number in range(48):
        payload_rows[number] = copy.deepcopy(frozen_rows[number])
    payload["vstrands"] = [
        payload_rows[number] for number in sorted(payload_rows)]
    rows = {int(row["num"]): row for row in payload["vstrands"]}
    reference_crossovers = _cross_helix_staple_edges(
        rows, seed_only=True)
    _crop_reference_staples(payload, capture_translation="origin")
    _normalize_capture_component_colors(rows, layout)
    half_positions = []
    layout = copy.deepcopy(layout)
    rows = {int(row["num"]): row for row in payload["vstrands"]}
    array_length = int(layout.get(
        "array_length", source.get("num_bases", payload.get(
            "num_bases", SST_ARRAY_LENGTH))))
    # The frozen Kagome catalogue is a 384-base topology reference, while
    # the accepted project canvas can be longer.  caDNAno requires scaf,
    # stap, loop and skip to have one common length on every virtual helix,
    # equal to the top-level ``num_bases`` value.  Padding empty records is a
    # serialization repair only: no occupied SST/capture edge or phase is
    # changed here.
    for row in rows.values():
        _resize_row(row, array_length)
    payload["num_bases"] = array_length
    # The Kagome catalogue helper may carry the review file's loop/skip
    # arrays.  Reset only indel arrays to the immutable template before the
    # deterministic allocation below, otherwise finalization would apply the
    # same requested mean twice.
    reference_rows = {
        int(row["num"]): row
        for row in _translated_frozen_reference(layout)["vstrands"]}
    for number in range(48):
        _resize_row(reference_rows[number], array_length)
        rows[number]["loop"] = copy.deepcopy(
            reference_rows[number]["loop"])
        rows[number]["skip"] = copy.deepcopy(
            reference_rows[number]["skip"])
    actual_scaffold_lengths = _apply_frozen_seed_z2_indels(payload, layout)
    layout.update({
        # Kagome changes only the SST/capture catalogue.  Its Seed scaffold
        # is still copied byte-for-byte from the same accepted frozen 2L
        # Square routing and must be validated under that routing contract.
        "seed_routing_is_frozen_reference": True,
        "seed_routing_reference": "Square_Seed_2L_newtemplate.json",
        "kagome_capture_topology_reference": (
            "Kagome_Seed_Ka-seed-pore_3L.json:first-two-layers/capture-only"),
        "seed_edge_stagger_limit_used_bp": 14,
        "seed_edge_selected_scaffold_count": 2,
        "capture_helices_define_maximum_actual_length": False,
        "seed_z1_requested_bp": 128,
        "seed_z3_requested_bp": 128,
        "seed_z1_actual_bp": int(layout.get("seed_z1_actual_bp", 128)),
        "seed_z3_actual_bp": int(layout.get("seed_z3_actual_bp", 128)),
        "seed_z1_edge_growth_bp": 0,
        "seed_z3_edge_growth_bp": 0,
        "seed_length_adjustment_enabled": False,
        "seed_geometry_policy": "immutable_2L_reference",
        "seed_staple_physical_range": [
            48 + _seed_canvas_shift(layout),
            335 + _seed_canvas_shift(layout)],
        "seed_staple_required_coverage_range": [
            64 + _seed_canvas_shift(layout),
            319 + _seed_canvas_shift(layout)],
    })
    metadata = copy.deepcopy(source.get("moire_structure_metadata", {}))
    seed_scaffold_count = len(_scaffold_components({
        number: rows[number] for number in range(48)}))
    final_lattice = str(layout.get("lattice_type", "kagome"))
    metadata.update({
        "stage": "structure_complete",
        # The same finalizer installs the Kagome layer of a mixed design,
        # but the file remains a true layer-specific Square--Kagome design.
        # Relabelling it as all-Kagome makes reopened files lose the Square
        # capture catalogue and causes structure exports to choose the wrong
        # ideal cross-section.
        "lattice_type": final_lattice,
        "lattice_by_layer": copy.deepcopy(layout.get(
            "lattice_by_layer", ["kagome", "kagome"])),
        "seed_cross_section": "S8-R4x4C (48 helices)",
        "seed_cross_section_preset": "s8_r4x4",
        "sst_layers": 2,
        "sst_scaffold_ranges": copy.deepcopy(layout.get(
            "scaffold_ranges", layout.get(
                "sst_scaffold_ranges", layout["layer_ranges"]))),
        "sst_staple_ranges": copy.deepcopy(layout["staple_ranges"]),
        "variable_length_layout": layout,
        "seed_scaffold_count": seed_scaffold_count,
        "seed_scaffold_lengths_after_indel": actual_scaffold_lengths,
        "seed_routing_source": "Square_Seed_2L_newtemplate.json",
        "kagome_capture_topology_source": (
            "Kagome_Seed_Ka-seed-pore_3L.json first two layers only"),
        "staple_generation": (
            "accepted 2L ordinary Seed topology plus Kagome-specific "
            "capture projection; no legacy Moiré post-optimizer"),
        "capture_half_crossover_fallback_positions": half_positions,
    })
    payload["name"] = Path(output).name
    payload["moire_structure_metadata"] = metadata
    payload["moire_structure_metadata"].update({
        "autobreak": {
            "supported": True,
            "already_applied": True,
            "source": "accepted 2L + Kagome capture catalogue",
            "rerun": False,
        },
        "capture_core_policy": (
            "immutable Kagome catalogue core; length is not evaluated"),
    })
    if _cross_helix_staple_edges(rows, seed_only=True) != reference_crossovers:
        raise RuntimeError("Kagome固定capture投影意外改动Seed普通staple。")
    payload = payload_to_sst_first_numbering(payload)
    removed_color_anchors = _remove_orphan_legacy_color_anchors(payload)
    payload["moire_structure_metadata"][
        "legacy_color_anchor_sanitization"] = removed_color_anchors
    Path(output).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return {
        "output": str(Path(output).resolve()),
        "stage": "structure_complete",
        "sst_first_base": 48,
        "capture_bridge_positions": list(layout["capture_positions"]),
        "legacy_color_anchor_sanitization": removed_color_anchors,
        "capture_half_crossover_fallback_positions": half_positions,
        "variable_length_layout": layout,
        "autobreak": payload["moire_structure_metadata"]["autobreak"],
    }


def finalize(input_json, output):
    """Finalize only through the frozen two-layer Seed template paths.

    SST generation remains lattice-specific.  Seed shrink/growth, edge seam
    synthesis, capture-core length repair and fallback AutoCS/Autobreak are
    intentionally absent from this dispatcher.
    """
    source = payload_to_internal_numbering(
        json.loads(Path(input_json).read_text(encoding="utf-8")))
    layout = source.get("moire_structure_metadata", {}).get(
        "variable_length_layout", structure_layout())
    layout["seed_z1_requested_bp"] = 128
    layout["seed_z3_requested_bp"] = 128
    layout.setdefault("seed_z1_actual_bp", 128)
    layout.setdefault("seed_z3_actual_bp", 128)
    layout["seed_z1_edge_growth_bp"] = 0
    layout["seed_z3_edge_growth_bp"] = 0
    layout["seed_length_adjustment_enabled"] = False
    layout["seed_geometry_policy"] = "immutable_2L_reference"
    layout["seed_routing_is_frozen_reference"] = True
    layout["seed_routing_reference"] = "Square_Seed_2L_newtemplate.json"
    lattice_type = layout.get("lattice_type", "square")
    if lattice_type == "square":
        return _finalize_frozen_square_reference(source, layout, output)
    if lattice_type == "kagome":
        return _finalize_frozen_kagome_reference(source, layout, output)
    if lattice_type == "square_kagome":
        return _finalize_frozen_kagome_reference(source, layout, output)
    raise ValueError("不支持的SST点阵类型：%s。" % lattice_type)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("capacity", "scaffold", "finalize"))
    parser.add_argument("output")
    parser.add_argument("--input")
    args = parser.parse_args()
    Path(args.output).expanduser().resolve().parent.mkdir(
        parents=True, exist_ok=True)
    if args.mode == "capacity":
        result = capacity(args.input)
    elif args.mode == "scaffold":
        result = scaffold(args.output, args.input)
    else:
        if not args.input:
            parser.error("--input is required for finalize")
        result = finalize(args.input, args.output)
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()

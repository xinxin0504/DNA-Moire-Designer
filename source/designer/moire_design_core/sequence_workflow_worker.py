#!/usr/bin/env python3
"""Headless sequence assignment, workbook and final-package worker."""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from moire_design_core import structure_worker as runtime
from moire_design_core.sequence_export_worker import (
    _add_capture_staple_columns,
    _apply_capture_sequence_rich_text,
    _capture_sequence_manifest,
    _capture_staple_row,
    _render_svg,
    _rows_for_variant,
    _sort_capture_manifest_by_template_color,
    _staple_base_sequences,
    _write_json,
)
from moire_design_core.structure import (
    CAPTURE_DIRECT_POSITIONS,
    CAPTURE_PHASE_MAPPINGS,
    CAPTURE_REFERENCE_COLUMN_BY_COLOR,
    SEED_CAPTURE_REFERENCE,
    _load_json,
    build_complete_sst_only_payload,
    capture_column_color,
    capture_column_index,
    capture_pair_index,
    payload_to_internal_numbering,
    payload_to_sst_first_numbering,
)
from cadnano2.model.enum import LatticeType
from cadnano2.model.io.legacydecoder import import_legacy_dict
from cadnano2.model.io.sequencexlsx import (
    HEADERS,
    _write_workbook,
    read_sequence_template,
    write_sequence_template,
)
from cadnano2.model.io.oxdnaexport import export_structure_bundle
from cadnano2.model.io import oxdnaexport as _oxdna
from cadnano2.data.dnasequences import sequences as CADNANO_SCAFFOLDS
from moire_designer.i18n import (localize_xlsx, set_language)


BASES = set("ACGT")
SCAFFOLD_COLORS = (
    "#1769aa", "#d1495b", "#2a9d8f", "#7b61b8",
    "#e17c05", "#348aa7", "#8f5d2f", "#5c946e",
)
STANDARD_SCAFFOLD_KEYS = (
    ("CS3L", "CS3L_7559"),
    ("CS4", "CS4_7557"),
    ("P7560", "p7560"),
)
MULTI_SCAFFOLD_NAMES = ("CS3L", "CS4", "P7560")
SCAFFOLD_NAME_ALIASES = {
    # Compatibility for projects saved before the display-name correction.
    "CS3": "CS3L",
    "CS4-L": "CS4",
}
_SEQUENCE_POSITION_RE = re.compile(
    r'^(?:[HS]\s*:\s*)?(-?\d+)\s*\[\s*(-?\d+)\s*\]$',
    re.IGNORECASE)

# The frozen Seed uses two non-black colours for ordinary support staples.
# They are visual categories in the reference design, not Capture products.
# The remaining non-black/non-gray components are the 128 immutable Capture
# cores; the sixteen gray components are the potential Z2 cores.
_SEED_NORMAL_TEMPLATE_COLORS = {
    "#000000", "#60c9f6", "#f49ae5",
}
_SEED_POTENTIAL_Z2_COLOR = "#999999"


def canonical_scaffold_name(name):
    """Return the cadnano-consistent public name for a scaffold."""
    value = str(name or "").strip()
    return SCAFFOLD_NAME_ALIASES.get(value, value)

# Keep pure-cylinder BILD exports visually identical to the 1.2 3D preview.
# Seed supports use the rod body colours; SST layers use their lattice-volume
# colours.  These values are deliberately centralized so an exported model
# cannot silently drift from the preview palette.
CYLINDER_PREVIEW_COLORS = {
    "seed_z1": "#4f8fce",
    "seed_z2": "#d9dee3",
    "seed_z3": "#d96a82",
    "sst_layer_1": "#2a78d1",
    "sst_layer_2": "#d65b74",
    "default": "#4682b4",
}


def _cylindrical_preview_color(name, parameter,
                               seed_support_ranges=None):
    """Return the 1.2-preview colour for one axial model position."""
    normalized_name = str(name).strip().lower().replace("-", "_")
    if "sst_layer_1" in normalized_name:
        return CYLINDER_PREVIEW_COLORS["sst_layer_1"]
    if "sst_layer_2" in normalized_name:
        return CYLINDER_PREVIEW_COLORS["sst_layer_2"]
    if normalized_name == "seed" and seed_support_ranges and \
            len(seed_support_ranges) >= 2:
        first, second = [tuple(map(float, item))
                         for item in seed_support_ranges[:2]]
        if float(parameter) <= first[1]:
            return CYLINDER_PREVIEW_COLORS["seed_z1"]
        if float(parameter) < second[0]:
            return CYLINDER_PREVIEW_COLORS["seed_z2"]
        return CYLINDER_PREVIEW_COLORS["seed_z3"]
    return CYLINDER_PREVIEW_COLORS["default"]


def _seed_preview_support_ranges(layout):
    """Return the actual Z1/Z3 intervals shown in the 1.2 preview.

    ``seed_layer_ranges`` is the frozen 128/32/128 reference support and is
    retained in structure metadata for Capture-template bookkeeping.  It is
    not the current Seed partition when layer spacing moves the Z2 interval.
    The preview and generated design both serialize that current partition as
    ``seed_partition_ranges`` (Z1, Z2, Z3), including any shared canvas
    translation.  Structure export must use the first and third ranges so its
    twist interpolation and region colours cannot drift from the preview.

    Older accepted designs may not contain ``seed_partition_ranges``.  Keep
    their established ``seed_layer_ranges`` as a compatibility fallback.
    """
    if not isinstance(layout, dict):
        return None
    partition = layout.get("seed_partition_ranges")
    if isinstance(partition, (list, tuple)) and len(partition) == 3:
        try:
            ranges = [tuple(map(int, item)) for item in partition]
        except (TypeError, ValueError):
            ranges = []
        if len(ranges) == 3 and all(
                len(item) == 2 and item[0] <= item[1]
                for item in ranges):
            return [list(ranges[0]), list(ranges[2])]
    legacy = layout.get("seed_layer_ranges")
    if isinstance(legacy, (list, tuple)) and len(legacy) >= 2:
        try:
            ranges = [tuple(map(int, item)) for item in legacy[:2]]
        except (TypeError, ValueError):
            return None
        if all(len(item) == 2 and item[0] <= item[1]
               for item in ranges):
            return [list(item) for item in ranges]
    return None


def _bild_color_command(hex_color):
    value = str(hex_color).strip().lstrip("#")
    if len(value) != 6:
        raise ValueError("BILD colour must use six hexadecimal digits")
    channels = [int(value[index:index+2], 16) / 255.0
                for index in (0, 2, 4)]
    return ".color %.5f %.5f %.5f" % tuple(channels)


def _sequence_position_sort_key(record):
    """Sort sequence rows by the full base integer, then helix integer."""
    start = record.get("start", "") if isinstance(record, dict) else record[0]
    text = str(start).strip()
    match = _SEQUENCE_POSITION_RE.match(text)
    if match is None:
        return (float("inf"), float("inf"), text)
    helix = int(match.group(1))
    base = int(match.group(2))
    return (base, helix, text)


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _sanitized_legacy_payload(payload):
    # Auxiliary h64--79 rows inherit the primary SST row's complete colour
    # list when the short-spacing detour is created.  Most of those colour
    # anchors are intentionally empty: only the colliding Layer-2 interval is
    # routed onto the auxiliary row.  cadnano's legacy decoder assumes every
    # colour anchor has a strand and dereferences ``None.oligo()`` otherwise.
    # Filter only those empty anchors in the in-memory import snapshot; the
    # accepted JSON topology and its visible colours remain unchanged.
    snapshot = deepcopy(payload)
    for row in snapshot.get("vstrands", []):
        staples = row.get("stap", [])
        row["stap_colors"] = [
            item for item in row.get("stap_colors", [])
            if (isinstance(item, (list, tuple)) and len(item) >= 2 and
                0 <= int(item[0]) < len(staples) and
                staples[int(item[0])] != [-1, -1, -1, -1])]
    return snapshot


def _document(payload):
    payload = _sanitized_legacy_payload(payload)
    document = runtime.Document()
    part = import_legacy_dict(
        document, payload, LatticeType.Square, forceLatticeType=True)
    return document, part


def _nodes(oligo):
    return [
        (int(strand.virtualHelix().number()), int(strand.lowIdx()),
         int(strand.highIdx()))
        for strand in oligo.strand5p().generator3pStrand()]


def _target_record(oligo, category, layer=None, color=None):
    strand = oligo.strand5p()
    start, end = oligo.sequenceEndpoints()
    sequence = oligo.sequence() or ""
    length = int(oligo.actualLength())
    return {
        "id": "%s:%d:%d:%d" % (
            category, strand.virtualHelix().number(),
            strand.idx5Prime(), length),
        "category": category,
        "layer": layer,
        "start": start,
        "end": end,
        "start_vh": int(strand.virtualHelix().number()),
        "start_idx": int(strand.idx5Prime()),
        "length": length,
        "color": color or str(oligo.color()),
        "sequence": sequence if sequence and set(sequence.upper()) <= BASES
        else "",
    }


def _sst_physical_layer(nodes, layer_ranges, routing):
    """Classify an SST oligo from the current physical SST layer ranges.

    ``layer_ranges`` contains the generated SST duplex intervals, whereas
    ``seed_layer_ranges`` is only the frozen Seed support/reference.  Using
    the nearest layer *centre* makes the Seed reference behave like a hidden
    boundary when the two SST layers have very different lengths: a strand
    that lies inside the long layer can then be mislabeled as the short layer.

    Assign ordinary strands by their actual base overlap with each current
    SST interval.  Boundary tails that do not overlap either duplex use the
    nearest interval edge, never a layer centre.  At 0/8 bp, the explicit
    auxiliary routing metadata remains authoritative because those h64--79
    strands are the physical detour of Layer 2 even when their base indices
    coincide with the Layer-1 boundary.
    """
    if isinstance(routing, dict) and routing.get("enabled"):
        auxiliary = set(map(int, routing.get(
            "auxiliary_internal_helices", range(64, 80))))
        if any(int(number) in auxiliary
               for number, unused_low, unused_high in nodes):
            return int(routing.get("layer", 2))

    normalized_ranges = [tuple(map(int, item)) for item in layer_ranges]
    overlaps = []
    for layer_low, layer_high in normalized_ranges:
        overlaps.append(sum(
            max(0, min(int(high), layer_high) -
                max(int(low), layer_low) + 1)
            for unused_number, low, high in nodes))
    maximum_overlap = max(overlaps) if overlaps else 0
    if maximum_overlap > 0:
        # A primary strand that is tied exactly across a zero-spacing
        # boundary belongs to Layer 1.  The competing Layer-2 strand is the
        # explicitly tagged auxiliary detour handled above.
        return overlaps.index(maximum_overlap) + 1

    bases = [int(index) for unused_number, low, high in nodes
             for index in (low, high)]
    centre = sum(bases) / max(1, len(bases))

    def distance_to_interval(interval):
        low_value, high_value = interval
        if centre < low_value:
            return low_value - centre
        if centre > high_value:
            return centre - high_value
        return 0.0

    return min(range(len(normalized_ranges)), key=lambda index: (
        distance_to_interval(normalized_ranges[index]), index)) + 1


def _raw_strand_chains(payload, field):
    """Return legacy caDNAno strands as explicit 5-prime-to-3-prime nodes."""
    rows = {int(row["num"]): row for row in payload.get("vstrands", [])}
    nodes = {}
    for number, row in rows.items():
        for index, record in enumerate(row.get(field, [])):
            if list(record) == [-1, -1, -1, -1]:
                continue
            nodes[(number, index)] = (
                (int(record[0]), int(record[1])),
                (int(record[2]), int(record[3])))
    starts = sorted(node for node, (previous, unused_next) in nodes.items()
                    if previous == (-1, -1))
    chains = []
    visited = set()
    for start in starts:
        chain = []
        current = start
        while current != (-1, -1):
            if current in visited:
                raise ValueError(
                    "Legacy %s topology revisits helix %d[%d]." %
                    (field, current[0], current[1]))
            if current not in nodes:
                raise ValueError(
                    "Legacy %s topology points to missing helix %d[%d]." %
                    (field, current[0], current[1]))
            visited.add(current)
            chain.append(current)
            current = nodes[current][1]
        chains.append(chain)
    remaining = sorted(set(nodes) - visited)
    if remaining:
        raise ValueError(
            "Legacy %s topology contains %d bases without a 5-prime start."
            % (field, len(remaining)))
    return rows, chains


def _row_coordinate(row):
    return (int(row["row"]), int(row["col"]))


def _scaffold_base_map(payload, allowed_coordinates=None):
    """Map assigned scaffold/input bases by lattice coordinate and index."""
    rows, chains = _raw_strand_chains(payload, "scaf")
    sequence_by_start = {
        (int(item["start_vh"]), int(item["start_idx"])):
            str(item.get("sequence", "")).upper()
        for item in payload.get("scaffold_sequences", [])}
    base_map = {}
    for chain in chains:
        # Complete-SST remapping supplies ``allowed_coordinates`` precisely
        # so Seed scaffolds are outside the transfer domain.  A Seed chain
        # containing deletions has fewer nucleotides than raw caDNAno base
        # slots, so validating its sequence before this coordinate filter
        # raises a false length error (notably at 8-bp spacing).  Ignore a
        # wholly out-of-domain chain before inspecting its sequence.
        if (allowed_coordinates is not None and not any(
                _row_coordinate(rows[number]) in allowed_coordinates
                for number, unused_index in chain)):
            continue
        sequence = sequence_by_start.get(chain[0])
        if not sequence:
            continue
        if len(sequence) != len(chain):
            raise ValueError(
                "Assigned input at helix %d[%d] has %d bases for a %d-base "
                "topology." % (chain[0][0], chain[0][1], len(sequence),
                               len(chain)))
        for (number, index), base in zip(chain, sequence):
            coordinate = _row_coordinate(rows[number])
            if (allowed_coordinates is not None and
                    coordinate not in allowed_coordinates):
                continue
            if base not in BASES:
                raise ValueError(
                    "Assigned input contains an unresolved base at helix "
                    "%d[%d]." % (number, index))
            key = (coordinate, int(index))
            prior = base_map.get(key)
            if prior is not None and prior != base:
                raise ValueError(
                    "Conflicting input bases at lattice coordinate %s[%d]."
                    % (coordinate, index))
            base_map[key] = base
    return base_map


def _mapped_complete_sst_input_records_from_payload(source, target,
                                                    target_name):
    """Coordinate-map every SST input onto a complete-SST raw topology."""
    target_rows, target_chains = _raw_strand_chains(target, "scaf")
    target_coordinates = {
        _row_coordinate(row) for row in target_rows.values()}
    source_bases = _scaffold_base_map(
        source, allowed_coordinates=target_coordinates)
    records = []
    for chain in target_chains:
        inherited = []
        for number, index in chain:
            coordinate = _row_coordinate(target_rows[number])
            base = source_bases.get((coordinate, int(index)))
            if not base:
                raise ValueError(
                    "%s lacks an accepted SST sublattice input base at lattice "
                    "coordinate %s[%d]." %
                    (target_name, coordinate, index))
            inherited.append(base)
        records.append({
            "start_vh": int(chain[0][0]),
            "start_idx": int(chain[0][1]),
            "sequence": "".join(inherited),
        })
    records.sort(key=lambda item: (
        int(item["start_idx"]), int(item["start_vh"])))
    if not records:
        raise ValueError(
            "%s contains no SST sublattice input strands." % target_name)
    return records


def _design_targets(payload):
    document, part = _document(payload)
    metadata = payload.get("moire_structure_metadata", {})
    sst_first = metadata.get("helix_numbering") == "sst_first"
    sst_helices = set(map(int, metadata.get(
        "sst_helix_numbers", range(16) if sst_first else range(48, 64))))
    seed_helices = set(map(int, metadata.get(
        "seed_helix_numbers", range(16, 64) if sst_first else range(48))))
    layout = metadata.get("variable_length_layout", {})
    layer_ranges = layout.get("layer_ranges", [[48, 175], [208, 335]])
    routing = metadata.get(
        "auxiliary_sst_routing", layout.get("auxiliary_sst_routing", {}))
    targets = {"seed_scaffold": [], "sst_input_layer_1": [],
               "sst_input_layer_2": []}
    scaffolds = []
    for oligo in part.oligos():
        strand = oligo.strand5p()
        if (strand is None or oligo.isLoop() or
                not strand.strandSet().isScaffold()):
            continue
        helices = {number for number, unused_low, unused_high in _nodes(oligo)}
        if helices and helices <= seed_helices:
            scaffolds.append(oligo)
            continue
        if not helices or not helices <= sst_helices:
            continue
        nodes = _nodes(oligo)
        layer = _sst_physical_layer(nodes, layer_ranges, routing)
        key = "sst_input_layer_%d" % layer
        targets[key].append(_target_record(oligo, key, layer))
    scaffolds.sort(key=lambda item: (
        item.strand5p().virtualHelix().number(),
        item.strand5p().idx5Prime()))
    for index, oligo in enumerate(scaffolds):
        targets["seed_scaffold"].append(_target_record(
            oligo, "seed_scaffold", color=SCAFFOLD_COLORS[
                index % len(SCAFFOLD_COLORS)]))
    for key in targets:
        targets[key].sort(key=lambda item: (
            item["start_vh"], item["start_idx"], item["length"]))
    return document, part, targets


def _complete_sst_output_groups(payload, complete_sst_payload=None):
    """Return intact SST-output oligos from the closed export snapshot.

    The accepted physical design deliberately opens every Seed--SST capture
    site.  Reading SST-only staple oligos from that design therefore turns a
    complete 32-nt U-shaped output into two 16-nt fragments.  It also leaves
    one fragment without the original colour marker, which cadnano renders
    with its gray fallback.  Sequence export must instead use the saved
    standalone complete-SST process design (or rebuild that snapshot for an
    old project): the live design remains capture-ready, while the purchasing
    table contains the strands that are actually synthesized.
    """
    # The staged workflow already saves a standalone, closed SST design.  It
    # is the authoritative output topology because it preserves the exact
    # routing that the user inspected and accepted.  Rebuilding is retained
    # only as a compatibility fallback for projects created before that stage
    # file was recorded.
    snapshot = (deepcopy(complete_sst_payload)
                if complete_sst_payload is not None else
                build_complete_sst_only_payload(
                    payload, "complete_sst_output_snapshot.json"))
    snapshot["scaffold_sequences"] = \
        _mapped_complete_sst_input_records_from_payload(
            payload, snapshot, str(snapshot.get(
                "name", "complete_sst_output_snapshot.json")))
    layer_ranges = snapshot.get("moire_structure_metadata", {}).get(
        "variable_length_layout", {}).get(
            "layer_ranges", [[48, 175], [208, 335]])
    metadata = snapshot.get("moire_structure_metadata", {})
    routing = metadata.get("auxiliary_sst_routing", metadata.get(
        "variable_length_layout", {}).get("auxiliary_sst_routing", {}))
    groups = {"sst_output_layer_1": [], "sst_output_layer_2": []}
    rows, chains = _raw_strand_chains(snapshot, "stap")
    scaffold_bases = _scaffold_base_map(snapshot)
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for chain in chains:
        nodes = [(number, index, index) for number, index in chain]
        layer = _sst_physical_layer(nodes, layer_ranges, routing)
        sequence = []
        for number, index in chain:
            row = rows[number]
            coordinate = _row_coordinate(row)
            paired = row.get("scaf", [])[index] != [-1, -1, -1, -1]
            base = scaffold_bases.get((coordinate, int(index)))
            if paired and not base:
                raise ValueError(
                    "Complete SST sublattice output helix %d[%d] lacks its paired "
                    "input base." % (number, index))
            sequence.append(complement[base] if base else "T")
        record = (
            "%d[%d]" % chain[0], "%d[%d]" % chain[-1],
            "".join(sequence), len(chain), "#000000")
        groups["sst_output_layer_%d" % layer].append(record)
    for rows in groups.values():
        rows.sort(key=_sequence_position_sort_key)
    return groups


def _summarize(items):
    counts = Counter(int(item["length"]) for item in items)
    return {
        "count": len(items),
        "lengths": {str(length): count
                    for length, count in sorted(counts.items())},
        "total_nt": sum(int(item["length"]) for item in items),
    }


def analyze(path):
    payload = _load(path)
    unused_document, unused_part, targets = _design_targets(payload)
    return {
        "path": str(Path(path).resolve()),
        "targets": targets,
        "summary": {key: _summarize(value)
                    for key, value in targets.items()},
        "metadata": payload.get("moire_structure_metadata", {}),
    }


def extract_scaffold(path, target):
    payload = _load(path)
    unused_document, unused_part, source_targets = _design_targets(payload)
    candidates = [item for group in source_targets.values() for item in group
                  if item.get("sequence")]
    exact = [item for item in candidates
             if item["start_vh"] == int(target["start_vh"]) and
             item["start_idx"] == int(target["start_idx"]) and
             item["length"] == int(target["length"])]
    if not exact:
        exact = [item for item in candidates
                 if item["length"] == int(target["length"])]
    if len(exact) != 1:
        raise ValueError(
            "所选JSON无法唯一匹配%s（需要%d nt；候选%d条）。" %
            (target.get("start", target.get("id", "scaffold")),
             int(target["length"]), len(exact)))
    selected = exact[0]
    return {
        "target_id": target["id"],
        "sequence": selected["sequence"].upper(),
        "length": selected["length"],
        "source": str(Path(path).resolve()),
        "source_start": selected["start"],
    }


def scaffold_catalog(target_length, multiple=False, used_names=()):
    target_length = int(target_length)
    used_names = {canonical_scaffold_name(item) for item in used_names}
    allowed = (set(MULTI_SCAFFOLD_NAMES) if multiple else
               {name for name, unused_key in STANDARD_SCAFFOLD_KEYS})
    rows = []
    for display_name, sequence_key in STANDARD_SCAFFOLD_KEYS:
        sequence = CADNANO_SCAFFOLDS.get(sequence_key, "")
        if (display_name not in allowed or display_name in used_names or
                len(sequence) < target_length):
            continue
        rows.append({
            "name": display_name,
            "cadnano_key": sequence_key,
            "length": len(sequence),
            "target_length": target_length,
        })
    rows.sort(key=lambda item: (item["length"], item["name"]))
    return {"scaffolds": rows, "multiple": bool(multiple),
            "target_length": target_length}


def assign_standard_scaffold(target, name, multiple=False, used_names=()):
    name = canonical_scaffold_name(name)
    available = scaffold_catalog(
        target["length"], multiple, used_names)["scaffolds"]
    selected = next((item for item in available
                     if item["name"] == name), None)
    if selected is None:
        raise ValueError(
            "The selected scaffold is unavailable, too short, or already "
            "assigned to another route.")
    sequence = CADNANO_SCAFFOLDS[selected["cadnano_key"]]
    target_length = int(target["length"])
    return {
        "target_id": target["id"],
        "start_vh": int(target["start_vh"]),
        "start_idx": int(target["start_idx"]),
        "sequence": sequence[:target_length].upper(),
        "length": target_length,
        "category": "seed_scaffold",
        "layer": None,
        "source": "caDNAno built-in: %s" % selected["cadnano_key"],
        "scaffold_name": selected["name"],
        "scaffold_source_length": selected["length"],
    }


def export_template(design, filename, identical):
    report = analyze(design)
    groups = ["sst_input_layer_1"]
    if not identical:
        groups.append("sst_input_layer_2")
    items = [item for key in groups for item in report["targets"][key]]
    rows = [(item["start"], item["end"], "", item["length"],
             item["color"]) for item in items]
    rows.sort(key=_sequence_position_sort_key)
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    write_sequence_template(str(filename), rows)
    return {"path": str(Path(filename).resolve()), "row_count": len(rows),
            "layers": [1] if identical else [1, 2]}


def _parse_template(design, filename, identical):
    report = analyze(design)
    all_targets = [item for layer in (1, 2)
                   for item in report["targets"][
                       "sst_input_layer_%d" % layer]]
    by_start = {item["start"]: item for item in all_targets}
    headers, rows = read_sequence_template(str(filename))
    if tuple(headers) != HEADERS:
        raise ValueError(
            "The input-template columns must be: %s." % ", ".join(HEADERS))
    assignments = []
    errors = []
    for row_number, values in rows:
        start, end, sequence, length, color = values
        sequence = str(sequence or "").strip().upper()
        if not sequence:
            continue
        target = by_start.get(str(start).strip())
        if target is None:
            errors.append(
                "The Start value in row %d does not belong to a current "
                "SST sublattice input: %s" % (row_number, start))
            continue
        try:
            recorded_length = int(length)
        except (TypeError, ValueError):
            recorded_length = -1
        if str(end).strip() != target["end"]:
            errors.append(
                "The End value in row %d does not match the current "
                "structure." % row_number)
        if recorded_length != target["length"] or \
                len(sequence) != target["length"]:
            errors.append(
                "Row %d requires %d nt; the imported sequence has %d nt." %
                (row_number, target["length"], len(sequence)))
        invalid = sorted(set(sequence) - BASES)
        if invalid:
            errors.append("第%d行包含非法碱基：%s。" %
                          (row_number, "".join(invalid)))
        if not errors or not any(item.startswith("第%d行" % row_number)
                                 for item in errors):
            assignments.append({
                "target_id": target["id"],
                "start_vh": target["start_vh"],
                "start_idx": target["start_idx"],
                "sequence": sequence,
                "length": target["length"],
                "category": target["category"],
                "layer": target["layer"],
                "source": str(Path(filename).resolve()),
            })
    if errors:
        raise ValueError("\n".join(errors))
    if identical:
        layer1 = sorted(
            report["targets"]["sst_input_layer_1"],
            key=lambda item: (item["length"], item["start_vh"],
                              item["start_idx"]))
        layer2 = sorted(
            report["targets"]["sst_input_layer_2"],
            key=lambda item: (item["length"], item["start_vh"],
                              item["start_idx"]))
        if [(item["length"]) for item in layer1] != \
                [(item["length"]) for item in layer2]:
            raise ValueError(
                "The two SST sublattice input topologies differ; sequences "
                "cannot be copied automatically between layers.")
        sequence_by_id = {item["target_id"]: item["sequence"]
                          for item in assignments}
        if len(sequence_by_id) != len(layer1):
            raise ValueError("The 1st-layer template is incomplete.")
        for first, second in zip(layer1, layer2):
            assignments.append({
                "target_id": second["id"],
                "start_vh": second["start_vh"],
                "start_idx": second["start_idx"],
                "sequence": sequence_by_id[first["id"]],
                "length": second["length"],
                "category": second["category"],
                "layer": 2,
                "source": str(Path(filename).resolve()),
                "copied_from_layer_1": True,
            })
    expected = sum(len(report["targets"]["sst_input_layer_%d" % layer])
                   for layer in (1, 2))
    if len(assignments) != expected:
        raise ValueError(
            "SST sublattice input assignment is incomplete: %d strands "
            "required, %d currently assigned." %
            (expected, len(assignments)))
    return assignments


def import_template(design, filename, identical):
    assignments = _parse_template(design, filename, identical)
    return {"assignments": assignments,
            "summary": _summarize(assignments),
            "source": str(Path(filename).resolve())}


def build_sequenced(design, output, assignments):
    payload = _load(design)
    unused_document, unused_part, targets = _design_targets(payload)
    expected = {item["id"]: item for group in targets.values()
                for item in group}
    supplied = {item["target_id"]: item for item in assignments}
    missing = sorted(set(expected) - set(supplied))
    if missing:
        raise ValueError(
            "%d scaffold or SST sublattice input sequences remain "
            "unassigned." % len(missing))
    records = []
    for target in sorted(expected.values(), key=_sequence_position_sort_key):
        target_id = target["id"]
        sequence = str(supplied[target_id]["sequence"]).upper()
        if len(sequence) != int(target["length"]) or \
                set(sequence) - BASES:
            raise ValueError("%s的序列与%d-nt结构不一致。" %
                             (target["start"], target["length"]))
        record = {"start_vh": target["start_vh"],
                  "start_idx": target["start_idx"],
                  "sequence": sequence}
        if target["category"] == "seed_scaffold":
            assignment = supplied[target_id]
            source = str(assignment.get("source", "")).strip()
            scaffold_name = canonical_scaffold_name(
                assignment.get("scaffold_name") or
                (Path(source).stem if source else "Scaffold"))
            total_length = int(
                assignment.get("scaffold_source_length") or
                assignment.get("scaffold_total_length") or
                len(sequence))
            record.update({
                "scaffold_name": scaffold_name,
                "scaffold_source_length": total_length,
                "scaffold_total_length": total_length,
                "scaffold_used_length": int(target["length"]),
            })
            if source:
                record["source"] = source
        records.append(record)
    payload["scaffold_sequences"] = records
    metadata = payload.setdefault("moire_structure_metadata", {})
    metadata["sequence_assignment"] = "accepted in DNA Moire Designer"
    metadata["sequence_assignment_count"] = len(records)
    metadata["sequence_assignment_source"] = str(Path(design).resolve())
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    _write_json(payload, Path(output))
    # Re-import once: this applies complement sequences and catches stale
    # start coordinates before the GUI enables final export.
    unused_document, part = _document(payload)
    output_rows = [
        oligo.sequenceRecord() for oligo in part.oligos()
        if oligo.strand5p() is not None and not oligo.isLoop() and
        oligo.strand5p().strandSet().isStaple()]
    defaulted = sum(str(row[2]).count("?") for row in output_rows)
    resolved_rows = [(
        row[0], row[1], str(row[2]).replace("?", "T"), row[3], row[4])
        for row in output_rows]
    # Bases without a paired imported input are real single-stranded output
    # tails.  They cannot be inferred by complementarity, so make cadnano's
    # existing atom/oxDNA policy explicit: assign T, record every resolved
    # output sequence, and never leak '?' into the final sequence package.
    payload["moire_staple_sequences"] = [
        {"start": row[0], "end": row[1], "sequence": row[2],
         "length": row[3], "color": row[4]}
        for row in resolved_rows]
    metadata["unpaired_output_base_policy"] = "explicit poly-T"
    metadata["defaulted_unpaired_output_bases"] = defaulted
    _write_json(payload, Path(output))
    return {"path": str(Path(output).resolve()),
            "input_sequence_count": len(records),
            "output_sequence_count": len(output_rows),
            "unresolved_output_bases": 0,
            "defaulted_unpaired_output_bases": defaulted}


def _empty_record():
    return [-1, -1, -1, -1]


def _subset_payload(source, keep_helices, base_range=None, name="subset.json"):
    payload = deepcopy(source)
    keep_helices = set(map(int, keep_helices))
    rows = []
    for original in payload.get("vstrands", []):
        number = int(original["num"])
        if number not in keep_helices:
            continue
        row = deepcopy(original)
        for field in ("scaf", "stap"):
            for index, record in enumerate(row.get(field, [])):
                if base_range is not None and not (
                        base_range[0] <= index <= base_range[1]):
                    row[field][index] = _empty_record()
                    continue
                for offset in (0, 2):
                    partner, partner_base = map(
                        int, record[offset:offset + 2])
                    if (partner not in keep_helices or
                            (base_range is not None and not
                             (base_range[0] <= partner_base <= base_range[1]))):
                        record[offset:offset + 2] = [-1, -1]
        row["stap_colors"] = [
            [int(index), int(color)]
            for index, color in row.get("stap_colors", [])
            if (base_range is None or
                base_range[0] <= int(index) <= base_range[1]) and
            int(index) < len(row.get("stap", [])) and
            row["stap"][int(index)] != _empty_record()]
        rows.append(row)
    payload["vstrands"] = rows
    payload["name"] = name
    payload["scaffold_sequences"] = [
        deepcopy(item) for item in payload.get("scaffold_sequences", [])
        if int(item.get("start_vh", -1)) in keep_helices and
        (base_range is None or
         base_range[0] <= int(item.get("start_idx", -1)) <= base_range[1])]
    return payload


def _write_design_pair(source_path, target_root, stem, assignments=None):
    payload = _load(source_path)
    no_sequence = deepcopy(payload)
    no_sequence.pop("scaffold_sequences", None)
    no_sequence_path = target_root / (stem + "_no_sequence.json")
    with_sequence_path = target_root / (stem + "_with_sequence.json")
    _write_json(no_sequence, no_sequence_path)
    if assignments is not None:
        payload["scaffold_sequences"] = deepcopy(assignments)
    _write_json(payload, with_sequence_path)
    no_sequence_svg = target_root / (stem + "_no_sequence.svg")
    with_sequence_svg = target_root / (stem + "_with_sequence.svg")
    _render_svg(no_sequence_path, no_sequence_svg)
    _render_svg(with_sequence_path, with_sequence_svg)
    return {"no_sequence": str(no_sequence_path),
            "with_sequence": str(with_sequence_path),
            "no_sequence_svg": str(no_sequence_svg),
            "with_sequence_svg": str(with_sequence_svg)}


def _mapped_input_records(source, target_path, groups):
    """Map accepted input sequences onto a stage's own helix numbering."""
    unused_document, unused_part, source_targets = _design_targets(source)
    target = _load(target_path)
    unused_document, unused_part, target_targets = _design_targets(target)
    sequence_by_start = {
        (int(item["start_vh"]), int(item["start_idx"])):
            str(item["sequence"]).upper()
        for item in source.get("scaffold_sequences", [])}
    records = []
    for group in groups:
        source_items = source_targets[group]
        target_items = target_targets[group]
        if [item["length"] for item in source_items] != [
                item["length"] for item in target_items]:
            raise ValueError("%s阶段与最终结构的%s input拓扑不一致。" %
                             (Path(target_path).name, group))
        for source_item, target_item in zip(source_items, target_items):
            sequence = sequence_by_start.get((
                source_item["start_vh"], source_item["start_idx"]))
            if not sequence:
                raise ValueError("%s缺少%s的已接受序列。" %
                                 (Path(target_path).name,
                                  source_item["start"]))
            records.append({
                "start_vh": target_item["start_vh"],
                "start_idx": target_item["start_idx"],
                "sequence": sequence})
    records.sort(key=lambda item: (
        int(item["start_idx"]), int(item["start_vh"])))
    return records


def _mapped_complete_sst_input_records(source, target_path):
    """Map final SST inputs onto the standalone complete-SST topology.

    At very short layer spacings (notably 0 bp), the production generator
    routes colliding layer-2 domains through auxiliary SST helices.  Closing
    the capture sites for the standalone SST file can consequently repartition
    those same physical bases into a different set of 32/48-nt oligos.  A
    positional zip of oligos is therefore invalid even though both designs are
    correct.

    Transfer each nucleotide by its lattice coordinate and base index instead.
    Coordinates are deliberately used rather than virtual-helix numbers:
    caDNAno's legacy importer compacts sparse standalone h64--79 identifiers to
    h16--31, while the row/column coordinates remain invariant.
    """
    target = _load(target_path)
    return _mapped_complete_sst_input_records_from_payload(
        source, target, Path(target_path).name)


def _seed_template_component_colors(payload):
    """Map every Seed base to its immutable reference-component colour."""
    frozen = payload_to_sst_first_numbering(
        _load_json(SEED_CAPTURE_REFERENCE))
    marker_colors = {
        (int(row["num"]), int(index)): "#%06x" % int(color)
        for row in frozen.get("vstrands", [])
        for index, color in row.get("stap_colors", [])
    }
    frozen_rows = {
        int(row["num"]): row for row in frozen.get("vstrands", [])}
    template_components, unused_labels, unused_adjacency = \
        runtime._staple_components_from_rows(frozen_rows)
    layout = payload.get("moire_structure_metadata", {}).get(
        "variable_length_layout", {})
    coordinate_shift = int(layout.get("coordinate_shift_bp", 0) or 0)
    if coordinate_shift < 0 or coordinate_shift % 32:
        raise ValueError(
            "Seed template colour mapping requires a non-negative 32-bp "
            "canvas shift; received %d bp." % coordinate_shift)

    result = {}
    for component in template_components:
        component_markers = {
            marker_colors[node]
            for node in component if node in marker_colors
        }
        if len(component_markers) != 1:
            raise ValueError(
                "The immutable Seed template must contain exactly one "
                "staple-colour marker per component; found %d." %
                len(component_markers))
        component_color = next(iter(component_markers)).lower()
        for helix, index in component:
            result[(int(helix), int(index) + coordinate_shift)] = \
                component_color
    return result


def _seed_template_capture_metadata(payload):
    """Map every immutable Capture-core base to its real Capture endpoint.

    An oligo's 5' coordinate is not necessarily the endpoint that can extend
    into SST.  In particular, the potential-Z2 component reported formerly
    as ``Seed 52[63]`` actually terminates on a Capture-face helix elsewhere
    in that component.  Resolve the endpoint from the frozen column phase
    and Seed face mapping, then translate only the base coordinate with the
    complete Seed canvas shift.
    """
    frozen = payload_to_sst_first_numbering(
        _load_json(SEED_CAPTURE_REFERENCE))
    marker_colors = {
        (int(row["num"]), int(index)): int(color)
        for row in frozen.get("vstrands", [])
        for index, color in row.get("stap_colors", [])}
    frozen_rows = {
        int(row["num"]): row for row in frozen.get("vstrands", [])}
    template_components, unused_labels, unused_adjacency = \
        runtime._staple_components_from_rows(frozen_rows)
    layout = payload.get("moire_structure_metadata", {}).get(
        "variable_length_layout", {})
    coordinate_shift = int(layout.get("coordinate_shift_bp", 0) or 0)
    if coordinate_shift < 0 or coordinate_shift % 32:
        raise ValueError(
            "Seed template Capture mapping requires a non-negative 32-bp "
            "canvas shift; received %d bp." % coordinate_shift)
    seed_helices = set(range(16, 64))
    result = {}
    for component in template_components:
        seed_component = {
            node for node in component if int(node[0]) in seed_helices}
        if not seed_component:
            continue
        component_markers = {
            marker_colors[node]
            for node in component if node in marker_colors}
        if len(component_markers) != 1:
            raise ValueError(
                "The immutable Seed template must contain exactly one "
                "staple-colour marker per component; found %d." %
                len(component_markers))
        reference_color = next(iter(component_markers))
        if "#%06x" % reference_color in _SEED_NORMAL_TEMPLATE_COLORS:
            continue
        logical_columns = ((184, 200) if reference_color == 0x999999 else
                           (CAPTURE_REFERENCE_COLUMN_BY_COLOR.get(
                               reference_color),))
        if logical_columns == (None,):
            raise ValueError(
                "Immutable Seed Capture colour #%06x has no column." %
                reference_color)
        matches = []
        for logical_column in logical_columns:
            unit_index = (int(logical_column) -
                          CAPTURE_DIRECT_POSITIONS[0]) // 16
            phase = "B" if unit_index % 2 == 0 else "A"
            # The reference payload is already SST-first: internal Seed
            # helix h is therefore displayed/exported as h+16.
            capture_helices = {
                int(seed_helix) + 16
                for cycle in (phase + "0", phase + "1")
                for unused_sst, seed_helix in
                CAPTURE_PHASE_MAPPINGS[cycle]}
            for helix in capture_helices:
                if (int(helix), int(logical_column)) in component:
                    matches.append((int(logical_column), int(helix)))
        if len(matches) != 1:
            raise ValueError(
                "Immutable Seed Capture core does not map to exactly one "
                "real endpoint: %s." % matches)
        logical_column, capture_helix = matches[0]
        capture_base = logical_column + coordinate_shift
        metadata = {
            "reference_color": "#%06x" % reference_color,
            "capture_seed_helix": capture_helix,
            "capture_base": capture_base,
            "column_base": capture_base,
            "column_index": capture_column_index(capture_base, layout),
        }
        if metadata["column_index"] is None:
            raise ValueError(
                "Immutable Seed Capture endpoint %d[%d] has no column." %
                (capture_helix, capture_base))
        for helix, index in seed_component:
            result[(int(helix), int(index) + coordinate_shift)] = metadata
    return result


def _template_capture_metadata_for_nodes(nodes, metadata_by_node):
    records = {
        (item["capture_seed_helix"], item["capture_base"]): item
        for helix, low, high in nodes
        for index in range(int(low), int(high) + 1)
        for item in [metadata_by_node.get((int(helix), index))]
        if item is not None}
    if len(records) == 1:
        return next(iter(records.values()))
    return None


def _template_color_for_nodes(nodes, template_component_colors):
    colors = {
        template_component_colors[(int(helix), int(index))]
        for helix, low, high in nodes
        for index in range(int(low), int(high) + 1)
        if (int(helix), int(index)) in template_component_colors
    }
    if len(colors) == 1:
        return next(iter(colors))
    return None


def _apply_seed_template_capture_colors(
        manifest, template_component_colors):
    """Use the accepted Seed core colour without losing extension colour."""
    for item in manifest:
        key = (int(item["seed_helix"]), int(item["base"]))
        template_color = template_component_colors.get(key)
        if (template_color is None or
                template_color in _SEED_NORMAL_TEMPLATE_COLORS):
            raise ValueError(
                "Capture at Seed %d[%d] does not map to an immutable "
                "Capture core in Square_Seed_2L_newtemplate.json." % key)
        # Gray #999999 components are immutable *potential* Z2 Capture
        # cores, not ordinary Seed staples.  When the current SST placement
        # physically reaches one of them it becomes an actual Capture and
        # must be exported here.  When it remains unconnected it is emitted
        # separately below as a potential-Z2 core.  Rejecting gray at this
        # point made valid short-spacing designs (for example Seed 45[200])
        # fail only during final export.
        extension_color = str(item.get(
            "capture_extension_color", item.get("capture_color", "#000000")))
        # The row identity and the coloured Capture extension must describe
        # the same current column.  The Seed core is still rendered black;
        # its historical template marker is used only to validate that this
        # is a real Capture component.
        item["color"] = extension_color
        item["staple_core_color"] = "#000000"
        item["capture_color"] = extension_color
        item["capture_extension_color"] = extension_color
        for run in item.get("sequence_color_runs", ()):
            if run.get("role") == "staple_core":
                run["color"] = "#000000"
            elif run.get("role") == "capture_extension":
                run["color"] = extension_color


def _apply_capture_column_colors(manifest, layout):
    """Synchronize row, map and extension colours by immutable column."""
    for item in manifest:
        column_base = int(item.get("column_base", item["base"]))
        column_index = capture_column_index(column_base, layout)
        color_value = capture_column_color(column_base, layout)
        if column_index is None or color_value is None:
            raise ValueError(
                "Capture endpoint Seed %d[%d] does not map to an immutable "
                "real-space column." % (
                    int(item.get("capture_seed_helix",
                                 item["seed_helix"])),
                    int(item.get("capture_base", item["base"]))))
        pair_index = item.get("pair_index")
        if pair_index is None:
            pair_index = capture_pair_index(column_base, layout)
        color = "#%06x" % int(color_value)
        item["column_base"] = column_base
        item["column_index"] = int(column_index)
        if pair_index is not None:
            item["pair_index"] = int(pair_index)
        item["color"] = color
        item["capture_color"] = color
        has_extension = False
        for run in item.get("sequence_color_runs", ()):
            if run.get("role") == "capture_extension":
                run["color"] = color
                has_extension = True
        if has_extension or item.get("capture_extension_color"):
            item["capture_extension_color"] = color


def _sequence_sheets(payload, complete_sst_payload=None):
    document, part, targets = _design_targets(payload)
    metadata = payload.get("moire_structure_metadata", {})
    sst_first = metadata.get("helix_numbering") == "sst_first"
    sst_helices = set(map(int, metadata.get(
        "sst_helix_numbers", range(16) if sst_first else range(48, 64))))
    seed_helices = set(map(int, metadata.get(
        "seed_helix_numbers", range(16, 64) if sst_first else range(48))))
    layout = metadata.get("variable_length_layout", {})
    lattice_type = metadata.get(
        "lattice_type", layout.get("lattice_type", "square"))
    layer_ranges = layout.get("layer_ranges", [[48, 175], [208, 335]])
    # Capture membership and Seed-core colour come only from the immutable
    # accepted template. Runtime/default caDNAno colours (notably #888888)
    # must never turn ordinary Seed staples into Capture products.
    template_component_colors = _seed_template_component_colors(payload)
    template_capture_metadata_by_node = \
        _seed_template_capture_metadata(payload)
    input_by_start = {
        (item["start_vh"], item["start_idx"]): item
        for group in targets.values() for item in group}
    input_groups = {"scaffold": [], "sst_input_layer_1": [],
                    "sst_input_layer_2": []}
    for oligo in part.oligos():
        strand = oligo.strand5p()
        if strand is None or not strand.strandSet().isScaffold():
            continue
        key = (int(strand.virtualHelix().number()), int(strand.idx5Prime()))
        target = input_by_start.get(key)
        if target is None:
            continue
        category = ("scaffold" if target["category"] == "seed_scaffold"
                    else target["category"])
        input_groups[category].append(tuple(oligo.sequenceRecord()))
    outputs = _complete_sst_output_groups(
        payload, complete_sst_payload=complete_sst_payload)
    outputs.update({"normal_staple": [], "staple_capture": []})
    z2_potential_capture_manifest = []
    unextended_capture_candidates = []
    for oligo in part.oligos():
        strand = oligo.strand5p()
        if strand is None or not strand.strandSet().isStaple() or oligo.isLoop():
            continue
        nodes = _nodes(oligo)
        helices = {item[0] for item in nodes}
        raw_record = tuple(oligo.sequenceRecord())
        record = (raw_record[0], raw_record[1],
                  str(raw_record[2]).replace("?", "T"),
                  raw_record[3], raw_record[4])
        if helices & seed_helices and helices & sst_helices:
            outputs["staple_capture"].append(record)
        elif helices and helices <= sst_helices:
            # SST outputs are read from the intact, capture-closed snapshot
            # above.  The physical design intentionally contains gaps and
            # must never be used as the purchasing-strand topology.
            continue
        elif helices and helices <= seed_helices:
            template_color = _template_color_for_nodes(
                nodes, template_component_colors)
            capture_metadata = _template_capture_metadata_for_nodes(
                nodes, template_capture_metadata_by_node)
            # The immutable Seed template contains sixteen gray (#999999)
            # Z2 cores reserved as future Capture sites.  They belong to the
            # staple-capture product even when no SST extension is physically
            # connected.  The other non-black Seed-only cores are the 32-bp
            # translated capture set and are exported below from the
            # translated manifest, so they must not be duplicated here.
            if template_color == _SEED_POTENTIAL_Z2_COLOR:
                if capture_metadata is None:
                    raise ValueError(
                        "Potential-Z2 core does not map to a real Capture "
                        "endpoint in the immutable Seed template.")
                seed_helix = int(
                    capture_metadata["capture_seed_helix"])
                seed_base = int(capture_metadata["capture_base"])
                z2_potential_capture_manifest.append({
                    # Private component geometry used below to suppress a
                    # potential-Z2 row when this same immutable core is
                    # already represented by a physical or translated
                    # Capture.  It is removed before the manifest leaves
                    # this function.
                    "_nodes": nodes,
                    "face": "potential_z2",
                    "seed_helix": seed_helix,
                    "seed_internal_helix": seed_helix - 16,
                    "sst_helix": -1,
                    "sst_internal_helix": -1,
                    "base": seed_base,
                    "capture_seed_helix": seed_helix,
                    "capture_base": seed_base,
                    "column_base": seed_base,
                    "column_index": int(
                        capture_metadata["column_index"]),
                    "pair_index": None,
                    "layer": None,
                    "phase": "Z2",
                    "translation": "potential",
                    "cycle": "Z2",
                    "export_only": True,
                    "connection_role": (
                        "immutable potential Z2 capture core; no current "
                        "SST sublattice extension"),
                    "start": record[0], "end": record[1],
                    "sequence": record[2], "length": record[3],
                    "color": "#999999",
                    "staple_core_color": "#000000",
                    "capture_color": "#999999",
                    "sequence_color_runs": [{
                        "start": 0, "end": int(record[3]),
                        "role": "potential_capture_core",
                        "color": "#000000",
                    }],
                })
            elif (template_color is None or
                  template_color in _SEED_NORMAL_TEMPLATE_COLORS):
                outputs["normal_staple"].append(
                    (record[0], record[1], record[2], record[3],
                     "#000000"))
            else:
                # Every non-black core in the immutable Seed belongs to the
                # staple-capture product family.  Usually it is represented
                # below by the 32-bp translated export payload.  Kagome SST
                # holes deliberately leave some of those frozen Seed cores
                # without a physical/translated extension; retain them as
                # bare cores instead of silently dropping them from export.
                if capture_metadata is None:
                    raise ValueError(
                        "Immutable Capture core does not map to a real "
                        "Capture endpoint.")
                unextended_capture_candidates.append({
                    "record": record,
                    "nodes": nodes,
                    "seed_helix": int(
                        capture_metadata["capture_seed_helix"]),
                    "seed_base": int(capture_metadata["capture_base"]),
                    "column_index": int(
                        capture_metadata["column_index"]),
                    "color": template_color,
                })
    translated = runtime.build_translated_capture_export_payload(payload)
    translated_manifest = _capture_sequence_manifest(
        _sanitized_legacy_payload(translated))
    physical_manifest = _capture_sequence_manifest(
        _sanitized_legacy_payload(payload))
    _apply_seed_template_capture_colors(
        physical_manifest + translated_manifest,
        template_component_colors)
    for item in physical_manifest + translated_manifest:
        item["sequence"] = str(item.get("sequence", "")).replace("?", "T")
    represented_capture_targets = {
        (int(item["seed_helix"]), int(item["base"]))
        for item in physical_manifest + translated_manifest}

    # A gray potential-Z2 core becomes an ordinary physical/translated
    # Capture as soon as the current SST placement reaches any base in that
    # immutable component.  Do not emit the same component once more as an
    # unconnected potential core.  Matching every base is essential because
    # the Capture coordinate is not necessarily the component's 5' start.
    filtered_z2_potential_manifest = []
    for item in z2_potential_capture_manifest:
        nodes = item.pop("_nodes", ())
        represented = any(
            (int(helix), int(base)) in represented_capture_targets
            for helix, low, high in nodes
            for base in range(int(low), int(high) + 1))
        if not represented:
            filtered_z2_potential_manifest.append(item)
    z2_potential_capture_manifest = filtered_z2_potential_manifest

    unextended_capture_manifest = []
    for candidate in unextended_capture_candidates:
        # A physical or translated target can sit anywhere inside the frozen
        # core, not necessarily at its 5' start.  Match against every Seed
        # segment of the oligo so an already represented core is never
        # duplicated.
        represented = any(
            (int(helix), int(base)) in represented_capture_targets
            for helix, low, high in candidate["nodes"]
            for base in range(int(low), int(high) + 1))
        if represented:
            continue
        record = candidate["record"]
        seed_helix = int(candidate["seed_helix"])
        seed_internal = seed_helix - 16 if sst_first else seed_helix
        color = str(candidate["color"])
        unextended_capture_manifest.append({
            "face": "face1" if seed_internal < 8 else "face2",
            "seed_helix": seed_helix,
            "seed_internal_helix": seed_internal,
            "sst_helix": -1,
            "sst_internal_helix": -1,
            "base": int(candidate["seed_base"]),
            "capture_seed_helix": seed_helix,
            "capture_base": int(candidate["seed_base"]),
            "column_base": int(candidate["seed_base"]),
            "column_index": int(candidate["column_index"]),
            "pair_index": None,
            "layer": None,
            "phase": "unextended",
            "translation": "unextended",
            "cycle": "unextended",
            "export_only": True,
            "connection_role": (
                "immutable unextended capture core; Kagome SST sublattice hole"
                if str(lattice_type).lower() == "kagome" else
                "immutable unextended capture core; no current SST sublattice "
                "extension"),
            "start": record[0], "end": record[1],
            "sequence": record[2], "length": record[3],
            "color": color,
            "staple_core_color": "#000000",
            "capture_color": color,
            "sequence_color_runs": [{
                "start": 0, "end": int(record[3]),
                "role": "unextended_capture_core",
                "color": "#000000",
            }],
        })
    manifest = (physical_manifest + translated_manifest +
                unextended_capture_manifest +
                z2_potential_capture_manifest)
    _apply_capture_column_colors(manifest, layout)
    manifest = _sort_capture_manifest_by_template_color(manifest)
    # The final capture sheet is deliberately rebuilt from the physical and
    # export-only capture manifest.  This both removes duplicate bridge rows
    # and preserves the product rule requested by the UI: template colour
    # first; numeric base then helix inside each colour block.
    capture_seen = set()
    capture_rows = []
    ordered_manifest = []
    for item in manifest:
        identity = (item["start"], item["end"], item["sequence"])
        if identity in capture_seen:
            continue
        capture_seen.add(identity)
        ordered_manifest.append(item)
        capture_rows.append(_capture_staple_row(item))
    outputs["staple_capture"] = capture_rows
    for key, rows in input_groups.items():
        rows.sort(key=_sequence_position_sort_key)
    for key, rows in outputs.items():
        if key != "staple_capture":
            rows.sort(key=_sequence_position_sort_key)
    return input_groups, outputs, ordered_manifest


def _scaffold_export_summary(payload, assignments=()):
    """Return user-facing metadata for every routed scaffold sequence."""
    unused_document, unused_part, targets = _design_targets(payload)
    saved_by_start = {
        (int(item.get("start_vh", -1)), int(item.get("start_idx", -1))): item
        for item in payload.get("scaffold_sequences", [])}
    assigned_by_start = {
        (int(item.get("start_vh", -1)), int(item.get("start_idx", -1))): item
        for item in assignments
        if item.get("category") == "seed_scaffold" or
        item.get("scaffold_name")}
    summary = []
    for index, target in enumerate(targets["seed_scaffold"], 1):
        key = (int(target["start_vh"]), int(target["start_idx"]))
        saved = saved_by_start.get(key, {})
        assignment = assigned_by_start.get(key, {})
        metadata = dict(saved)
        metadata.update(assignment)
        source = str(metadata.get("source", "")).strip()
        name = canonical_scaffold_name(
            metadata.get("scaffold_name") or
            (Path(source).stem if source and not source.startswith(
                "caDNAno built-in:") else "Scaffold %d" % index))
        sequence = str(saved.get("sequence") or
                       assignment.get("sequence") or
                       target.get("sequence") or "").upper()
        used_length = int(target["length"])
        total_length = int(
            metadata.get("scaffold_source_length") or
            metadata.get("scaffold_total_length") or
            max(used_length, len(sequence)))
        summary.append({
            "start": target["start"],
            "end": target["end"],
            "scaffold_name": name,
            "total_scaffold_length_nt": total_length,
            "length_used_in_structure_nt": used_length,
        })
    summary.sort(key=_sequence_position_sort_key)
    return summary


def _xlsx_cell(reference, value, header=False):
    style = ' s="1"' if header else ""
    if isinstance(value, int):
        return '<c r="%s"%s><v>%d</v></c>' % (reference, style, value)
    return ('<c r="%s" t="inlineStr"%s><is><t>%s</t></is></c>' %
            (reference, style, escape(str(value))))


def _add_scaffold_metadata_columns(workbook, scaffold_summary):
    """Append scaffold identity and length metadata to the scaffold sheet."""
    path = Path(workbook)
    temporary_path = None
    headers = (
        "Scaffold Name",
        "Total Scaffold Length (nt)",
        "Length Used in Structure (nt)",
    )
    try:
        with tempfile.NamedTemporaryFile(
                dir=str(path.parent), prefix=path.stem + "-scaffold-",
                suffix=path.suffix, delete=False) as temporary:
            temporary_path = Path(temporary.name)
        with ZipFile(path, "r") as source_archive, ZipFile(
                temporary_path, "w", ZIP_DEFLATED) as target_archive:
            for member in source_archive.infolist():
                payload = source_archive.read(member.filename)
                if member.filename == "xl/worksheets/sheet1.xml":
                    xml = payload.decode("utf-8")
                    last_row = max(1, len(scaffold_summary) + 1)
                    xml = re.sub(
                        r'<dimension ref="A1:E\d+"\s*/>',
                        '<dimension ref="A1:H%d"/>' % last_row, xml)
                    xml = re.sub(
                        r'<autoFilter ref="A1:E\d+"\s*/>',
                        '<autoFilter ref="A1:H%d"/>' % last_row, xml)
                    xml = xml.replace(
                        "</cols>",
                        '<col min="6" max="6" width="20" customWidth="1"/>'
                        '<col min="7" max="8" width="29" customWidth="1"/>'
                        "</cols>", 1)
                    row_values = [headers] + [(
                        item["scaffold_name"],
                        item["total_scaffold_length_nt"],
                        item["length_used_in_structure_nt"],
                    ) for item in scaffold_summary]
                    for row_number, values in enumerate(row_values, 1):
                        extra = "".join(
                            _xlsx_cell("%s%d" % (column, row_number), value,
                                       header=row_number == 1)
                            for column, value in zip("FGH", values))
                        pattern = r'(<row\b[^>]*\br="%d"[^>]*>)(.*?)(</row>)' % row_number
                        xml, count = re.subn(
                            pattern, lambda match: (
                                match.group(1) + match.group(2) + extra +
                                match.group(3)), xml, count=1)
                        if count != 1:
                            raise ValueError(
                                "Scaffold workbook row %d is missing" %
                                row_number)
                    payload = xml.encode("utf-8")
                target_archive.writestr(member, payload)
        temporary_path.replace(path)
    except Exception:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except OSError:
                pass
        raise
    return str(path)


def _project_auxiliary_geometry(records, payload):
    """Project real auxiliary sequence channels onto ideal SST axes.

    caDNAno h64--79 remain the authoritative helix/base identifiers in every
    sequence and mapping record.  Only their spatial coordinates are moved,
    matching the ideal 3D preview where those channels are routing detours
    rather than sixteen extra duplex rods.
    """
    metadata = payload.get("moire_structure_metadata", {})
    routing = metadata.get("auxiliary_sst_routing", {})
    if not (isinstance(routing, dict) and routing.get("enabled")):
        return 0
    sst_first = metadata.get("helix_numbering") == "sst_first"
    row_coords = {
        int(row["num"]): (int(row["row"]), int(row["col"]))
        for row in payload.get("vstrands", [])}
    number_by_coord = {coord: number for number, coord in row_coords.items()}
    projected = 0
    for actual in range(64, 80):
        ideal = actual - (64 if sst_first else 16)
        if actual not in row_coords or ideal not in row_coords:
            continue
        for record in records:
            strand = record.get("strand")
            if strand is None:
                continue
            vh = strand.virtualHelix()
            record_actual = number_by_coord.get(tuple(vh.coord()))
            if record_actual != actual:
                continue
            part = vh.part()
            actual_xy = part.latticeCoordToPositionXY(*row_coords[actual])
            ideal_xy = part.latticeCoordToPositionXY(*row_coords[ideal])
            scale = 2.8 / (2.0 * float(part.radius())) / \
                _oxdna.OXDNA_LENGTH_NM
            delta = ((actual_xy[0] - ideal_xy[0]) * scale,
                     (actual_xy[1] - ideal_xy[1]) * scale, 0.0)
            for coordinate_key in ("pos", "output_pos"):
                if coordinate_key in record:
                    record[coordinate_key] = tuple(
                        record[coordinate_key][axis] - delta[axis]
                        for axis in range(3))
            record["cadnano_helix"] = actual
            record["ideal_geometry_helix"] = ideal
            projected += 1
    return projected


def _export_structure(payload, root, name, seed_twist_angle_deg=None,
                      seed_support_ranges=None):
    document, unused_part = _document(payload)
    # Legacy PDB has only 62 portable chain IDs; force mmCIF beyond that.
    strand_count = len([oligo for oligo in document.oligos()
                        if oligo.strand5p() is not None])
    atom_limit = 99999 if strand_count <= 62 else 0
    routing = payload.get("moire_structure_metadata", {}).get(
        "auxiliary_sst_routing", {})
    has_auxiliary_geometry = bool(
        isinstance(routing, dict) and routing.get("enabled"))
    if (seed_twist_angle_deg is None or not seed_support_ranges) and not \
            has_auxiliary_geometry:
        records, unused_strands, unused_assigned, unused_residual = \
            _oxdna._collect(document, 2.8)
        report = export_structure_bundle(
            document, str(root), name, spacing_nm=2.8,
            include_oxdna=True, pdb_atom_limit=atom_limit)
        return _add_cylindrical_model(
            report, records, root, name,
            seed_support_ranges=seed_support_ranges)

    # Reproduce the 1.2 preview geometry in the exported Seed: the first
    # support remains fixed, the Seed Z2 spacer rotates continuously, and the
    # second support remains at the requested relative twist.  Applying the
    # same transform to positions and local nucleotide frames keeps PDB/CIF
    # and oxDNA DAT geometrically consistent.
    records, strands, assigned, residual = _oxdna._collect(document, 2.8)

    # h64--79 are real caDNAno sequence channels, but not extra physical
    # helices.  Translate their nucleotide coordinates onto the ideal SST
    # h0--15 (public SST-first) / h48--63 (internal) axes.  Strand/vh/base
    # identities are deliberately untouched, so PDB/CIF mapping and every
    # sequence table retain the exact caDNAno numbering.
    projected_auxiliary = _project_auxiliary_geometry(records, payload)

    total_angle = None
    if seed_twist_angle_deg is not None and seed_support_ranges:
        first, second = [tuple(map(float, item))
                         for item in seed_support_ranges[:2]]
        lower_end = first[1]
        upper_start = second[0]
        total_angle = float(seed_twist_angle_deg)
        for record in records:
            index = float(record.get("idx", 0))
            if index <= lower_end:
                fraction = 0.0
            elif index >= upper_start:
                fraction = 1.0
            else:
                fraction = ((index - lower_end) /
                            max(1.0, upper_start - lower_end))
            angle = total_angle * fraction * 3.141592653589793 / 180.0
            cosine, sine = math.cos(angle), math.sin(angle)

            def rotate_xy(vector):
                x, y, z = vector
                return (cosine*x - sine*y, sine*x + cosine*y, z)

            record["pos"] = rotate_xy(record["pos"])
            if "output_pos" in record:
                record["output_pos"] = rotate_xy(record["output_pos"])
            record["a1"] = rotate_xy(record["a1"])
            record["a2"] = rotate_xy(record["a2"])
            record["a3"] = rotate_xy(record["a3"])

    ordered = _oxdna._number_records(strands)
    dat_text = _oxdna._dat_text(ordered)
    atoms, unused_mapping, remarks = _oxdna._all_atom_structure(
        strands, 2.8, residual)
    if total_angle is not None:
        remarks.append(
            "MOIRE SEED RELATIVE TWIST %.6f DEG APPLIED ACROSS SEED Z2" %
            total_angle)
    if projected_auxiliary:
        remarks.append(
            "AUXILIARY caDNAno HELICES 64-79 PROJECTED TO IDEAL SST SUBLATTICE "
            "GEOMETRY; ORIGINAL HELIX/BASE IDS RETAINED")
    paths = _oxdna.structure_bundle_paths(str(root), name, True)
    Path(paths["root"]).mkdir(parents=True, exist_ok=True)
    cylindrical_path = Path(paths["root"]) / (
        name + "_cylindrical_model.bild")
    paths["cylindrical_model"] = str(cylindrical_path)
    if len(atoms) <= atom_limit:
        structure_format = "PDB"
        structure_path = paths["all_atom"]["pdb"]
        files = {structure_path: _oxdna._pdb_text_from_atoms(atoms, remarks)}
    else:
        structure_format = "mmCIF"
        structure_path = paths["all_atom"]["cif"]
        files = {structure_path: _oxdna._mmcif_text(
            atoms, name + "_all_atom",
            "DNA Moire Seed with preview-matched relative twist")}
    files[paths["oxdna"]["top"]] = _oxdna._top_text(
        ordered, len(strands))
    files[paths["oxdna"]["dat"]] = dat_text
    files[str(cylindrical_path)] = _cylindrical_model_bild(
        records, name, seed_support_ranges=seed_support_ranges)
    _oxdna._atomic_write_many(files)
    return {"paths": paths, "nucleotides": len(ordered),
            "strands": len(strands), "assigned_bases": assigned,
            "all_atom_count": len(atoms),
            "structure_path": structure_path,
            "structure_format": structure_format,
            "include_oxdna": True, "hybrid_residual": residual,
            "cylindrical_model_path": str(cylindrical_path),
            "moire_seed_twist_angle_deg": total_angle,
            "moire_seed_twist_profile": (
                "fixed / linear Z2 / fixed" if total_angle is not None
                else None),
            "auxiliary_geometry_projected": bool(projected_auxiliary),
            "auxiliary_projected_nucleotide_count": projected_auxiliary,
            "sequence_numbering_policy": (
                "retain actual cadnano helix/base ids")}


def _cylindrical_model_bild(records, name, radius_nm=1.125,
                            radial_subdivisions=48,
                            seed_support_ranges=None):
    """Render occupied helix axes as watertight high-resolution BILD tubes.

    Axis points are reconstructed from the already positioned oxDNA records.
    Consequently, any Seed twist applied to ``pos`` and ``a1`` is inherited
    instead of being independently recalculated for this model.  BILD's
    native ``.cylinder`` primitive has a fixed, coarse circumference and
    gives every base-to-base segment its own end caps.  Those caps produce
    visible seams wherever a twisted axis bends.  The exporter therefore
    writes a continuous polygon tube with shared rings and end caps only at
    the ends of each genuinely connected run.
    """
    radial_subdivisions = max(16, int(radial_subdivisions))
    helix_points = {}
    scale_nm = float(_oxdna.OXDNA_LENGTH_NM)
    com_radius = float(_oxdna.COM_RADIUS_OX)
    for record in records:
        strand = record.get("strand")
        if strand is None:
            continue
        # Auxiliary caDNAno channels remain distinct in every sequence and
        # mapping table, but they are routing detours rather than additional
        # physical duplexes.  Group their projected coordinates with the
        # corresponding ideal SST axis so the cylindrical model is exactly
        # the same compact geometry shown in the 3D preview.
        helix = int(record.get(
            "ideal_geometry_helix", strand.virtualHelix().number()))
        count = max(1, int(record.get("count", 1)))
        direction = int(record.get("direction", 1))
        parameter = (float(record.get("idx", 0)) + direction *
                     float(record.get("sub", 0)) / count)
        position = tuple(map(float, record["pos"]))
        a1 = tuple(map(float, record["a1"]))
        axis = tuple(
            (position[index] + com_radius * a1[index]) * scale_nm
            for index in range(3))
        key = round(parameter, 8)
        helix_points.setdefault(helix, {})[key] = (parameter, axis)

    def subtract(first, second):
        return tuple(first[index] - second[index] for index in range(3))

    def add(first, second):
        return tuple(first[index] + second[index] for index in range(3))

    def multiply(vector, scalar):
        return tuple(value * scalar for value in vector)

    def dot(first, second):
        return sum(first[index] * second[index] for index in range(3))

    def cross(first, second):
        return (
            first[1]*second[2] - first[2]*second[1],
            first[2]*second[0] - first[0]*second[2],
            first[0]*second[1] - first[1]*second[0],
        )

    def norm(vector):
        return math.sqrt(dot(vector, vector))

    def normalized(vector):
        length = norm(vector)
        if length <= 1e-12:
            return None
        return multiply(vector, 1.0/length)

    def line_distance(point, first, second):
        delta = subtract(second, first)
        length_squared = dot(delta, delta)
        if length_squared <= 1e-16:
            return norm(subtract(point, first))
        fraction = max(0.0, min(
            1.0, dot(subtract(point, first), delta)/length_squared))
        nearest = add(first, multiply(delta, fraction))
        return norm(subtract(point, nearest))

    def simplify(run, tolerance_nm=0.001):
        """Keep the twist path while collapsing exactly straight spans."""
        if len(run) <= 2:
            return run

        def recurse(first_index, last_index, retained):
            greatest = -1.0
            greatest_index = None
            first = run[first_index][1]
            last = run[last_index][1]
            for index in range(first_index+1, last_index):
                distance = line_distance(run[index][1], first, last)
                if distance > greatest:
                    greatest = distance
                    greatest_index = index
            if greatest_index is not None and greatest > tolerance_nm:
                recurse(first_index, greatest_index, retained)
                retained.add(greatest_index)
                recurse(greatest_index, last_index, retained)

        retained = {0, len(run)-1}
        # Preserve both rings at every preview-colour boundary. Otherwise a
        # perfectly straight Seed axis could simplify across Z1/Z2 or Z2/Z3
        # and place the colour transition at the wrong physical position.
        for index in range(1, len(run)):
            previous_color = _cylindrical_preview_color(
                name, run[index-1][0], seed_support_ranges)
            current_color = _cylindrical_preview_color(
                name, run[index][0], seed_support_ranges)
            if previous_color != current_color:
                retained.update((index-1, index))
        recurse(0, len(run)-1, retained)
        return [run[index] for index in sorted(retained)]

    def tangent(points, index):
        if index == 0:
            delta = subtract(points[1], points[0])
        elif index == len(points)-1:
            delta = subtract(points[-1], points[-2])
        else:
            delta = subtract(points[index+1], points[index-1])
        return normalized(delta)

    def initial_normal(axis_tangent):
        references = ((1.0, 0.0, 0.0),
                      (0.0, 1.0, 0.0),
                      (0.0, 0.0, 1.0))
        reference = min(references,
                        key=lambda candidate: abs(dot(
                            axis_tangent, candidate)))
        return normalized(cross(axis_tangent, reference))

    def tube_rings(points):
        tangents = [tangent(points, index) for index in range(len(points))]
        if any(item is None for item in tangents):
            return None
        normals = [initial_normal(tangents[0])]
        for axis_tangent in tangents[1:]:
            previous = normals[-1]
            transported = subtract(
                previous, multiply(axis_tangent,
                                   dot(previous, axis_tangent)))
            current = normalized(transported)
            if current is None:
                current = initial_normal(axis_tangent)
            normals.append(current)
        rings = []
        for center, axis_tangent, normal in zip(points, tangents, normals):
            binormal = normalized(cross(axis_tangent, normal))
            ring = []
            for division in range(radial_subdivisions):
                angle = 2.0*math.pi*division/radial_subdivisions
                offset = add(multiply(normal, math.cos(angle)),
                             multiply(binormal, math.sin(angle)))
                ring.append(add(center, multiply(offset, float(radius_nm))))
            rings.append(ring)
        return rings

    def polygon(vertices):
        coordinates = " ".join(
            "%.5f %.5f %.5f" % tuple(vertex) for vertex in vertices)
        return ".polygon " + coordinates

    lines = [
        ".comment DNA Moire Designer pure cylindrical model",
        ".comment Continuous tubes trace the exported helix axes; units "
        "are nm",
        ".comment Geometry, spacing, and twist match the paired PDB/mmCIF "
        "and "
        "oxDNA files",
        ".comment Shared cross-section rings remove base-to-base caps and "
        "seams",
        ".comment radial_subdivisions %d" % radial_subdivisions,
        ".comment preview_color Seed_Z1 %s" %
        CYLINDER_PREVIEW_COLORS["seed_z1"],
        ".comment preview_color Seed_Z2 %s" %
        CYLINDER_PREVIEW_COLORS["seed_z2"],
        ".comment preview_color Seed_Z3 %s" %
        CYLINDER_PREVIEW_COLORS["seed_z3"],
        ".comment preview_color SST_first_layer %s" %
        CYLINDER_PREVIEW_COLORS["sst_layer_1"],
        ".comment preview_color SST_second_layer %s" %
        CYLINDER_PREVIEW_COLORS["sst_layer_2"],
    ]
    active_color = None

    def select_color(parameter):
        nonlocal active_color
        color = _cylindrical_preview_color(
            name, parameter, seed_support_ranges)
        if color != active_color:
            lines.append(_bild_color_command(color))
            active_color = color

    axial_segment_count = 0
    polygon_count = 0
    tube_count = 0
    helix_count = 0
    for helix in sorted(helix_points):
        sampled = helix_points[helix]
        points = [item for unused_key, item in sorted(sampled.items())]
        if len(points) < 2:
            continue
        runs = []
        run = [points[0]]
        for item in points[1:]:
            # A one-index gap is a normal base step; a two-index gap can be
            # produced by one legal deletion.  Larger gaps represent a
            # genuinely disconnected axial region and must remain open.
            if abs(item[0] - run[-1][0]) > 2.01:
                if len(run) >= 2:
                    runs.append(run)
                run = [item]
            else:
                run.append(item)
        if len(run) >= 2:
            runs.append(run)
        if not runs:
            continue
        contour_length = 0.0
        local_segments = 0
        local_tubes = 0
        for connected_run in runs:
            simplified = simplify(connected_run)
            centers = [point for unused_parameter, point in simplified]
            rings = tube_rings(centers)
            if rings is None:
                continue
            for segment_index, (first_ring, second_ring) in enumerate(
                    zip(rings, rings[1:])):
                first_parameter = simplified[segment_index][0]
                second_parameter = simplified[segment_index+1][0]
                select_color((first_parameter + second_parameter) / 2.0)
                for division in range(radial_subdivisions):
                    following = (division+1) % radial_subdivisions
                    lines.append(polygon((
                        first_ring[division], first_ring[following],
                        second_ring[following], second_ring[division])))
                    polygon_count += 1
                local_segments += 1
            # Only the two physical ends of a connected tube are capped.
            select_color(simplified[0][0])
            lines.append(polygon(tuple(reversed(rings[0]))))
            select_color(simplified[-1][0])
            lines.append(polygon(tuple(rings[-1])))
            polygon_count += 2
            local_tubes += 1
            contour_length += sum(norm(subtract(second, first))
                                  for first, second in zip(
                                      centers, centers[1:]))
        if local_segments:
            helix_count += 1
            axial_segment_count += local_segments
            tube_count += local_tubes
            lines.append(
                ".comment helix %d actual_bases %d tube_spans %d "
                "mesh_axis_points %d contour_length_nm %.5f" %
                (helix, len(points), local_tubes,
                 local_segments+local_tubes, contour_length))
    lines.extend((".comment helix_count %d" % helix_count,
                  ".comment tube_spans %d" % tube_count,
                  ".comment mesh_axial_segments %d" % axial_segment_count,
                  ".comment surface_polygons %d" % polygon_count, ""))
    return "\n".join(lines)


def _add_cylindrical_model(report, records, root, name,
                           seed_support_ranges=None):
    """Add the BILD model beside an existing structure-bundle report."""
    path = Path(root) / (name + "_cylindrical_model.bild")
    _oxdna._atomic_write_many({
        str(path): _cylindrical_model_bild(
            records, name, seed_support_ranges=seed_support_ranges)})
    report.setdefault("paths", {})["cylindrical_model"] = str(path)
    report["cylindrical_model_path"] = str(path)
    return report


def _build_seed_sst_assembled_payload(source):
    """Materialize both 4x4 SST translations around the accepted 8x8 Seed.

    The accepted cadnano file intentionally contains only the origin A0/B0
    capture bridges.  Sequence export derives the A1/B1 translation in
    memory.  For the assembled atomistic model we need both physical copies:
    retain the accepted SST unit, duplicate the translated unit four lattice
    columns to its right, and apply the translated capture routing only to
    the second four helices on each Seed face.  The two face groups are
    disjoint, so this preserves the origin capture staples and yields all
    eight capture-face helices without changing the accepted cadnano JSON.
    """
    result = payload_to_internal_numbering(source)
    translated = payload_to_internal_numbering(
        runtime.build_translated_capture_export_payload(source))
    rows = {int(row["num"]): row for row in result["vstrands"]}
    translated_rows = {
        int(row["num"]): row for row in translated["vstrands"]}
    metadata = result.setdefault("moire_structure_metadata", {})
    layout = metadata.get("variable_length_layout", runtime.structure_layout())
    original_sst = sorted(map(int, metadata.get(
        "sst_helix_numbers", range(48, 64))))
    if len(original_sst) != 16:
        raise ValueError(
            "The assembled seed + SST sublattice model currently requires "
            "one 16-helix (4×4) SST sublattice unit.")
    translated_seed = sorted({
        int(bridge["seed_helix"])
        for assignment in runtime.capture_site_assignments(
            layout, "translated")
        for bridge in assignment.get("bridges", [])})
    for number in translated_seed:
        rows[number]["stap"] = deepcopy(translated_rows[number]["stap"])
        rows[number]["stap_colors"] = deepcopy(
            translated_rows[number].get("stap_colors", []))

    first_new = max(rows) + 1
    helix_map = {number: first_new + offset
                 for offset, number in enumerate(original_sst)}
    duplicate_rows = []
    for number in original_sst:
        row = deepcopy(translated_rows[number])
        row["num"] = helix_map[number]
        row["col"] = int(row["col"]) + 4
        for field in ("scaf", "stap"):
            for record in row.get(field, []):
                for slot in (0, 2):
                    if int(record[slot]) in helix_map:
                        record[slot] = helix_map[int(record[slot])]
        duplicate_rows.append(row)
    result["vstrands"].extend(duplicate_rows)

    # The translated Seed endpoints still name the source 4x4 unit; redirect
    # only those cross-interface partners to the duplicated unit.
    for number in translated_seed:
        for record in rows[number].get("stap", []):
            for slot in (0, 2):
                if int(record[slot]) in helix_map:
                    record[slot] = helix_map[int(record[slot])]

    for assignment in list(result.get("scaffold_sequences", [])):
        number = int(assignment.get("start_vh", -1))
        if number in helix_map:
            duplicate = deepcopy(assignment)
            duplicate["start_vh"] = helix_map[number]
            result["scaffold_sequences"].append(duplicate)
    result.pop("moire_staple_sequences", None)
    metadata.update({
        "stage": "seed_sst_assembled_structure",
        "sst_helix_numbers": original_sst + [
            helix_map[number] for number in original_sst],
        "assembled_sst_units_per_face": 2,
        "assembled_sst_unit_shape": "4x4 helices",
        "assembled_seed_face_coverage": "all 8 helices on both faces",
        "assembled_translation_rule":
            "origin 4x4 unit plus legal +4-column translated 4x4 unit",
        "assembled_duplicate_helix_map": {
            str(key): value for key, value in helix_map.items()},
    })
    result["name"] = "seed_sst_assembled.json"
    return result


def final_export(project_file, design, output_directory, workflow):
    project_payload = {}
    if project_file and Path(project_file).is_file():
        project_payload = _load(project_file)
        set_language(project_payload.get("settings", {}).get(
            "interface_language", "en"))
    source = _load(design)
    workflow_assignments = workflow.get("sequence_assignments", [])
    root = Path(output_directory).expanduser().resolve()
    if root.name != "final_export":
        root = root / "final_export"
    input_root = root / "Input parameters"
    design_root = root / "caDNAno design files"
    # Store the sequence deliverables directly in one top-level directory.
    # Do not create a nested SST/STT-input directory.
    sequence_root = root / "Oligonucleotide sequences"
    # U+2215 is a filesystem-safe slash on both macOS and Windows while
    # preserving the requested Finder/Explorer presentation.
    structure_root = root / "PDB∕oxView files"
    for folder in (input_root, design_root, sequence_root, structure_root):
        folder.mkdir(parents=True, exist_ok=True)
    if project_file and Path(project_file).is_file():
        shutil.copy2(project_file, input_root / Path(project_file).name)
    # Recent projects generate the complete SST file as part of the automatic
    # design and no longer require a separate SST-acceptance step.  Resolve
    # both schemas so final export always retains and uses that process file.
    automatic_exports = workflow.get("automatic_design_exports", {})
    complete_sst_path = (
        workflow.get("sst_accepted") or workflow.get("sst_review") or
        workflow.get("sst_two_layer") or automatic_exports.get("sst"))
    designs = {}
    for key, stem, path in (
            ("sst_accepted", "complete_sst", complete_sst_path),
            ("scaffold_accepted", "sst_scaffold_routing",
             workflow.get("scaffold_accepted")),
            ("structure_complete", "sst_scaffold_staple_capture",
             workflow.get("structure_complete"))):
        if path and Path(path).is_file():
            assignments = (
                _mapped_complete_sst_input_records(source, path)
                if key == "sst_accepted" else
                _mapped_input_records(
                    source, path,
                    ("seed_scaffold", "sst_input_layer_1",
                     "sst_input_layer_2")))
            designs[key] = _write_design_pair(
                path, design_root, stem, assignments)
    complete_sst_payload = None
    if complete_sst_path and Path(complete_sst_path).is_file():
        complete_sst_payload = _load(complete_sst_path)
        # The standalone process design owns the intact SST-output routing,
        # while the final accepted design owns the user's SST-input sequences.
        # Map those sequences back to the process design before reading its
        # complementary purchasing strands.
        complete_sst_payload["scaffold_sequences"] = \
            _mapped_complete_sst_input_records(source, complete_sst_path)
    input_groups, output_groups, manifest = _sequence_sheets(
        source, complete_sst_payload=complete_sst_payload)
    scaffold_summary = _scaffold_export_summary(
        source, workflow_assignments)
    workbook = sequence_root / "all_sequences.xlsx"
    project_settings = project_payload.get("settings", {})
    identical_layers = bool(
        project_settings.get("layers_design_sequence_identical", False) and
        project_settings.get("lattice_symmetry") != "square_kagome"
    ) if project_file and Path(project_file).is_file() else False
    input_l1 = "sst_input_L1 (=L2)" if identical_layers else "sst_input_L1"
    input_l2 = "sst_input_L2 (=L1)" if identical_layers else "sst_input_L2"
    output_l1 = "sst_output_L1 (=L2)" if identical_layers else "sst_output_L1"
    output_l2 = "sst_output_L2 (=L1)" if identical_layers else "sst_output_L2"
    _write_workbook(str(workbook), (
        ("scaffold", input_groups["scaffold"]),
        (input_l1, input_groups["sst_input_layer_1"]),
        (input_l2, input_groups["sst_input_layer_2"]),
        (output_l1, output_groups["sst_output_layer_1"]),
        (output_l2, output_groups["sst_output_layer_2"]),
        ("normal_staple", output_groups["normal_staple"]),
        ("staple_capture", output_groups["staple_capture"]),
    ), use_row_colors=True)
    _add_scaffold_metadata_columns(workbook, scaffold_summary)
    _add_capture_staple_columns(workbook, 7)
    localize_xlsx(workbook)
    _apply_capture_sequence_rich_text(workbook, 7, manifest)
    _write_json({
        "scaffold": [list(item) for item in input_groups["scaffold"]],
        "scaffold_summary": scaffold_summary,
        "sst_input_layer_1": [list(item) for item in
                              input_groups["sst_input_layer_1"]],
        "sst_input_layer_2": [list(item) for item in
                              input_groups["sst_input_layer_2"]],
        "sst_output_layer_1": [list(item) for item in
                               output_groups["sst_output_layer_1"]],
        "sst_output_layer_2": [list(item) for item in
                               output_groups["sst_output_layer_2"]],
        "normal_staple": [list(item) for item in
                          output_groups["normal_staple"]],
        "staple_capture": [list(item) for item in
                           output_groups["staple_capture"]],
        "capture_map": manifest,
    }, sequence_root / "all_sequences.json")

    internal = payload_to_internal_numbering(source)
    metadata = internal.get("moire_structure_metadata", {})
    layout = metadata.get("variable_length_layout", {})
    seed_helices = set(range(48))
    sst_helices = set(map(int, metadata.get(
        "sst_helix_numbers", range(48, 64))))
    complete_sst = payload_to_internal_numbering(
        build_complete_sst_only_payload(source, "complete_sst.json"))
    seed = _subset_payload(internal, seed_helices, name="seed.json")
    layer_payloads = []
    for layer_index, base_range in enumerate(
            layout.get("layer_ranges", [[48, 175], [208, 335]]), 1):
        # Include the 8-base staple overhang at each end of a complete layer.
        layer_payloads.append(_subset_payload(
            complete_sst, sst_helices,
            (int(base_range[0]) - 8, int(base_range[1]) + 8),
            "sst_layer_%d.json" % layer_index))
    seed_twist = (project_payload.get("prediction", {}).get(
        "reported_angle_deg", 0.0) if project_file and
        Path(project_file).is_file() else 0.0)
    structures = {
        "seed": _export_structure(
            seed, structure_root / "seed", "seed",
            seed_twist_angle_deg=seed_twist,
            seed_support_ranges=_seed_preview_support_ranges(layout)),
        "sst_layer_1": _export_structure(
            layer_payloads[0], structure_root / "sst_layer_1", "sst_layer_1"),
        "sst_layer_2": _export_structure(
            layer_payloads[1], structure_root / "sst_layer_2", "sst_layer_2"),
    }
    manifest_path = root / "manifest.json"
    manifest_payload = {
        "format": "DNA Moire Designer Final Export",
        "helix_center_spacing_nm": 2.8,
        "designs": designs,
        "sequence_workbook": str(workbook),
        "sequence_counts": {
            **{key: len(value) for key, value in input_groups.items()},
            **{key: len(value) for key, value in output_groups.items()}},
        "scaffolds": scaffold_summary,
        "structures": structures,
    }
    _write_json(manifest_payload, manifest_path)
    manifest_payload["root"] = str(root)
    manifest_payload["manifest"] = str(manifest_path)
    return manifest_payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=(
        "analyze", "extract-scaffold", "list-scaffolds",
        "assign-scaffold", "export-input-template",
        "import-input-template", "build-sequenced", "final-export"))
    parser.add_argument("arguments", nargs="*")
    args = parser.parse_args()
    values = args.arguments
    if len(values) == 2 and values[0] == "@arguments-file":
        values = json.loads(Path(values[1]).read_text(encoding="utf-8"))
        if not isinstance(values, list) or not all(
                isinstance(item, str) for item in values):
            raise ValueError(
                "The sequence-workflow arguments file is invalid.")
    if args.mode == "analyze":
        result = analyze(values[0])
    elif args.mode == "extract-scaffold":
        result = extract_scaffold(values[0], json.loads(values[1]))
    elif args.mode == "list-scaffolds":
        result = scaffold_catalog(
            int(values[0]), values[1] == "1", json.loads(values[2]))
    elif args.mode == "assign-scaffold":
        result = assign_standard_scaffold(
            json.loads(values[0]), values[1], values[2] == "1",
            json.loads(values[3]))
    elif args.mode == "export-input-template":
        result = export_template(values[0], values[1], values[2] == "1")
    elif args.mode == "import-input-template":
        result = import_template(values[0], values[1], values[2] == "1")
    elif args.mode == "build-sequenced":
        result = build_sequenced(values[0], values[1], json.loads(values[2]))
    else:
        result = final_export(
            values[0], values[1], values[2], json.loads(values[3]))
    # Worker stdout is a machine-readable transport channel.  Keep it ASCII
    # safe because frozen Windows GUI processes can inherit a legacy code page
    # such as cp1252; json.loads restores the original Unicode text exactly.
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()

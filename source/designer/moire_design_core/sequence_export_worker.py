#!/usr/bin/env python3
"""Create sequence and vector exports from isolated design states."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
import re
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from moire_runtime import worker_command

from moire_design_core import structure_worker as runtime
from moire_design_core.structure import (
    SST_LAYER_RANGES,
    build_complete_sst_only_payload,
    capture_column_color,
    capture_column_index,
    capture_pair_index,
)
from cadnano2.model.enum import LatticeType
from cadnano2.model.io.legacydecoder import import_legacy_dict
from cadnano2.model.io.sequencexlsx import (
    _write_workbook,
    write_sequence_workbook,
)
from moire_designer.i18n import localize_svg, localize_xlsx


def _oligo_nodes(oligo):
    return [
        (strand.virtualHelix().number(), strand.lowIdx(), strand.highIdx())
        for strand in oligo.strand5p().generator3pStrand()]


def _staple_base_sequences(payload):
    """Return saved/derived staple bases keyed by the full-design coordinate."""
    document = runtime.Document()
    part = import_legacy_dict(
        document, payload, LatticeType.Square, forceLatticeType=True)
    result = {}
    for oligo in part.oligos():
        if not oligo.isStaple() or oligo.strand5p() is None:
            continue
        for strand in oligo.strand5p().generator3pStrand():
            sequence = strand.sequence(forExport=True) or ""
            indices = (range(strand.lowIdx(), strand.highIdx() + 1)
                       if strand.isDrawn5to3() else
                       range(strand.highIdx(), strand.lowIdx() - 1, -1))
            for index, base in zip(indices, sequence):
                result[(strand.virtualHelix().number(), index)] = base
    return result


def _rows_for_variant(payload, variant, source_base_sequences=None,
                      source_uses_sst_first=True):
    document = runtime.Document()
    part = import_legacy_dict(
        document, payload, LatticeType.Square, forceLatticeType=True)
    input_rows = part.getInputSequenceRows()
    output_rows = []
    empty_sequences = 0
    layer_ranges = payload.get("moire_structure_metadata", {}).get(
        "variable_length_layout", {}).get(
            "layer_ranges", [list(item) for item in SST_LAYER_RANGES])
    metadata = payload.get("moire_structure_metadata", {})
    payload_uses_sst_first = metadata.get("helix_numbering") == "sst_first"
    default_sst = range(16) if payload_uses_sst_first else range(48, 64)
    default_seed = range(16, 64) if payload_uses_sst_first else range(48)
    sst_helices = set(map(int, metadata.get(
        "sst_helix_numbers", default_sst)))
    seed_helices = set(map(int, metadata.get(
        "seed_helix_numbers", default_seed)))
    for oligo in part.oligos():
        if oligo.strand5p() is None or not oligo.isStaple() or oligo.isLoop():
            continue
        nodes = _oligo_nodes(oligo)
        helices = {item[0] for item in nodes}
        if variant == "capture":
            selected = bool(helices & sst_helices and helices & seed_helices)
        elif variant == "complete_sst":
            # A Seed-free legacy JSON is renumbered 0..15 by cadnano's
            # importer; every staple in this isolated document is SST.
            selected = bool(nodes and all(
                any(low - 8 <= node_low and node_high <= high + 8
                    for low, high in layer_ranges)
                for unused_helix, node_low, node_high in nodes))
        else:
            selected = bool(nodes and min(helices) >= 48 and all(
                any(low - 8 <= node_low and node_high <= high + 8
                    for low, high in layer_ranges)
                for unused_helix, node_low, node_high in nodes))
        if not selected:
            continue
        record = list(oligo.sequenceRecord())
        if variant == "complete_sst" and source_base_sequences is not None:
            inherited = []
            for strand in oligo.strand5p().generator3pStrand():
                indices = (range(strand.lowIdx(), strand.highIdx() + 1)
                           if strand.isDrawn5to3() else
                           range(strand.highIdx(), strand.lowIdx() - 1, -1))
                inherited.extend(
                    source_base_sequences.get(
                        (strand.virtualHelix().number() +
                         (0 if source_uses_sst_first else 48), index), "?")
                    for index in indices)
            record[2] = "".join(inherited)
        if not record[2] or set(record[2]) <= {"?"}:
            empty_sequences += 1
        output_rows.append(tuple(record))
    from cadnano2.model.parts.part import _sequenceRowSortKey
    output_rows.sort(key=_sequenceRowSortKey)
    return input_rows, output_rows, empty_sequences


def _capture_sequence_manifest(payload):
    """Describe every exported Seed–SST bridge in structure coordinates."""
    document = runtime.Document()
    part = import_legacy_dict(
        document, payload, LatticeType.Square, forceLatticeType=True)
    metadata = payload.get("moire_structure_metadata", {})
    sst_first = metadata.get("helix_numbering") == "sst_first"
    default_sst = range(16) if sst_first else range(48, 64)
    default_seed = range(16, 64) if sst_first else range(48)
    sst_helices = set(map(int, metadata.get(
        "sst_helix_numbers", default_sst)))
    seed_helices = set(map(int, metadata.get(
        "seed_helix_numbers", default_seed)))
    assignments = metadata.get("capture_site_assignments_internal") or \
        metadata.get("variable_length_layout", {}).get(
            "capture_site_assignments", [])
    assignment_by_edge = {}
    for assignment in assignments:
        for bridge in assignment.get("bridges", []):
            assignment_by_edge[(
                int(bridge["seed_helix"]),
                int(bridge["sst_helix"]),
                int(assignment["position"]),
            )] = (assignment, bridge)
    rows_by_number = {
        int(row["num"]): row for row in payload.get("vstrands", [])}

    def internal_number(number, kind):
        if not sst_first:
            return int(number)
        return int(number) - 16 if kind == "seed" else int(number) + 48

    manifest = []
    seen = set()
    for oligo in part.oligos():
        if oligo.strand5p() is None or not oligo.isStaple() or oligo.isLoop():
            continue
        nodes = _oligo_nodes(oligo)
        helices = {item[0] for item in nodes}
        if not (helices & seed_helices and helices & sst_helices):
            continue
        record = list(oligo.sequenceRecord())
        for strand in oligo.strand5p().generator3pStrand():
            number = int(strand.virtualHelix().number())
            row = rows_by_number[number]
            for index in range(strand.lowIdx(), strand.highIdx() + 1):
                # Read the serialized row because it exposes both neighboring
                # endpoints without relying on private Strand methods.
                link = row["stap"][index]
                for offset in (0, 2):
                    partner = int(link[offset])
                    partner_index = int(link[offset + 1])
                    if partner < 0 or partner_index < 0:
                        continue
                    if number in seed_helices and partner in sst_helices:
                        seed_number, sst_number = number, partner
                        position = index
                    elif number in sst_helices and partner in seed_helices:
                        seed_number, sst_number = partner, number
                        position = partner_index
                    else:
                        continue
                    edge = (seed_number, sst_number, position)
                    if edge in seen:
                        continue
                    seen.add(edge)
                    seed_internal = internal_number(seed_number, "seed")
                    sst_internal = internal_number(sst_number, "sst")
                    assignment, bridge = assignment_by_edge.get(
                        (seed_internal, sst_internal, position), ({}, {}))
                    face = bridge.get(
                        "face", "face1" if seed_internal < 8 else "face2")
                    pair_index = (capture_pair_index(
                        position, metadata["variable_length_layout"])
                        if metadata.get("variable_length_layout") else None)
                    column_index = (capture_column_index(
                        position, metadata["variable_length_layout"])
                        if metadata.get("variable_length_layout") else None)
                    capture_color = (capture_column_color(
                        position, metadata["variable_length_layout"])
                        if metadata.get("variable_length_layout") else None)
                    if capture_color is None:
                        capture_color = 0x000000
                    sequence_runs = []
                    sequence_offset = 0
                    for component_strand in \
                            oligo.strand5p().generator3pStrand():
                        component_sequence = (
                            component_strand.sequence(forExport=True) or "")
                        component_length = len(component_sequence)
                        run_color = (capture_color if int(
                            component_strand.virtualHelix().number()) in
                            sst_helices else 0x000000)
                        run_role = ("capture_extension" if run_color else
                                    "staple_core")
                        if component_length:
                            sequence_runs.append({
                                "start": sequence_offset,
                                "end": sequence_offset + component_length,
                                "role": run_role,
                                "color": "#%06x" % run_color,
                            })
                        sequence_offset += component_length
                    manifest.append({
                        "face": face,
                        "seed_helix": seed_number,
                        "seed_internal_helix": seed_internal,
                        "sst_helix": sst_number,
                        "sst_internal_helix": sst_internal,
                        "base": position,
                        "capture_base": position,
                        "capture_seed_helix": seed_number,
                        "column_base": position,
                        "column_index": column_index,
                        "pair_index": pair_index,
                        "layer": assignment.get("layer"),
                        "phase": assignment.get("phase"),
                        "translation": assignment.get("translation"),
                        "cycle": assignment.get("cycle"),
                        "export_only": bool(
                            assignment.get("export_only", False)),
                        "connection_role": (
                            "export-only translated capture" if
                            assignment.get("export_only", False) else
                            "physical Seed-SST connection"),
                        "start": record[0],
                        "end": record[1],
                        "sequence": record[2],
                        "length": record[3],
                        "color": "#000000",
                        "staple_core_color": "#000000",
                        "capture_color": "#%06x" % capture_color,
                        "sequence_color_runs": sequence_runs,
                    })
    return _sort_capture_manifest_by_template_color(manifest)


def _sort_capture_manifest_by_template_color(manifest):
    """Sort Capture products by real-space column, then actual start.

    The function name is retained for API compatibility.  Colour is now a
    consequence of column identity and is never used as the sort key.
    """
    records = list(manifest)

    def sort_key(item):
        return (
            int(item.get("column_base", item["base"])),
            int(item.get("capture_base", item["base"])),
            int(item.get("capture_seed_helix", item["seed_helix"])),
            int(item.get("sst_helix", -1)), str(item.get("face", "")),
            bool(item.get("export_only")))

    return sorted(records, key=sort_key)


def _capture_map_start_end(item):
    seed_helix = int(item.get(
        "capture_seed_helix", item["seed_helix"]))
    capture_base = int(item.get("capture_base", item["base"]))
    if int(item.get("sst_helix", -1)) < 0:
        label = ("Potential Z2" if item.get("phase") == "Z2" else
                 "Candidate")
        return (
            "%s / Seed %d[%d]" % (
                label, seed_helix, capture_base),
            "No SST extension in current design",
        )
    phase = item.get("phase") or "?"
    # A/B remain distinct internal Square routing phases, but the exported
    # coordinate label only needs to identify the lattice family.  Present
    # both Square phases as S; retain K for Kagome.
    phase_label = "S" if phase in ("A", "B") else phase
    return (
        "%s / Seed %d[%d]" % (
            item["face"].replace("face", "Face "),
            seed_helix, capture_base),
        "SST %d[%d] / %s / %s / %s" % (
            item["sst_helix"], capture_base,
            phase_label, item.get("translation") or "?",
            "export-only" if item.get("export_only") else "physical"),
    )


def _capture_manifest_rows(manifest):
    return [(
        *_capture_map_start_end(item),
        item["sequence"], item["length"],
        item.get("color", item["capture_color"]))
        for item in manifest]


def _capture_staple_row(item):
    """Sequence row plus its corresponding capture-map coordinates."""
    capture_start, capture_end = _capture_map_start_end(item)
    return (
        item["start"], item["end"], item["sequence"], item["length"],
        item.get("color", item["capture_color"]),
        capture_start, capture_end)


def _add_capture_staple_columns(workbook, sheet_index):
    """Expose capture-map Start/End after Color on a sequence sheet."""
    workbook = Path(workbook)
    target_name = "xl/worksheets/sheet%d.xml" % int(sheet_index)
    temporary = workbook.with_name(workbook.name + ".capture-columns.tmp")
    with ZipFile(str(workbook), "r") as source, ZipFile(
            str(temporary), "w", ZIP_DEFLATED) as target:
        for info in source.infolist():
            data = source.read(info.filename)
            if info.filename == target_name:
                xml = data.decode("utf-8")
                xml = re.sub(
                    r'<dimension ref="A1:E(\d+)"\s*/>',
                    r'<dimension ref="A1:G\1"/>', xml, count=1)
                xml = re.sub(
                    r'<autoFilter ref="A1:E(\d+)"\s*/>',
                    r'<autoFilter ref="A1:G\1"/>', xml, count=1)
                xml = xml.replace(
                    "</cols>",
                    '<col min="6" max="6" width="30" customWidth="1"/>'
                    '<col min="7" max="7" width="52" customWidth="1"/>'
                    "</cols>", 1)
                header_cells = (
                    '<c r="F1" t="inlineStr" s="1"><is><t>'
                    'Capture Start</t></is></c>'
                    '<c r="G1" t="inlineStr" s="1"><is><t>'
                    'Capture End</t></is></c>')
                xml, count = re.subn(
                    r'(<row\b[^>]*\br="1"[^>]*>)(.*?)(</row>)',
                    lambda match: (
                        match.group(1) + match.group(2) + header_cells +
                        match.group(3)), xml, count=1)
                if count != 1:
                    raise RuntimeError(
                        "Could not add capture Start/End headers.")
                data = xml.encode("utf-8")
            target.writestr(info, data)
    temporary.replace(workbook)


def _apply_capture_sequence_rich_text(
        workbook, sheet_index, manifest):
    """Color only capture-extension letters in a generated XLSX sheet."""
    workbook = Path(workbook)
    if not manifest:
        return
    target_name = "xl/worksheets/sheet%d.xml" % int(sheet_index)
    temporary = workbook.with_name(workbook.name + ".capture-colors.tmp")
    with ZipFile(str(workbook), "r") as source, ZipFile(
            str(temporary), "w", ZIP_DEFLATED) as target:
        for info in source.infolist():
            data = source.read(info.filename)
            if info.filename == target_name:
                xml = data.decode("utf-8")
                for row_index, item in enumerate(manifest, 2):
                    sequence = str(item.get("sequence", ""))
                    rich_runs = []
                    for run in item.get("sequence_color_runs", ()):
                        start = max(0, int(run.get("start", 0)))
                        end = min(len(sequence), int(run.get("end", 0)))
                        if end <= start:
                            continue
                        color = str(run.get("color", "#000000")).lstrip(
                            "#").upper()
                        rich_runs.append(
                            '<r><rPr><rFont val="Menlo"/><color rgb="FF%s"/>'
                            '</rPr><t>%s</t></r>' % (
                                color, escape(sequence[start:end])))
                    if not rich_runs:
                        continue
                    reference = "C%d" % row_index
                    replacement = (
                        '<c r="%s" t="inlineStr"><is>%s</is></c>' %
                        (reference, "".join(rich_runs)))
                    xml, count = re.subn(
                        r'<c r="%s"[^>]*>.*?</c>' % reference,
                        replacement, xml, count=1)
                    if count != 1:
                        raise RuntimeError(
                            "Could not apply capture-segment color to %s." %
                            reference)
                data = xml.encode("utf-8")
            target.writestr(info, data)
    temporary.replace(workbook)


def _write_json(payload, path):
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                    encoding="utf-8")


def _render_svg(json_path, svg_path):
    completed = subprocess.run(
        worker_command("vector-export", str(json_path), str(svg_path)),
        check=False, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError("SVG生成失败：%s" % detail)
    localize_svg(svg_path)


def export_variants(source_path, output_directory, base_name=None):
    source_path = Path(source_path).expanduser().resolve()
    root = Path(output_directory).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if not source.get("scaffold_sequences"):
        raise ValueError(
            "所选JSON没有保存scaffold序列。请在cadnano中使用"
            "Save as with Sequences后再载入。")
    base = base_name or source_path.stem
    variants = {
        "capture": json.loads(json.dumps(source)),
        "complete_sst": build_complete_sst_only_payload(
            source, "%s_complete_sst_design.json" % base),
    }
    report = {"source": str(source_path), "source_unchanged": True,
              "variants": {}}
    source_base_sequences = _staple_base_sequences(source)
    source_uses_sst_first = source.get(
        "moire_structure_metadata", {}).get(
            "helix_numbering") == "sst_first"
    for variant, payload in variants.items():
        payload["name"] = "%s_%s_design.json" % (base, variant)
        metadata = payload.setdefault("moire_structure_metadata", {})
        metadata["export_role"] = variant
        json_path = root / payload["name"]
        xlsx_path = root / ("%s_%s_sequences.xlsx" % (base, variant))
        svg_path = root / ("%s_%s_design.svg" % (base, variant))
        input_rows, output_rows, empty_count = _rows_for_variant(
            payload, variant, source_base_sequences,
            source_uses_sst_first=source_uses_sst_first)
        if not output_rows:
            raise ValueError("%s状态没有找到可导出的序列链。" % variant)
        capture_manifest = []
        if variant == "capture":
            physical_manifest = _capture_sequence_manifest(payload)
            translated_payload = \
                runtime.build_translated_capture_export_payload(payload)
            unused_inputs, translated_rows, translated_empty = \
                _rows_for_variant(
                    translated_payload, variant, source_base_sequences,
                    source_uses_sst_first=source_uses_sst_first)
            # Origin and translated designs have disjoint Seed start/end
            # coordinates.  Keep one record per designed capture oligo while
            # guarding against malformed duplicate topology.
            combined = {}
            for row in list(output_rows) + list(translated_rows):
                combined[(str(row[0]), str(row[1]))] = tuple(row)
            output_rows = list(combined.values())
            from cadnano2.model.parts.part import _sequenceRowSortKey
            output_rows.sort(key=_sequenceRowSortKey)
            empty_count += translated_empty
            translated_manifest = _capture_sequence_manifest(
                translated_payload)
            capture_manifest = _sort_capture_manifest_by_template_color(
                physical_manifest + translated_manifest)
            capture_seen = set()
            output_rows = []
            ordered_manifest = []
            for item in capture_manifest:
                identity = (item["start"], item["end"], item["sequence"])
                if identity in capture_seen:
                    continue
                capture_seen.add(identity)
                ordered_manifest.append(item)
                output_rows.append(_capture_staple_row(item))
            capture_manifest = ordered_manifest
            expected_assignments = metadata.get(
                "capture_export_site_assignments_internal") or \
                metadata.get("variable_length_layout", {}).get(
                    "capture_export_site_assignments", [])
            expected_edges = {
                (int(bridge["seed_helix"]), int(bridge["sst_helix"]),
                 int(assignment["position"]))
                for assignment in expected_assignments
                for bridge in assignment.get("bridges", [])}
            manifest_edges = {
                (int(item["seed_internal_helix"]),
                 int(item["sst_internal_helix"]), int(item["base"]))
                for item in capture_manifest}
            if manifest_edges != expected_edges:
                raise ValueError(
                    "Capture序列平移不完整：缺少%d个位点，出现%d个意外位点。" %
                    (len(expected_edges - manifest_edges),
                     len(manifest_edges - expected_edges)))
            if len(output_rows) != len(capture_manifest):
                raise ValueError(
                    "Capture序列链数量与位点清单不一致：%d条序列 / %d个位点。" %
                    (len(output_rows), len(capture_manifest)))
            payload["capture_sequence_manifest"] = capture_manifest
            face_counts = {}
            for item in capture_manifest:
                face_counts.setdefault(item["face"], set()).add(
                    item["seed_helix"])
            metadata["capture_export_face_helix_counts"] = {
                face: len(helices) for face, helices in face_counts.items()}
            metadata["capture_physical_manifest_count"] = len(
                physical_manifest)
            metadata["capture_export_only_manifest_count"] = len(
                translated_manifest)
            metadata["capture_sequence_policy"] = (
                "A0/B0 are physical Seed-SST connections; A1/B1 are "
                "translated alternatives added only to sequence exports")
        unresolved_count = sum(str(row[2]).count("?") for row in output_rows)
        if variant == "complete_sst":
            payload["moire_staple_sequences"] = [
                {"start": row[0], "end": row[1], "sequence": row[2],
                 "length": row[3], "color": row[4]}
                for row in output_rows]
            payload["moire_structure_metadata"][
                "unresolved_single_strand_bases"] = unresolved_count
        _write_json(payload, json_path)
        if variant == "capture":
            _write_workbook(
                str(xlsx_path),
                (("input", input_rows), ("output", output_rows),
                 ("capture_map", _capture_manifest_rows(capture_manifest))),
                use_row_colors=True)
            _add_capture_staple_columns(xlsx_path, 2)
        else:
            write_sequence_workbook(str(xlsx_path), input_rows, output_rows)
        localize_xlsx(xlsx_path)
        if variant == "capture":
            _apply_capture_sequence_rich_text(
                xlsx_path, 2, capture_manifest)
        _render_svg(json_path, svg_path)
        report["variants"][variant] = {
            "json": str(json_path),
            "xlsx": str(xlsx_path),
            "svg": str(svg_path),
            "sequence_count": len(output_rows),
            "empty_sequence_count": empty_count,
            "unresolved_base_count": unresolved_count,
            "capture_manifest_count": len(capture_manifest),
        }
    print(json.dumps(report, ensure_ascii=True))
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output_directory")
    parser.add_argument("--name")
    arguments = parser.parse_args()
    export_variants(arguments.source, arguments.output_directory,
                    arguments.name)


if __name__ == "__main__":
    main()

"""Export cadnano DNA as an all-atom PDB/mmCIF and oxDNA input files."""

import csv
import colorsys
import io
import json
import math
import os
import tempfile

from .athena import decode_geometry_payload


OXDNA_LENGTH_NM = 0.8518
AXIAL_RISE_NM = 0.34
COM_RADIUS_OX = 0.60
BACKBONE_A1 = -0.34
BACKBONE_A2 = 0.3408
BASE_A1 = 0.40
BASES = set("ACGT")
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _vadd(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def _vsub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def _vscale(a, value):
    return tuple(value * a[i] for i in range(3))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(sum(value * value for value in a))


def _oligo_sort_key(oligo):
    strand = oligo.strand5p()
    if strand is None:
        return (9, 10 ** 9, 10 ** 9)
    part = strand.part()
    lattice = 0 if getattr(part, "_step", None) == 21 else 1
    return (lattice, strand.virtualHelix().number(), strand.idx5Prime(),
            1 if oligo.isStaple() else 0)


def _part_label(part):
    return "honeycomb" if getattr(part, "_step", None) == 21 else "square"


def _iter_strand_nucleotides(strand):
    """Yield (design index, insertion sub-index, sequence character)."""
    idx5 = strand.idx5Prime()
    idx3 = strand.idx3Prime()
    direction = 1 if idx3 >= idx5 else -1
    insertions = {item.idx(): item.length()
                  for item in strand.insertionsOnStrand()}
    sequence = strand.sequence() or ""
    seq_pos = 0
    for idx in range(idx5, idx3 + direction, direction):
        insertion = insertions.get(idx, 0)
        count = 0 if insertion == -1 else 1 + max(0, insertion)
        for sub in range(count):
            char = sequence[seq_pos].upper() if seq_pos < len(sequence) else "?"
            seq_pos += 1
            yield idx, sub, count, direction, char


def _collect(document, spacing_nm):
    oligos = sorted(document.oligos(), key=_oligo_sort_key)
    records = []
    strands = []
    by_strand = {}
    part_records = {}

    for oligo in oligos:
        first = oligo.strand5p()
        if first is None:
            continue
        oligo_records = []
        for strand in first.generator3pStrand():
            segment_records = []
            part = strand.part()
            vh = strand.virtualHelix()
            row, col = vh.coord()
            x_model, y_model = part.latticeCoordToPositionXY(row, col)
            scale = spacing_nm / (2.0 * float(part.radius())) / OXDNA_LENGTH_NM
            axis_x, axis_y = x_model * scale, y_model * scale
            is_staple = strand.isStaple()
            for idx, sub, count, direction, char in _iter_strand_nucleotides(strand):
                effective_idx = idx + direction * (float(sub) / max(1, count))
                theta = math.radians(getattr(part, "_twistOffset", 0.0) +
                                     effective_idx * part._twistPerBase)
                radial = (math.cos(theta), math.sin(theta), 0.0)
                if is_staple:
                    radial = _vscale(radial, -1.0)
                pos = (axis_x + COM_RADIUS_OX * radial[0],
                       axis_y + COM_RADIUS_OX * radial[1],
                       effective_idx * AXIAL_RISE_NM / OXDNA_LENGTH_NM)
                # a1 points from the base-pair centre towards the base site.
                a1 = _vscale(radial, -1.0)
                a3 = (0.0, 0.0, float(direction))
                a2 = _cross(a3, a1)
                rec = {"part": part, "strand": strand, "oligo": oligo,
                       "is_staple": is_staple, "idx": idx, "sub": sub,
                       "count": count, "direction": direction, "base": char,
                       "pos": pos, "a1": a1, "a2": a2, "a3": a3}
                records.append(rec)
                oligo_records.append(rec)
                segment_records.append(rec)
                part_records.setdefault(part, []).append(rec)
            by_strand[strand] = segment_records
        if oligo_records:
            strands.append({"oligo": oligo, "records": oligo_records,
                            "loop": bool(oligo.isLoop())})

    # Centre every lattice independently before hybrid fitting/placement.
    for part, recs in part_records.items():
        axes = []
        seen = set()
        for rec in recs:
            vh = rec["strand"].virtualHelix()
            key = vh.coord()
            if key not in seen:
                seen.add(key)
                row, col = key
                x, y = part.latticeCoordToPositionXY(row, col)
                scale = spacing_nm / (2.0 * float(part.radius())) / OXDNA_LENGTH_NM
                axes.append((x * scale, y * scale))
        cx = sum(x for x, _ in axes) / max(1, len(axes))
        cy = sum(y for _, y in axes) / max(1, len(axes))
        for rec in recs:
            rec["pos"] = _vsub(rec["pos"], (cx, cy, 0.0))

    hybrid_residual = _place_parts(part_records, by_strand, spacing_nm)
    assigned = _resolve_bases(records)
    return records, strands, assigned, hybrid_residual


def _place_parts(part_records, by_strand, spacing_nm):
    """Place Hybrid lattices exactly like the interactive 3D view.

    Cross-lattice endpoints cannot in general all be satisfied by one rigid
    translation: averaging those endpoint translations can put the complete
    honeycomb and square bodies on top of each other.  The 3D view instead
    lays each lattice out independently from its minimum X/Y axis and places
    the next lattice after the previous width plus a six-radius gap.  Apply
    the same deterministic layout here so PDB/mmCIF and oxDNA agree with the
    on-screen structure and can never interpenetrate merely because links
    occur at different positions.
    """
    parts = list(part_records)
    if len(parts) < 2:
        return 0.0
    residuals = []
    lattice_offset_x = 0.0
    for part in parts:
        axes = []
        seen = set()
        for rec in part_records[part]:
            vh = rec["strand"].virtualHelix()
            if vh.coord() in seen:
                continue
            seen.add(vh.coord())
            axes.append(part.latticeCoordToPositionXY(*vh.coord()))
        if not axes:
            continue
        scale = (spacing_nm / (2.0 * float(part.radius())) /
                 OXDNA_LENGTH_NM)
        scaled_axes = [(x * scale, y * scale) for x, y in axes]
        min_x = min(x for x, unused_y in scaled_axes)
        max_x = max(x for x, unused_y in scaled_axes)
        min_y = min(y for unused_x, y in scaled_axes)
        center_x = sum(x for x, unused_y in scaled_axes) / len(scaled_axes)
        center_y = sum(y for unused_x, y in scaled_axes) / len(scaled_axes)

        # Records were independently centred immediately before this call.
        # Shift the minimum helix axis to the same origin used by ThreeDView.
        centered_min_x = min_x - center_x
        centered_min_y = min_y - center_y
        translation = (lattice_offset_x - centered_min_x,
                       -centered_min_y, 0.0)
        for rec in part_records[part]:
            rec["pos"] = _vadd(rec["pos"], translation)

        minimum_width = 2.0 * float(part.radius()) * scale
        part_width = max(max_x - min_x, minimum_width)
        gap = max(float(part.radius()) * 6.0, 6.0) * scale
        lattice_offset_x += part_width + gap

    # Report the worst physical connection span offset after 3D placement.
    for strand, recs in by_strand.items():
        connected = strand.connection3p()
        target_recs = by_strand.get(connected, [])
        if connected is None or not recs or not target_recs:
            continue
        if connected.part() is strand.part():
            continue
        residuals.append(abs(_norm(_vsub(target_recs[0]["pos"], recs[-1]["pos"]))
                             - 0.76))
    return max(residuals) if residuals else 0.0


def _resolve_bases(records):
    paired = {}
    for rec in records:
        key = (rec["part"], rec["strand"].virtualHelix().coord(),
               rec["idx"], rec["sub"], rec["count"])
        paired.setdefault(key, {})[bool(rec["is_staple"])] = rec
    assigned = 0
    for pair in paired.values():
        scaffold = pair.get(False)
        staple = pair.get(True)
        scaf_known = scaffold is not None and scaffold["base"] in BASES
        stap_known = staple is not None and staple["base"] in BASES
        if scaffold is not None and not scaf_known and stap_known:
            scaffold["base"] = COMPLEMENT[staple["base"]]
            assigned += 1
            scaf_known = True
        if staple is not None and not stap_known and scaf_known:
            staple["base"] = COMPLEMENT[scaffold["base"]]
            assigned += 1
            stap_known = True
        if scaffold is not None and not scaf_known:
            scaffold["base"] = "A"
            assigned += 1
            scaf_known = True
        if staple is not None and not stap_known:
            staple["base"] = (COMPLEMENT[scaffold["base"]]
                              if scaffold is not None else "T")
            assigned += 1
    return assigned


def _number_records(strands):
    ordered = []
    for strand_id, strand_info in enumerate(strands, 1):
        # Classic oxDNA topology lists each strand in 3' -> 5' order.
        recs = list(reversed(strand_info["records"]))
        start = len(ordered)
        for rec in recs:
            rec["strand_id"] = strand_id
            rec["global_index"] = len(ordered)
            ordered.append(rec)
        for pos, rec in enumerate(recs):
            if strand_info["loop"]:
                rec["three_neighbor"] = start + ((pos - 1) % len(recs))
                rec["five_neighbor"] = start + ((pos + 1) % len(recs))
            else:
                rec["three_neighbor"] = -1 if pos == 0 else start + pos - 1
                rec["five_neighbor"] = (-1 if pos == len(recs) - 1
                                        else start + pos + 1)
    return ordered


def _top_text(ordered, strand_count):
    lines = ["%d %d" % (len(ordered), strand_count)]
    for rec in ordered:
        lines.append("%d %s %d %d" %
                     (rec["strand_id"], rec["base"],
                      rec["three_neighbor"], rec["five_neighbor"]))
    return "\n".join(lines) + "\n"


def _dat_text(ordered):
    if ordered:
        mins = [min(rec["pos"][i] for rec in ordered) for i in range(3)]
        maxs = [max(rec["pos"][i] for rec in ordered) for i in range(3)]
    else:
        mins, maxs = [0.0] * 3, [1.0] * 3
    margin = 5.0
    shift = tuple(margin - value for value in mins)
    box = [max(20.0, maxs[i] - mins[i] + 2.0 * margin)
           for i in range(3)]
    # Hybrid parts are deliberately displayed side by side.  Their physical
    # cross-lattice bonds can therefore be much longer than a normal
    # backbone bond.  If such a bond approaches half a periodic box, oxView
    # unwraps different sections of one circular scaffold into different
    # periodic images, making an intact lattice look torn in two even though
    # the raw coordinates (and PDB) are correct.  Keep every cross-part bond
    # component below one quarter of its box dimension to make imaging
    # unambiguous while preserving the actual .top connectivity.
    by_index = dict((rec["global_index"], rec) for rec in ordered)
    max_cross_component = [0.0, 0.0, 0.0]
    for rec in ordered:
        neighbor = by_index.get(rec.get("five_neighbor"))
        if neighbor is None or neighbor["part"] is rec["part"]:
            continue
        for axis in range(3):
            max_cross_component[axis] = max(
                max_cross_component[axis],
                abs(neighbor["pos"][axis] - rec["pos"][axis]))
    for axis in range(3):
        if max_cross_component[axis]:
            box[axis] = max(
                box[axis],
                4.0 * max_cross_component[axis] + 2.0 * margin)
    box = tuple(box)
    lines = ["t = 0", "b = %.8f %.8f %.8f" % box, "E = 0 0 0"]
    for rec in ordered:
        rec["output_pos"] = _vadd(rec["pos"], shift)
        values = (rec["output_pos"] + rec["a1"] + rec["a3"] +
                  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        lines.append(" ".join("%.8f" % value for value in values))
    return "\n".join(lines) + "\n"


CHAIN_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
MAPPING_FIELDS = (
    "mode", "lattice", "strand_type", "oligo_id", "chain_5to3",
    "base_order_5to3", "helix", "base_index", "insertion_index",
    "nucleotide", "oxdna_index", "pdb_chain", "pdb_residue_id",
    "atom_serial_start", "atom_serial_end", "center_x_angstrom",
    "center_y_angstrom", "center_z_angstrom")


def _structure_records(strands):
    """Yield records in conventional per-chain 5' -> 3' order."""
    for strand_id, strand_info in enumerate(strands, 1):
        recs = strand_info["records"]
        for residue_id, rec in enumerate(recs, 1):
            yield strand_id, residue_id, rec, strand_info


def _hy36(value, width):
    """Encode large PDB serial/residue identifiers using hybrid-36."""
    if value < 10 ** width:
        return ("%%%dd" % width) % value
    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    value = value - 10 ** width + 10 * 36 ** (width - 1)
    chars = []
    for unused in range(width):
        chars.append(digits[value % 36])
        value //= 36
    return "".join(reversed(chars))


def _pdb_atom_line(atom):
    return ("%-6s%s %-4s %-3s %1s%s    %8.3f%8.3f%8.3f"
            "  1.00  0.00          %2s" %
            (atom["group"], _hy36(atom["serial"], 5), atom["name"][:4],
             atom["resname"][:3], atom["pdb_chain"],
             _hy36(atom["residue_id"], 4), atom["xyz"][0], atom["xyz"][1],
             atom["xyz"][2], atom["element"][:2].rjust(2)))


def _pdb_text_from_atoms(atoms, remarks, connections=None):
    lines = ["REMARK 900 " + item for item in remarks]
    previous = None
    for atom in atoms:
        if previous is not None and atom["strand_id"] != previous:
            lines.append("TER")
        previous = atom["strand_id"]
        lines.append(_pdb_atom_line(atom))
    if atoms:
        lines.append("TER")
    for first, second in connections or []:
        lines.append("CONECT%s%s" % (_hy36(first, 5), _hy36(second, 5)))
    lines.append("END")
    return "\n".join(lines) + "\n"


def _cif_quote(value):
    value = str(value)
    if not value or any(char.isspace() for char in value) or \
            value[0] in "_#$;'\"":
        return "\"%s\"" % value.replace('"', "'")
    return value


def _mmcif_text(atoms, name, description):
    lines = ["data_%s" % ''.join(c if c.isalnum() else '_'
                                 for c in name),
             "_entry.id %s" % _cif_quote(name),
             "_struct.title %s" % _cif_quote(description), "#", "loop_"]
    columns = ("group_PDB", "id", "type_symbol", "label_atom_id",
               "label_comp_id", "label_asym_id", "label_entity_id",
               "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
               "occupancy", "B_iso_or_equiv", "auth_asym_id",
               "auth_seq_id", "pdbx_PDB_model_num")
    lines.extend("_atom_site." + column for column in columns)
    for atom in atoms:
        chain = "S%d" % atom["strand_id"]
        x, y, z = atom["xyz"]
        values = (atom["group"], atom["serial"], atom["element"],
                  atom["name"], atom["resname"], chain,
                  atom["strand_id"], atom["residue_id"], "%.5f" % x,
                  "%.5f" % y, "%.5f" % z, "1.00", "0.00", chain,
                  atom["residue_id"], 1)
        lines.append(" ".join(_cif_quote(value) for value in values))
    lines.extend(("#", ""))
    return "\n".join(lines)


def _mapping_row(mode, strand_id, residue_id, rec, start, end):
    center = tuple(value * OXDNA_LENGTH_NM * 10.0
                   for value in rec["output_pos"])
    return {"mode": mode, "lattice": _part_label(rec["part"]),
            "strand_type": "staple" if rec["is_staple"] else "scaffold",
            "oligo_id": strand_id, "chain_5to3": "S%d" % strand_id,
            "base_order_5to3": residue_id,
            "helix": rec["strand"].virtualHelix().number(),
            "base_index": rec["idx"], "insertion_index": rec["sub"],
            "nucleotide": rec["base"], "oxdna_index": rec["global_index"],
            "pdb_chain": CHAIN_CHARS[(strand_id - 1) % len(CHAIN_CHARS)],
            "pdb_residue_id": residue_id, "atom_serial_start": start,
            "atom_serial_end": end, "center_x_angstrom": "%.5f" % center[0],
            "center_y_angstrom": "%.5f" % center[1],
            "center_z_angstrom": "%.5f" % center[2]}


def _mapping_csv(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=MAPPING_FIELDS,
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _coarse_structure(strands, spacing_nm, hybrid_residual):
    """Build a one-backbone-bead PDB/mmCIF without a false centre helix."""
    atoms, mapping, connections = [], [], []
    serial = 1
    strand_serials = {}
    for strand_id, residue_id, rec, strand_info in _structure_records(strands):
        backbone = _vadd(rec["output_pos"],
                         _vadd(_vscale(rec["a1"], BACKBONE_A1),
                               _vscale(rec["a2"], BACKBONE_A2)))
        xyz = tuple(value * OXDNA_LENGTH_NM * 10.0 for value in backbone)
        atom = {"group": "HETATM", "serial": serial, "name": "BB",
                "resname": "D" + rec["base"], "strand_id": strand_id,
                "pdb_chain": CHAIN_CHARS[(strand_id - 1) % len(CHAIN_CHARS)],
                "residue_id": residue_id, "xyz": xyz, "element": "C"}
        atoms.append(atom)
        strand_serials.setdefault(strand_id, []).append(serial)
        mapping.append(_mapping_row("Coarse-grained", strand_id, residue_id,
                                    rec, serial, serial))
        serial += 1
    for strand_id, serials in strand_serials.items():
        for first, second in zip(serials, serials[1:]):
            connections.append((first, second))
        strand_info = strands[strand_id - 1]
        if strand_info["loop"] and len(serials) > 1:
            connections.append((serials[-1], serials[0]))
    remarks = ["CADNANO OXDNA-ALIGNED ONE-BEAD COARSE-GRAINED MODEL",
               "ONE OUTER BACKBONE BEAD PER NUCLEOTIDE; NO CENTRE BASE HELIX",
               "HELIX CENTER SPACING %.4f NM" % spacing_nm,
               "USE THE MATCHING TOP/DAT FILES FOR OXDNA SIMULATION",
               "INITIAL RELAXATION IS REQUIRED BEFORE PRODUCTION DYNAMICS"]
    if hybrid_residual:
        remarks.append("HYBRID 3D-LAYOUT MAX LINK SPAN OFFSET %.5f OXDNA UNITS" %
                       hybrid_residual)
    return atoms, mapping, remarks, connections


def _load_atom_templates():
    path = os.path.join(os.path.dirname(__file__), "dna_atom_templates.json")
    with open(path, encoding="utf-8") as source:
        return json.load(source)


def _all_atom_structure(strands, spacing_nm, hybrid_residual):
    templates = _load_atom_templates()
    atoms, mapping = [], []
    serial = 1
    for strand_id, residue_id, rec, strand_info in _structure_records(strands):
        if strand_info["loop"]:
            suffix = ""
        elif residue_id == 1:
            suffix = "5"
        elif residue_id == len(strand_info["records"]):
            suffix = "3"
        else:
            suffix = ""
        template = templates["D" + rec["base"] + suffix]
        start = serial
        origin = tuple(value * OXDNA_LENGTH_NM * 10.0
                       for value in rec["output_pos"])
        for atom_name, element, x, y, z in template:
            # Canonical template axes are x=a1, y=a2, z=a3.
            offset = _vadd(_vscale(rec["a1"], x),
                           _vadd(_vscale(rec["a2"], y),
                                 _vscale(rec["a3"], z)))
            xyz = _vadd(origin, offset)
            atoms.append({"group": "ATOM", "serial": serial,
                          "name": atom_name,
                          "resname": "D" + rec["base"] + suffix,
                          "strand_id": strand_id,
                          "pdb_chain": CHAIN_CHARS[
                              (strand_id - 1) % len(CHAIN_CHARS)],
                          "residue_id": residue_id, "xyz": xyz,
                          "element": element})
            serial += 1
        mapping.append(_mapping_row("All-atom", strand_id, residue_id, rec,
                                    start, serial - 1))
    remarks = ["CADNANO OXDNA-BACKMAPPED ALL-ATOM DNA MODEL",
               "HELIX CENTER SPACING %.4f NM" % spacing_nm,
               "SUGAR PHOSPHATE BASE AND HYDROGEN ATOMS ARE INCLUDED",
               "ENERGY MINIMIZATION IS REQUIRED BEFORE ALL-ATOM SIMULATION"]
    if hybrid_residual:
        remarks.append("HYBRID 3D-LAYOUT MAX LINK SPAN OFFSET %.5f OXDNA UNITS" %
                       hybrid_residual)
    return atoms, mapping, remarks


def _atomic_write_many(files):
    temporary = []
    try:
        for path, content in files.items():
            directory = os.path.dirname(path) or os.curdir
            fd, temp_path = tempfile.mkstemp(prefix=".cadnano-export-",
                                             dir=directory, text=True)
            temporary.append(temp_path)
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(content)
        for temp_path, path in zip(temporary, files):
            os.replace(temp_path, path)
        temporary = []
    finally:
        for path in temporary:
            try:
                os.unlink(path)
            except OSError:
                pass


def export_oxdna_bundle(document, pdb_path, spacing_nm=2.8):
    """Backward-compatible coarse export used by older callers/tests."""
    if spacing_nm <= 0:
        raise ValueError("Helix center spacing must be greater than zero.")
    root, extension = os.path.splitext(os.path.abspath(pdb_path))
    if extension.lower() != ".pdb":
        root = os.path.abspath(pdb_path)
    paths = {"pdb": root + ".pdb", "top": root + ".top", "dat": root + ".dat"}
    records, strands, assigned, residual = _collect(document, float(spacing_nm))
    if not records:
        raise ValueError("The design does not contain any DNA strands to export.")
    ordered = _number_records(strands)
    dat = _dat_text(ordered)
    atoms, unused_mapping, remarks, connections = _coarse_structure(
        strands, spacing_nm, residual)
    files = {paths["pdb"]: _pdb_text_from_atoms(
                                atoms, remarks, connections),
             paths["top"]: _top_text(ordered, len(strands)),
             paths["dat"]: dat}
    _atomic_write_many(files)
    return {"paths": paths, "nucleotides": len(ordered),
            "strands": len(strands), "assigned_bases": assigned,
            "hybrid_residual": residual}


def _safe_design_name(design_name):
    value = ''.join(char if char.isalnum() or char in "-_" else "_"
                    for char in str(design_name)).strip("_")
    return value or "cadnano-design"


def structure_bundle_paths(output_root, design_name, include_oxdna=True):
    """Return the possible paths for the compact structure export.

    Exactly one of ``pdb`` and ``cif`` is written.  The oxDNA pair is
    optional so Guided Design can let the user choose it explicitly.
    """
    root = os.path.abspath(output_root)
    design_name = _safe_design_name(design_name)
    paths = {
        "root": root,
        "all_atom": {
            "pdb": os.path.join(root, design_name + ".pdb"),
            "cif": os.path.join(root, design_name + ".cif")}}
    if include_oxdna:
        paths["oxdna"] = {
            "top": os.path.join(root, design_name + ".top"),
            "dat": os.path.join(root, design_name + ".dat")}
    return paths


def export_structure_bundle(document, output_root, design_name,
                            spacing_nm=2.8, include_oxdna=True,
                            pdb_atom_limit=99999):
    """Export one all-atom structure and, optionally, an oxDNA pair.

    PDB's conventional five-column atom serial field tops out at 99,999.
    Larger models are therefore written only as mmCIF.  Mapping tables and
    coarse-grained PDB/mmCIF files are deliberately not part of this compact
    export.
    """
    if spacing_nm <= 0:
        raise ValueError("Helix center spacing must be greater than zero.")
    design_name = _safe_design_name(design_name)
    paths = structure_bundle_paths(output_root, design_name, include_oxdna)
    records, strands, assigned, residual = _collect(
        document, float(spacing_nm))
    if not records:
        raise ValueError("The design does not contain any DNA strands to export.")
    ordered = _number_records(strands)
    # Besides producing the optional configuration text, this establishes the
    # shared shifted output coordinates used by the all-atom backmapping.
    dat = _dat_text(ordered)
    atom_atoms, unused_mapping, atom_remarks = _all_atom_structure(
        strands, spacing_nm, residual)

    os.makedirs(paths["root"], exist_ok=True)
    if len(atom_atoms) <= int(pdb_atom_limit):
        structure_format = "PDB"
        structure_path = paths["all_atom"]["pdb"]
        files = {structure_path: _pdb_text_from_atoms(
                                                atom_atoms, atom_remarks)}
    else:
        structure_format = "mmCIF"
        structure_path = paths["all_atom"]["cif"]
        files = {structure_path: _mmcif_text(
                    atom_atoms, design_name + "_all_atom",
                    "cadnano oxDNA-backmapped all-atom DNA model")}
    if include_oxdna:
        files[paths["oxdna"]["top"]] = _top_text(ordered, len(strands))
        files[paths["oxdna"]["dat"]] = dat
    _atomic_write_many(files)
    return {"paths": paths, "nucleotides": len(ordered),
            "strands": len(strands), "assigned_bases": assigned,
            "all_atom_count": len(atom_atoms),
            "structure_path": structure_path,
            "structure_format": structure_format,
            "include_oxdna": bool(include_oxdna),
            "hybrid_residual": residual}


def _apply_athena_frames(records, metadata):
    """Replace lattice coordinates with ATHENA target-space frames."""
    geometry = decode_geometry_payload(metadata)
    frames = geometry.get("frames", {})
    used = set()
    missing = []
    def frame_for(strand_type, vh_num, index):
        return frames.get("%s:%d:%d" % (strand_type, vh_num, index))

    def normalized_lerp(first, second, fraction):
        vector = tuple((1.0 - fraction) * float(first[i]) +
                       fraction * float(second[i]) for i in range(3))
        length = _norm(vector)
        if length <= 1e-12:
            return tuple(float(value) for value in first)
        return tuple(value / length for value in vector)

    for rec in records:
        strand_type = "staple" if rec["is_staple"] else "scaffold"
        vh_num = rec["strand"].virtualHelix().number()
        key = "%s:%d:%d" % (strand_type, vh_num, rec["idx"])
        frame = frames.get(key)
        if frame is None:
            missing.append(key)
            continue

        # Inserted nucleotides have no independent frame in the saved target
        # geometry.  Place them between this design position and the next
        # native position along the strand, preserving the curved path.  At a
        # circular boundary, use the opposite neighbouring native frame and
        # extrapolate in the strand direction rather than falling back to a
        # parallel-lattice coordinate.
        fraction = float(rec["sub"]) / max(1, int(rec["count"]))
        if rec["sub"]:
            next_index = rec["idx"] + rec["direction"]
            neighbour = frame_for(strand_type, vh_num, next_index)
            if neighbour is None:
                previous = frame_for(
                    strand_type, vh_num, rec["idx"] - rec["direction"])
                if previous is not None:
                    neighbour = {
                        "pos": _vadd(frame["pos"],
                                     _vsub(frame["pos"], previous["pos"])),
                        "a1": frame["a1"], "a3": frame["a3"]}
            if neighbour is None:
                missing.append(key + ":insertion")
                continue
            rec["pos"] = tuple(
                (1.0 - fraction) * float(frame["pos"][i]) +
                fraction * float(neighbour["pos"][i]) for i in range(3))
            rec["a1"] = normalized_lerp(
                frame["a1"], neighbour["a1"], fraction)
            rec["a3"] = normalized_lerp(
                frame["a3"], neighbour["a3"], fraction)
        else:
            rec["pos"] = tuple(float(value) for value in frame["pos"])
            rec["a1"] = tuple(float(value) for value in frame["a1"])
            rec["a3"] = tuple(float(value) for value in frame["a3"])
        rec["a2"] = _cross(rec["a3"], rec["a1"])
        used.add(key)
    if missing:
        examples = ", ".join(sorted(set(missing))[:8])
        raise ValueError(
            "The current design contains nucleotide positions that are not "
            "present in the saved ATHENA geometry: %s%s" %
            (examples, " ..." if len(set(missing)) > 8 else ""))
    return geometry, len(frames) - len(used)


def _bild_xyz(position):
    """Convert stored oxDNA coordinates to the Angstrom frame used by PDB."""
    scale = OXDNA_LENGTH_NM * 10.0
    return tuple(float(value) * scale for value in position)


def _bild_nm(position):
    """Convert stored oxDNA coordinates to the nm frame used by ATHENA BILD."""
    return tuple(float(value) * OXDNA_LENGTH_NM for value in position)


def _bild_color(index, saturation=0.58, value=0.86):
    """Return a deterministic, well-separated RGB color."""
    hue = (0.618033988749895 * (int(index) + 1)) % 1.0
    return colorsys.hsv_to_rgb(hue, saturation, value)


def _bild_color_line(color):
    return ".color %.5f %.5f %.5f" % tuple(color)


def _bild_cylinder_line(first, second, radius):
    return (".cylinder %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %
            (first[0], first[1], first[2],
             second[0], second[1], second[2], float(radius)))


def _bild_sphere_line(position, radius):
    return ".sphere %.5f %.5f %.5f %.5f" % (
        position[0], position[1], position[2], float(radius))


def _is_frame_geometry(metadata):
    """Return whether the embedded target is a locally bent Frame path."""
    return (str(metadata.get("format", "")).startswith("cadnano-frame-") or
            bool(metadata.get("frame_plan")))


def _curved_axis_geometry(metadata):
    """Return twist-free helix axes and optional fitted Curved rings in nm.

    A Curved Design axis is intentionally regularized to a circle.  A Frame
    axis must never pass through that fit: its straight edges and local bend
    windows are already encoded in ``geometry_data`` and are returned
    verbatim (after complementary-base centre averaging).
    """
    geometry = decode_geometry_payload(metadata)
    centres = {}
    for key, frame in geometry.get("frames", {}).items():
        pieces = str(key).split(":")
        if len(pieces) != 3:
            continue
        try:
            helix = int(pieces[1])
            index = int(pieces[2])
        except ValueError:
            continue
        centres.setdefault((helix, index), []).append(
            tuple(float(value) for value in frame["pos"]))

    by_helix = {}
    for (helix, index), positions in centres.items():
        count = float(len(positions))
        centre = tuple(sum(position[axis] for position in positions) / count
                       for axis in range(3))
        by_helix.setdefault(helix, {})[index] = _bild_nm(centre)

    if _is_frame_geometry(metadata):
        # New Frame files carry the exact twist-free duplex axes explicitly.
        # Use them instead of the phase-bearing nucleotide coordinates.
        stored_axes = (geometry.get("frame_geometry") or {}).get(
            "helix_axes") or {}
        exact = {}
        for key, position in stored_axes.items():
            pieces = str(key).split(":")
            if len(pieces) != 2:
                continue
            try:
                helix, index = int(pieces[0]), int(pieces[1])
            except ValueError:
                continue
            exact.setdefault(helix, {})[index] = _bild_nm(position)
        if exact:
            return exact, {}

        # Backward-compatible recovery for Frame JSON files written before
        # explicit axes were embedded.  Reconstruct the rounded-polygon path
        # from frame_plan and fit one constant cross-section offset per helix.
        # Averaging removes the stored duplex phase, so old exports also become
        # smooth rods rather than nucleotide-scale helices.
        plan = metadata.get("frame_plan") or {}
        try:
            from .frame import _point_on_path, _rounded_segments
            path = _rounded_segments(plan)
            nominal = float(plan["nominal_perimeter_bp"])
            origin_shift = float(plan.get("native_origin_shift_bp", 0.0))
        except (ImportError, KeyError, TypeError, ValueError):
            return by_helix, {}
        regularized = {}
        for helix, indexed in by_helix.items():
            samples = []
            for index in sorted(indexed):
                point, tangent = _point_on_path(
                    path, (index + 0.5 + origin_shift) / nominal)
                outward = (tangent[1], -tangent[0])
                raw = indexed[index]
                normal_offset = ((raw[0]-point[0])*outward[0] +
                                 (raw[1]-point[1])*outward[1])
                samples.append((index, point, outward, normal_offset,
                                raw[2]))
            if not samples:
                continue
            divisor = float(len(samples))
            normal_offset = sum(item[3] for item in samples) / divisor
            z_offset = sum(item[4] for item in samples) / divisor
            regularized[helix] = dict(
                (index, (point[0] + outward[0]*normal_offset,
                         point[1] + outward[1]*normal_offset, z_offset))
                for index, point, outward, unused_normal, unused_z in samples)
        return regularized or by_helix, {}

    fitted = {}
    rings = {}
    for helix, indexed in by_helix.items():
        indices = sorted(indexed)
        if not indices:
            continue
        points = [indexed[index] for index in indices]
        divisor = float(len(points))
        centre_x = sum(point[0] for point in points) / divisor
        centre_y = sum(point[1] for point in points) / divisor
        centre_z = sum(point[2] for point in points) / divisor
        radius = sum(math.hypot(point[0] - centre_x,
                                point[1] - centre_y)
                     for point in points) / divisor
        axis = {}
        for index, point in zip(indices, points):
            angle = math.atan2(point[1] - centre_y,
                               point[0] - centre_x)
            axis[index] = (centre_x + radius * math.cos(angle),
                           centre_y + radius * math.sin(angle), centre_z)
        fitted[helix] = axis
        rings[helix] = (centre_x, centre_y, centre_z, radius)
    return fitted, rings


def _axis_tangent(indexed, index):
    """Return the increasing-base tangent of a closed stored helix axis."""
    ordered = sorted(indexed)
    if len(ordered) < 2:
        return (1.0, 0.0, 0.0)
    try:
        position = ordered.index(index)
    except ValueError:
        position = min(range(len(ordered)),
                       key=lambda item: abs(ordered[item]-index))
    previous = indexed[ordered[(position-1) % len(ordered)]]
    following = indexed[ordered[(position+1) % len(ordered)]]
    vector = _vsub(following, previous)
    length = _norm(vector)
    if length <= 1e-12:
        following = indexed[ordered[(position+1) % len(ordered)]]
        current = indexed[ordered[position]]
        vector, length = _vsub(following, current), _norm(
            _vsub(following, current))
    return (_vscale(vector, 1.0/length) if length > 1e-12 else
            (1.0, 0.0, 0.0))


def _closed_axis_sampler(indexed, start_index):
    """Return an arclength sampler anchored at ``start_index``."""
    ordered = sorted(indexed)
    if not ordered:
        return None
    start_position = min(range(len(ordered)),
                         key=lambda item: abs(ordered[item]-start_index))
    ordered = ordered[start_position:] + ordered[:start_position]
    points = [indexed[index] for index in ordered]
    segments, total = [], 0.0
    for first, second in zip(points, points[1:] + points[:1]):
        vector = _vsub(second, first)
        length = _norm(vector)
        if length <= 1e-12:
            continue
        tangent = _vscale(vector, 1.0/length)
        segments.append((total, total+length, first, second, tangent))
        total += length
    if not segments or total <= 1e-12:
        return None

    def sample(fraction):
        distance = (float(fraction) % 1.0) * total
        for lower, upper, first, second, tangent in segments:
            if distance <= upper + 1e-12:
                local = (distance-lower) / max(1e-12, upper-lower)
                return (tuple(first[axis] +
                              (second[axis]-first[axis])*local
                              for axis in range(3)), tangent)
        unused_lower, unused_upper, first, unused_second, tangent = \
            segments[0]
        return first, tangent
    return sample


def _frame_region(metadata, native_coordinate):
    """Classify a nominal Frame coordinate as straight or bend."""
    plan = metadata.get("frame_plan") or {}
    nominal = float(plan.get("nominal_perimeter_bp", 0.0) or 0.0)
    centres = plan.get("vertex_native_centres") or ()
    lengths = plan.get("bend_length_bp") or ()
    if nominal <= 0.0:
        return "straight"
    coordinate = float(native_coordinate) % nominal
    for centre, length in zip(centres, lengths):
        distance = abs((coordinate-float(centre)+nominal/2.0) %
                       nominal-nominal/2.0)
        if distance <= float(length)/2.0 + 1e-9:
            return "bend"
    return "straight"


def _curved_axis_position(record, axes, rings):
    """Place a native or inserted nucleotide on its fitted curved axis."""
    redistributed = record.get("curved_axis_pos_nm")
    if redistributed is not None:
        return tuple(float(value) for value in redistributed)
    helix = record["strand"].virtualHelix().number()
    indexed = axes.get(helix, {})
    if not indexed:
        return _bild_nm(record["pos"])
    index = int(record["idx"])
    point = indexed.get(index)
    if point is None:
        nearest = min(indexed, key=lambda candidate: abs(candidate - index))
        point = indexed[nearest]
    if not record.get("sub"):
        return point
    neighbour_index = index + int(record.get("direction", 1))
    neighbour = indexed.get(neighbour_index)
    if neighbour is None:
        ordered = sorted(indexed)
        neighbour = indexed[ordered[0] if record.get("direction", 1) > 0
                            else ordered[-1]]
    fraction = float(record["sub"]) / max(1, int(record["count"]))
    ring = rings.get(helix)
    if ring is None:
        return tuple((1.0 - fraction) * point[axis] +
                     fraction * neighbour[axis] for axis in range(3))
    centre_x, centre_y, centre_z, radius = ring
    first_angle = math.atan2(point[1] - centre_y, point[0] - centre_x)
    second_angle = math.atan2(neighbour[1] - centre_y,
                              neighbour[0] - centre_x)
    delta = (second_angle - first_angle + math.pi) % (2.0 * math.pi) - math.pi
    angle = first_angle + fraction * delta
    return (centre_x + radius * math.cos(angle),
            centre_y + radius * math.sin(angle), centre_z)


def _rotate_about_z(vector, angle):
    cosine, sine = math.cos(angle), math.sin(angle)
    return (cosine * vector[0] - sine * vector[1],
            sine * vector[0] + cosine * vector[1], vector[2])


def _redistribute_frame_indels(records, metadata, axes):
    """Redistribute actual bases on the embedded rounded-polygon axes.

    This is the Frame counterpart of circular ring redistribution.  It uses
    closed-polyline arclength, so indels change the actual base spacing while
    straight edges remain straight and only the planned corner windows bend.
    """
    by_helix = {}
    for record in records:
        helix = record["strand"].virtualHelix().number()
        by_helix.setdefault(helix, []).append(record)
    nominal = float((metadata.get("frame_plan") or {}).get(
        "nominal_perimeter_bp", 0.0) or 0.0)
    summaries = {}
    for helix, helix_records in by_helix.items():
        indexed = axes.get(helix)
        if not indexed:
            continue
        count_by_index = {}
        for record in helix_records:
            index = int(record["idx"])
            count_by_index[index] = max(
                count_by_index.get(index, 0), int(record.get("count", 1)))
        occupied = sorted(index for index, count in count_by_index.items()
                          if count > 0)
        if not occupied:
            continue
        prefix, total = {}, 0
        for index in occupied:
            prefix[index] = total
            total += count_by_index[index]
        if total < 2:
            continue
        sampler = _closed_axis_sampler(indexed, occupied[0])
        if sampler is None:
            continue
        effective_nominal = nominal or float(len(indexed))
        for record in helix_records:
            index = int(record["idx"])
            count = max(1, int(record.get("count", 1)))
            sub = int(record.get("sub", 0))
            insertion_slot = (sub if record.get("direction", 1) > 0 or
                              sub == 0 else count-sub)
            slot = prefix[index] + insertion_slot
            new_axis_nm, new_tangent = sampler(slot/float(total))
            old_axis_nm = _curved_axis_position(record, axes, {})
            old_tangent = _axis_tangent(indexed, index)
            old_angle = math.atan2(old_tangent[1], old_tangent[0])
            new_angle = math.atan2(new_tangent[1], new_tangent[0])
            rotation = ((new_angle-old_angle+math.pi) %
                        (2.0*math.pi)-math.pi)
            old_axis_ox = tuple(value/OXDNA_LENGTH_NM
                                for value in old_axis_nm)
            new_axis_ox = tuple(value/OXDNA_LENGTH_NM
                                for value in new_axis_nm)
            relative = _vsub(record["pos"], old_axis_ox)
            record["pos"] = _vadd(
                new_axis_ox, _rotate_about_z(relative, rotation))
            record["a1"] = _rotate_about_z(record["a1"], rotation)
            record["a3"] = _rotate_about_z(record["a3"], rotation)
            record["a2"] = _cross(record["a3"], record["a1"])
            record["curved_axis_pos_nm"] = new_axis_nm
            record["curved_actual_slot"] = slot
            record["curved_actual_count"] = total
            record["frame_tangent_nm"] = new_tangent
            coordinate = ((slot/float(total))*effective_nominal +
                          occupied[0])
            record["frame_region"] = _frame_region(metadata, coordinate)
        summaries[helix] = total
    return summaries


def _redistribute_curved_indels(records, metadata):
    """Map each ring to uniform cumulative coordinates after all indels.

    DNAxiS creates the ideal ring with a complete native lattice period and
    cadnano subsequently encodes the circumference correction as insertions
    and deletions.  Keeping native angular coordinates would make an inserted
    base merely subdivide one old interval while leaving every crossover at
    its pre-indel angle.  This pass instead gives every *actual* base one
    angular slot and therefore moves crossover positions by the cumulative
    indel count.  The same transformed records feed DAT and PDB/mmCIF; BILD
    reads the stored axis positions, so all exported formats stay aligned.
    """
    axes, rings = _curved_axis_geometry(metadata)
    if _is_frame_geometry(metadata):
        return _redistribute_frame_indels(records, metadata, axes)
    by_helix = {}
    for record in records:
        helix = record["strand"].virtualHelix().number()
        by_helix.setdefault(helix, []).append(record)

    summaries = {}
    for helix, helix_records in by_helix.items():
        indexed = axes.get(helix)
        ring = rings.get(helix)
        if not indexed or ring is None:
            continue
        # One native index contributes either zero bases (deletion) or
        # 1 + insertion_count actual angular slots.  Both strands carry the
        # same edit, so use the maximum count observed at each index.
        count_by_index = {}
        for record in helix_records:
            index = int(record["idx"])
            count_by_index[index] = max(
                count_by_index.get(index, 0), int(record.get("count", 1)))
        occupied = sorted(index for index, count in count_by_index.items()
                          if count > 0)
        if not occupied:
            continue
        prefix = {}
        total = 0
        for index in occupied:
            prefix[index] = total
            total += count_by_index[index]
        if total < 2:
            continue

        first_index = occupied[0]
        centre_x, centre_y, centre_z, radius = ring
        first_point = indexed.get(first_index)
        if first_point is None:
            first_point = indexed[min(indexed,
                                      key=lambda value: abs(value-first_index))]
        phase = math.atan2(first_point[1] - centre_y,
                           first_point[0] - centre_x)
        direction = 1.0
        if len(occupied) > 1:
            second_point = indexed.get(occupied[1])
            if second_point is not None:
                second_angle = math.atan2(second_point[1] - centre_y,
                                          second_point[0] - centre_x)
                native_delta = ((second_angle - phase + math.pi) %
                                (2.0 * math.pi) - math.pi)
                direction = 1.0 if native_delta >= 0.0 else -1.0

        for record in helix_records:
            index = int(record["idx"])
            count = max(1, int(record.get("count", 1)))
            sub = int(record.get("sub", 0))
            # Complementary reverse-direction insertion records enumerate
            # the same physical inserted bases in reverse order.
            insertion_slot = (sub if record.get("direction", 1) > 0 or
                              sub == 0 else count - sub)
            slot = prefix[index] + insertion_slot
            new_angle = phase + direction * 2.0 * math.pi * slot / total
            new_axis_nm = (centre_x + radius * math.cos(new_angle),
                           centre_y + radius * math.sin(new_angle), centre_z)
            old_axis_nm = _curved_axis_position(record, axes, rings)
            old_angle = math.atan2(old_axis_nm[1] - centre_y,
                                   old_axis_nm[0] - centre_x)
            rotation = ((new_angle - old_angle + math.pi) %
                        (2.0 * math.pi) - math.pi)

            centre_ox = (centre_x / OXDNA_LENGTH_NM,
                         centre_y / OXDNA_LENGTH_NM,
                         centre_z / OXDNA_LENGTH_NM)
            relative = _vsub(record["pos"], centre_ox)
            record["pos"] = _vadd(centre_ox,
                                  _rotate_about_z(relative, rotation))
            record["a1"] = _rotate_about_z(record["a1"], rotation)
            record["a3"] = _rotate_about_z(record["a3"], rotation)
            record["a2"] = _cross(record["a3"], record["a1"])
            record["curved_axis_pos_nm"] = new_axis_nm
            record["curved_actual_slot"] = slot
            record["curved_actual_count"] = total
        summaries[helix] = total
    return summaries


def _curved_cylindrical_bild(metadata, strands):
    """Trace each duplex centre axis as a continuous 2-nm BILD rod."""
    axes, rings = _curved_axis_geometry(metadata)
    is_frame = _is_frame_geometry(metadata)
    scaffold_points = {}
    fallback_points = {}
    for strand_info in strands:
        for record in strand_info.get("records", []):
            helix = record["strand"].virtualHelix().number()
            fraction = (float(record.get("sub", 0)) /
                        max(1, int(record.get("count", 1))))
            parameter = float(record.get(
                "curved_actual_slot",
                float(record["idx"]) +
                record.get("direction", 1) * fraction))
            item = (parameter, _curved_axis_position(record, axes, rings),
                    record.get("frame_region", "straight"))
            fallback_points.setdefault(helix, {})[round(parameter, 8)] = item
            if not record["is_staple"]:
                scaffold_points.setdefault(helix, {})[
                    round(parameter, 8)] = item

    lines = [
        "# Curved / Frame Design cylindrical model",
        "# Smooth dsDNA centre-axis rods (no double-helix trace); units are nm",
        "# dsDNA_rod_diameter_nm 2.00000"]
    if is_frame:
        lines.extend(("# straight_region_color light gray (#748493)",
                      "# bend_region_color orange (#ff710d)"))
    else:
        lines.append(".color 0.45500 0.51800 0.57600")
    segment_count = 0
    for helix in sorted(axes):
        sampled = scaffold_points.get(helix) or fallback_points.get(helix, {})
        samples = [item for unused_key, item in sorted(sampled.items())]
        if len(samples) < 2:
            samples = [(float(index), axes[helix][index], "straight")
                       for index in sorted(axes[helix])]
        points = [item[1] for item in samples]
        if len(points) < 2:
            continue
        contour_length = 0.0
        following_samples = samples[1:] + samples[:1]
        for first_sample, second_sample in zip(samples, following_samples):
            first, second = first_sample[1], second_sample[1]
            if _norm(_vsub(first, second)) <= 1e-8:
                continue
            if is_frame:
                bend = (first_sample[2] == "bend" or
                        second_sample[2] == "bend")
                lines.append(
                    ".color 1.00000 0.44300 0.05100" if bend else
                    ".color 0.45500 0.51800 0.57600")
            lines.append(_bild_cylinder_line(first, second, 1.000))
            contour_length += _norm(_vsub(first, second))
            segment_count += 1
        lines.append("# helix %d actual_bases %d contour_length_nm %.5f" %
                     (helix, len(points), contour_length))
    lines.extend(("# helix_count %d" % len(axes),
                  "# cylinder_segments %d" % segment_count, ""))
    return "\n".join(lines)


def _curved_routing_bild(strands, metadata):
    """Render current routing on the saved Curved/Frame helix axes."""
    axes, rings = _curved_axis_geometry(metadata)
    is_frame = _is_frame_geometry(metadata)

    def routing_position(record):
        point = _curved_axis_position(record, axes, rings)
        if not record["is_staple"]:
            return point
        # ATHENA's routing_model_multi draws scaffold on the cylindrical
        # centreline and offsets staples by 0.5 nm so the complementary paths
        # remain separately visible.  Apply the same convention radially on
        # each Curved ring.
        helix = record["strand"].virtualHelix().number()
        ring = rings.get(helix)
        if ring is None:
            tangent = record.get("frame_tangent_nm")
            if tangent is None:
                tangent = _axis_tangent(axes.get(helix, {}),
                                        int(record["idx"]))
            normal = (-float(tangent[1]), float(tangent[0]), 0.0)
            return tuple(point[axis]-0.5*normal[axis]
                         for axis in range(3))
        centre_x, centre_y, unused_z, unused_radius = ring
        radial_x = point[0] - centre_x
        radial_y = point[1] - centre_y
        radial_length = math.hypot(radial_x, radial_y)
        if radial_length <= 1e-12:
            return point
        offset = 0.5 / radial_length
        return (point[0] - offset * radial_x,
                point[1] - offset * radial_y, point[2])

    lines = [
        "# Curved / Frame Design routing model",
        "# Wireframe-compatible routing on target helix axes; units are nm"]
    if is_frame:
        lines.extend(("# straight_region_color light gray",
                      "# bend_region_color orange"))
    segment_count = 0
    for strand_id, strand_info in enumerate(strands, 1):
        records = strand_info.get("records", [])
        if not records:
            continue
        is_staple = bool(records[0]["is_staple"])
        if is_frame:
            pass
        elif is_staple:
            lines.append(_bild_color_line(
                _bild_color(strand_id, 0.58, 0.86)))
        else:
            lines.append(".color steel blue")
        points = [routing_position(record) for record in records]
        pairs = list(zip(points, points[1:]))
        if strand_info.get("loop") and len(points) > 1:
            pairs.append((points[-1], points[0]))
        record_pairs = list(zip(records, records[1:]))
        if strand_info.get("loop") and len(records) > 1:
            record_pairs.append((records[-1], records[0]))
        for (first, second), (first_record, second_record) in zip(
                pairs, record_pairs):
            if _norm(_vsub(first, second)) <= 1e-8:
                continue
            if is_frame:
                bend = (first_record.get("frame_region") == "bend" or
                        second_record.get("frame_region") == "bend")
                lines.append(
                    ".color 1.00000 0.44300 0.05100" if bend else
                    ".color 0.45500 0.51800 0.57600")
            lines.append(_bild_cylinder_line(first, second, 0.100))
            segment_count += 1
        # Match Wireframe routing conventions: red 5' endpoint sphere and a
        # directional arrow at the 3' end.  Crossover endpoints get small
        # spheres naturally visible at the junction between helix axes.
        for before, after, before_record, after_record in zip(
                points, points[1:], records, records[1:]):
            before_helix = before_record["strand"].virtualHelix().number()
            after_helix = after_record["strand"].virtualHelix().number()
            if before_helix != after_helix:
                lines.append(_bild_sphere_line(before, 0.100))
                lines.append(_bild_sphere_line(after, 0.100))
        if not strand_info.get("loop"):
            lines.append(".color red")
            lines.append(_bild_sphere_line(points[0], 0.400))
            if is_frame:
                lines.append(
                    ".color 1.00000 0.44300 0.05100" if
                    records[-1].get("frame_region") == "bend" else
                    ".color 0.45500 0.51800 0.57600")
            else:
                lines.append(_bild_color_line(
                    _bild_color(strand_id, 0.58, 0.86)) if is_staple else
                    ".color steel blue")
            if len(points) > 1:
                start = tuple(0.70 * points[-2][axis] +
                              0.30 * points[-1][axis] for axis in range(3))
                vector = _vsub(points[-1], points[-2])
                end = _vadd(points[-1], tuple(0.10 * value
                                              for value in vector))
                lines.append(
                    ".arrow %.5f %.5f %.5f %.5f %.5f %.5f "
                    "0.100 0.300 0.300" %
                    (start[0], start[1], start[2],
                     end[0], end[1], end[2]))
    lines.extend(("# strand_count %d" % len(strands),
                  "# routing_segments %d" % segment_count, ""))
    return "\n".join(lines)


def athena_structure_bundle_paths(output_root, design_name):
    root = os.path.abspath(output_root)
    design_name = _safe_design_name(design_name)
    return {
        "root": root,
        "pdb": os.path.join(root, design_name + ".pdb"),
        "cif": os.path.join(root, design_name + ".cif"),
        "top": os.path.join(root, design_name + ".top"),
        "dat": os.path.join(root, design_name + ".dat"),
        "info": os.path.join(root, design_name + "_export_info.txt")}


def export_athena_structure_bundle(document, metadata, output_root,
                                   design_name, pdb_atom_limit=99999):
    """Export the current topology/sequence in saved ATHENA 3D coordinates."""
    if not metadata:
        raise ValueError("The current document is not an ATHENA design.")
    design_name = _safe_design_name(design_name)
    paths = athena_structure_bundle_paths(output_root, design_name)
    # Spacing is irrelevant because every record is replaced by a stored
    # target-space frame.  A conventional value keeps _collect deterministic.
    records, strands, assigned, unused_residual = _collect(document, 2.8)
    if not records:
        raise ValueError("The design contains no DNA strands.")
    geometry, removed_count = _apply_athena_frames(records, metadata)
    ordered = _number_records(strands)
    dat = _dat_text(ordered)
    atom_atoms, unused_mapping, atom_remarks = _all_atom_structure(
        strands, 2.8, 0.0)
    atom_remarks = [
        "ATHENA TARGET-GEOMETRY ALL-ATOM DNA MODEL",
        "SEQUENCE AND TOPOLOGY TAKEN FROM THE CURRENT CADNANO DOCUMENT",
        "UNRELAXED IDEAL INITIAL CONFIGURATION"] + atom_remarks[1:]

    os.makedirs(paths["root"], exist_ok=True)
    if len(atom_atoms) <= int(pdb_atom_limit):
        structure_format = "PDB"
        structure_path = paths["pdb"]
        files = {
            structure_path: _pdb_text_from_atoms(atom_atoms, atom_remarks)}
    else:
        structure_format = "mmCIF"
        structure_path = paths["cif"]
        files = {
            structure_path: _mmcif_text(
                atom_atoms, design_name,
                "ATHENA target-geometry unrelaxed all-atom DNA model")}
    files[paths["top"]] = _top_text(ordered, len(strands))
    files[paths["dat"]] = dat
    files[paths["info"]] = (
        "Wireframe 3D Export\n"
        "================\n"
        "Engine: %s\n"
        "Dimension: %s\n"
        "Edge type: %s\n"
        "Edge length: %s bp\n"
        "Nucleotides: %d\n"
        "Strands: %d\n"
        "All-atom format: %s\n"
        "All-atom count: %d\n"
        "Coordinate source: %s\n"
        "Relaxation: not performed\n"
        "Note: PDB/mmCIF and DAT preserve the ideal ATHENA target geometry; "
        "they are not equilibrium conformations.\n" %
        (metadata.get("engine", "unknown"),
         metadata.get("dimension", "unknown"),
         metadata.get("edge_type", "unknown"),
         metadata.get("edge_length_bp", "unknown"),
         len(ordered), len(strands), structure_format, len(atom_atoms),
         geometry.get("source", "ATHENA")))
    _atomic_write_many(files)
    return {
        "paths": paths,
        "structure_path": structure_path,
        "structure_format": structure_format,
        "all_atom_count": len(atom_atoms),
        "nucleotides": len(ordered),
        "strands": len(strands),
        "assigned_bases": assigned,
        "removed_geometry_bases": removed_count}


def curved_structure_bundle_paths(output_root, design_name):
    """Return the final sequence-accurate Curved / Frame export paths."""
    paths = athena_structure_bundle_paths(output_root, design_name)
    paths.update({
        "cylindrical_model": os.path.join(
            paths["root"], design_name + "_cylindrical_model.bild"),
        "routing_model_multi": os.path.join(
            paths["root"], design_name + "_routing_model_multi.bild")})
    return paths


def export_curved_structure_bundle(document, metadata, output_root,
                                   design_name, pdb_atom_limit=99999):
    """Export current sequence/topology in its saved Curved/Frame geometry."""
    if not metadata:
        raise ValueError("The current document is not a Curved Design.")
    design_name = _safe_design_name(design_name)
    paths = curved_structure_bundle_paths(output_root, design_name)
    records, strands, assigned, unused_residual = _collect(document, 2.8)
    if not records:
        raise ValueError("The design contains no DNA strands.")
    geometry, removed_count = _apply_athena_frames(records, metadata)
    actual_ring_counts = _redistribute_curved_indels(records, metadata)
    ordered = _number_records(strands)
    dat = _dat_text(ordered)
    atom_atoms, unused_mapping, atom_remarks = _all_atom_structure(
        strands, 2.8, 0.0)
    frame_export = _is_frame_geometry(metadata)
    target_label = ("FRAME ROUNDED-POLYGON" if frame_export else
                    "DNAXIS CURVED")
    atom_remarks = [
        "%s TARGET-GEOMETRY ALL-ATOM DNA MODEL" % target_label,
        "SEQUENCE AND TOPOLOGY TAKEN FROM THE CURRENT CADNANO DOCUMENT",
        "UNRELAXED IDEAL INITIAL CONFIGURATION"] + atom_remarks[1:]

    os.makedirs(paths["root"], exist_ok=True)
    if len(atom_atoms) <= int(pdb_atom_limit):
        structure_format = "PDB"
        structure_path = paths["pdb"]
        files = {
            structure_path: _pdb_text_from_atoms(atom_atoms, atom_remarks)}
    else:
        structure_format = "mmCIF"
        structure_path = paths["cif"]
        files = {
            structure_path: _mmcif_text(
                atom_atoms, design_name,
                "%s target-geometry unrelaxed all-atom DNA model" %
                target_label)}
    files[paths["top"]] = _top_text(ordered, len(strands))
    files[paths["dat"]] = dat
    files[paths["cylindrical_model"]] = \
        _curved_cylindrical_bild(metadata, strands)
    files[paths["routing_model_multi"]] = \
        _curved_routing_bild(strands, metadata)
    files[paths["info"]] = (
        "Curved / Frame 3D Export\n"
        "================\n"
        "Shape: %s\n"
        "Height: %s nm\n"
        "Maximum diameter: %s nm\n"
        "Minimum diameter: %s nm\n"
        "Layers: %s\n"
        "Nucleotides: %d\n"
        "Strands: %d\n"
        "All-atom format: %s\n"
        "All-atom count: %d\n"
        "Coordinate source: %s\n"
        "Indel-aware target-path coordinates: yes (%d helices)\n"
        "Target-path type: %s\n"
        "Cylindrical model: %s\n"
        "Routing model (multi-color): %s\n"
        "Relaxation: not performed\n"
        "Note: PDB/mmCIF and DAT preserve the ideal Curved or Frame "
        "target geometry; they are not parallel-lattice coordinates or "
        "equilibrium conformations.\n" %
        (metadata.get("shape", "unknown"),
         metadata.get("height_nm", "unknown"),
         metadata.get("maximum_diameter_nm", "unknown"),
         metadata.get("minimum_diameter_nm", "unknown"),
         metadata.get("layers", "unknown"), len(ordered), len(strands),
         structure_format, len(atom_atoms),
         geometry.get("source", "DNAxiS"),
         len(actual_ring_counts),
         ("rounded polygon with straight edges and local bend windows"
          if frame_export else "closed curved ring"),
         os.path.basename(paths["cylindrical_model"]),
         os.path.basename(paths["routing_model_multi"])))
    _atomic_write_many(files)
    return {
        "paths": paths,
        "structure_path": structure_path,
        "structure_format": structure_format,
        "all_atom_count": len(atom_atoms),
        "nucleotides": len(ordered),
        "strands": len(strands),
        "assigned_bases": assigned,
        "removed_geometry_bases": removed_count,
        "model_paths": [paths["cylindrical_model"],
                        paths["routing_model_multi"]]}

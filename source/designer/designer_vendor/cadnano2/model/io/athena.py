"""ATHENA wireframe-design integration.

The original ATHENA Qt5 GUI embeds PERDIX, DAEDALUS2, METIS and TALOS as
command-line backends.  cadnano calls those backends out of process so a
backend failure cannot take down the main Qt application.  The target-space
nucleotide frames are embedded (gzip + base64) in the resulting JSON, which
makes later ATHENA 3D exports independent of the original project folder.
"""

import base64
import gzip
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ATHENA_METADATA_VERSION = 1
OXDNA_LENGTH_ANGSTROM = 8.518
OXDNA_AXIAL_STEP = 3.4 / OXDNA_LENGTH_ANGSTROM


def athena_root():
    return os.path.abspath(os.path.join(
        os.path.dirname(__file__), os.pardir, os.pardir,
        "third_party", "athena"))


def preset_shapes():
    """Return all bundled ATHENA PLY examples in display order."""
    result = []
    root = os.path.join(athena_root(), "sample_inputs")
    for dimension in ("2D", "3D"):
        folder = os.path.join(root, dimension)
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            if not name.lower().endswith(".ply"):
                continue
            stem = os.path.splitext(name)[0]
            label = stem
            if len(stem) > 3 and stem[:2].isdigit() and stem[2] == "_":
                label = stem[3:]
            result.append({
                "dimension": dimension,
                "label": label.replace("_", " ").title(),
                "path": os.path.join(folder, name)})
    return result


def inspect_ply_dimension(path):
    """Classify an ASCII PLY as 2D (all z=0) or 3D."""
    vertex_count = None
    vertex_properties = []
    in_vertex = False
    with open(path, "r", encoding="utf-8", errors="strict") as source:
        first = source.readline().strip()
        if first != "ply":
            raise ValueError("The selected file is not a PLY mesh.")
        format_line = source.readline().strip().lower()
        if "ascii" not in format_line:
            raise ValueError(
                "ATHENA currently accepts an ASCII PLY mesh for dimension "
                "detection.")
        while True:
            line = source.readline()
            if not line:
                raise ValueError("The PLY header is incomplete.")
            fields = line.strip().split()
            if fields[:1] == ["element"]:
                in_vertex = len(fields) >= 3 and fields[1] == "vertex"
                if in_vertex:
                    vertex_count = int(fields[2])
                    vertex_properties = []
            elif fields[:1] == ["property"] and in_vertex:
                vertex_properties.append(fields[-1])
            elif fields[:1] == ["end_header"]:
                break
        if vertex_count is None:
            raise ValueError("The PLY mesh has no vertex element.")
        try:
            z_index = vertex_properties.index("z")
        except ValueError:
            z_index = None
        is_3d = False
        for unused in range(vertex_count):
            fields = source.readline().strip().split()
            if len(fields) < len(vertex_properties):
                raise ValueError("The PLY vertex table is incomplete.")
            if z_index is not None and abs(float(fields[z_index])) > 1e-9:
                is_3d = True
    return "3D" if is_3d else "2D"


def read_ply_preview_mesh(path):
    """Return vertices and polygon faces for a lightweight UI preview."""
    with open(path, "r", encoding="utf-8", errors="strict") as source:
        lines = source.read().splitlines()
    vertex_count = None
    face_count = 0
    data_start = None
    for index, line in enumerate(lines):
        fields = line.split()
        if fields[:2] == ["element", "vertex"]:
            vertex_count = int(fields[2])
        elif fields[:2] == ["element", "face"]:
            face_count = int(fields[2])
        elif fields[:1] == ["end_header"]:
            data_start = index + 1
            break
    if vertex_count is None or data_start is None:
        raise ValueError("PLY vertex table is missing.")
    vertices = []
    for line in lines[data_start:data_start + vertex_count]:
        values = line.split()
        if len(values) < 3:
            raise ValueError("PLY vertex table is incomplete.")
        vertices.append(tuple(float(value) for value in values[:3]))
    faces = []
    face_start = data_start + vertex_count
    for line in lines[face_start:face_start + face_count]:
        values = [int(value) for value in line.split()]
        if values and len(values) >= values[0] + 1:
            face = values[1:values[0] + 1]
            if len(face) >= 2 and all(
                    0 <= value < len(vertices) for value in face):
                faces.append(face)
    if not vertices or not faces:
        raise ValueError("PLY contains no drawable polygon faces.")
    return vertices, faces


def estimate_scaffold_minimum(path, edge_type, shortest_edge_bp):
    """Estimate the paired scaffold nt required by the target wireframe.

    DX edges contain two scaffolded duplexes and 6HB edges contain six.  The
    backend can add a small number of unpaired vertex nucleotides, so this is
    deliberately reported as a hard geometric minimum rather than an exact
    final scaffold length.
    """
    vertices, faces = read_ply_preview_mesh(path)
    edges = set()
    for face in faces:
        for first, second in zip(face, face[1:] + face[:1]):
            edges.add(tuple(sorted((first, second))))
    lengths = []
    for first, second in edges:
        delta = [vertices[first][axis] - vertices[second][axis]
                 for axis in range(3)]
        length = math.sqrt(sum(value * value for value in delta))
        if length > 1e-12:
            lengths.append(length)
    if not lengths:
        raise ValueError("PLY contains no non-zero wireframe edges.")
    shortest = min(lengths)
    scaled_lengths = [max(1, int(round(
        float(shortest_edge_bp) * length / shortest)))
                      for length in lengths]
    edge_type = str(edge_type).upper()
    if edge_type == "DX":
        duplexes = 2
    elif edge_type == "6HB":
        duplexes = 6
    else:
        raise ValueError("Unsupported wireframe edge type: %s" % edge_type)
    return duplexes * sum(scaled_lengths)


def recommended_engine(dimension, edge_type):
    table = {
        ("2D", "DX"): "PERDIX",
        ("3D", "DX"): "DAEDALUS2",
        ("2D", "6HB"): "METIS",
        ("3D", "6HB"): "TALOS"}
    try:
        return table[(str(dimension).upper(), str(edge_type).upper())]
    except KeyError:
        raise ValueError("Unsupported ATHENA dimension/edge combination.")


def safe_name(value):
    name = "".join(char if char.isalnum() or char in "-_" else "_"
                   for char in str(value)).strip("_")
    return name or "wireframe-design"


def wireframe_output_name(shape_name, edge_type):
    """Return a filename stem that always identifies shape and edge type."""
    shape = safe_name(shape_name)
    edge = safe_name(str(edge_type).upper())
    suffix = "_" + edge
    if shape.upper().endswith(suffix.upper()):
        return shape
    return shape + suffix


def backend_command(engine, output_dir, ply_path, scaffold_path,
                    edge_length):
    """Build the official nine-argument ATHENA backend command."""
    engine = str(engine).upper()
    settings = {
        "PERDIX": (1, 1),
        "DAEDALUS2": (1, 2),
        "METIS": (3, 2),
        "TALOS": (2, 1)}
    if engine not in settings:
        raise ValueError("Unknown ATHENA engine: %s" % engine)
    sections, vertex_design = settings[engine]
    suffix = ".exe" if platform.system() == "Windows" else ""
    executable = os.path.join(
        athena_root(), "tools", engine, engine + suffix)
    if not os.path.isfile(executable):
        raise FileNotFoundError(
            "The bundled %s backend is missing: %s" %
            (engine, executable))
    command = [
        executable, str(output_dir), str(ply_path), str(scaffold_path),
        str(sections), str(vertex_design), "0", str(int(edge_length)),
        "0.0", "s"]
    if platform.system() == "Darwin" and platform.machine() == "arm64":
        # ATHENA's published macOS backends are x86_64.  Rosetta runs them
        # reliably without loading their Qt5 GUI into cadnano's Qt process.
        command = ["/usr/bin/arch", "-x86_64"] + command
    return command


def _run_backend(command, progress=None, cancelled=None):
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1)
    lines = []
    while True:
        line = process.stdout.readline()
        if line:
            lines.append(line)
            if progress is not None:
                progress(line.rstrip())
        if cancelled is not None and cancelled():
            process.terminate()
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                process.kill()
            raise RuntimeError("ATHENA design was cancelled.")
        if line == "" and process.poll() is not None:
            break
    return_code = process.wait()
    output = "".join(lines)
    if return_code != 0:
        tail = "\n".join(output.splitlines()[-20:])
        raise RuntimeError(
            "The %s backend failed with exit code %d.\n\n%s" %
            (os.path.basename(command[-10] if command[:2] ==
                              ["/usr/bin/arch", "-x86_64"] else command[0]),
             return_code, tail))
    return output


def _json_nodes(obj):
    nodes = {}
    for virtual_helix in obj.get("vstrands", []):
        vh_num = int(virtual_helix["num"])
        for strand_type in ("scaf", "stap"):
            for idx, record in enumerate(virtual_helix.get(strand_type, [])):
                if record == [-1, -1, -1, -1]:
                    continue
                key = (strand_type, vh_num, idx)
                nodes[key] = {
                    "prev": (None if int(record[0]) < 0 else
                             (strand_type, int(record[0]), int(record[1]))),
                    "next": (None if int(record[2]) < 0 else
                             (strand_type, int(record[2]), int(record[3]))),
                    "across": None}
    for key, record in nodes.items():
        other_type = "stap" if key[0] == "scaf" else "scaf"
        other = (other_type, key[1], key[2])
        if other in nodes:
            record["across"] = other
    return nodes


def _read_cndo(path):
    sections = {}
    current = None
    with open(path, "r", encoding="utf-8") as source:
        for raw_line in source:
            line = raw_line.strip()
            if not line:
                current = None
                continue
            if line.startswith("dnaTop,"):
                current = "dnaTop"
                sections[current] = []
                continue
            if line.startswith("dNode,"):
                current = "dNode"
                sections[current] = []
                continue
            if line.startswith("triad,"):
                current = "triad"
                sections[current] = []
                continue
            if line.startswith("id_nt,"):
                current = "id_nt"
                sections[current] = []
                continue
            if current is not None:
                sections[current].append(line.split(","))
    topology = {}
    for row in sections.get("dnaTop", []):
        nt_id = int(row[1])
        topology[nt_id] = {
            "prev": int(row[2]) if int(row[2]) > 0 else None,
            "next": int(row[3]) if int(row[3]) > 0 else None,
            "across": int(row[4]) if int(row[4]) > 0 else None,
            "base": row[5].upper()}
    dnodes = dict((int(row[0]), tuple(float(value) for value in row[1:4]))
                  for row in sections.get("dNode", []))
    triads = dict((int(row[0]), tuple(float(value) for value in row[1:10]))
                  for row in sections.get("triad", []))
    pairs = [(int(row[0]), int(row[1]), int(row[2]))
             for row in sections.get("id_nt", [])]
    return topology, dnodes, triads, pairs


def _component_orders(nodes):
    unseen = set(nodes)
    orders = []
    while unseen:
        seed = next(iter(unseen))
        component = {seed}
        pending = [seed]
        while pending:
            node = pending.pop()
            for neighbor in (nodes[node]["prev"], nodes[node]["next"]):
                if neighbor is not None and neighbor not in component:
                    component.add(neighbor)
                    pending.append(neighbor)
        starts = [node for node in component
                  if nodes[node]["prev"] is None]
        start = starts[0] if starts else min(component)
        order = []
        node = start
        while node is not None and node not in order:
            order.append(node)
            node = nodes[node]["next"]
        if len(order) != len(component):
            raise ValueError("ATHENA produced an invalid strand topology.")
        unseen.difference_update(component)
        orders.append(order)
    return orders


def _topology_mapping(json_obj, cndo_topology):
    """Map every CanDo nucleotide id to its caDNAno (type, vh, idx)."""
    json_nodes = _json_nodes(json_obj)
    json_orders = _component_orders(json_nodes)
    cndo_orders = _component_orders(cndo_topology)
    json_scaffolds = [order for order in json_orders
                      if order and order[0][0] == "scaf"]
    if not json_scaffolds:
        raise ValueError("ATHENA JSON contains no scaffold strand.")
    json_scaffold = max(json_scaffolds, key=len)
    cndo_scaffold = max(cndo_orders, key=len)
    if len(json_scaffold) != len(cndo_scaffold):
        raise ValueError("ATHENA JSON/CanDo scaffold lengths disagree.")

    json_position = {
        node: (len(order), position)
        for order in json_orders for position, node in enumerate(order)}
    cndo_position = {
        node: (len(order), position)
        for order in cndo_orders for position, node in enumerate(order)}

    def signatures(nodes, positions, order):
        return [positions[nodes[node]["across"]]
                if nodes[node]["across"] is not None else None
                for node in order]

    json_signature = signatures(
        json_nodes, json_position, json_scaffold)
    cndo_signature = signatures(
        cndo_topology, cndo_position, cndo_scaffold)
    if json_signature != cndo_signature:
        reversed_scaffold = list(reversed(cndo_scaffold))
        reversed_signature = signatures(
            cndo_topology, cndo_position, reversed_scaffold)
        if json_signature != reversed_signature:
            raise ValueError(
                "ATHENA JSON and CanDo topology could not be aligned.")
        cndo_scaffold = reversed_scaffold

    cndo_to_json = dict(zip(cndo_scaffold, json_scaffold))
    for cndo_node, json_node in zip(cndo_scaffold, json_scaffold):
        cndo_across = cndo_topology[cndo_node]["across"]
        json_across = json_nodes[json_node]["across"]
        if cndo_across is not None and json_across is not None:
            cndo_to_json[cndo_across] = json_across

    cndo_order_for_node = {
        node: order for order in cndo_orders for node in order}
    for json_order in json_orders:
        if not json_order or json_order[0][0] != "stap":
            continue
        anchors = []
        for json_pos, json_node in enumerate(json_order):
            for cndo_node, mapped_json in cndo_to_json.items():
                if mapped_json == json_node:
                    cndo_order = cndo_order_for_node[cndo_node]
                    anchors.append(
                        (json_pos, cndo_order.index(cndo_node), cndo_order))
                    break
        if not anchors:
            raise ValueError("An ATHENA staple has no paired anchor.")
        cndo_order = anchors[0][2]
        if len(cndo_order) != len(json_order):
            raise ValueError("ATHENA staple lengths disagree.")
        direct = all(json_pos == cndo_pos
                     for json_pos, cndo_pos, unused in anchors)
        reverse = all(json_pos == len(json_order) - 1 - cndo_pos
                      for json_pos, cndo_pos, unused in anchors)
        if not direct and not reverse:
            raise ValueError("An ATHENA staple could not be aligned.")
        if reverse:
            cndo_order = list(reversed(cndo_order))
        for cndo_node, json_node in zip(cndo_order, json_order):
            existing = cndo_to_json.get(cndo_node)
            if existing is not None and existing != json_node:
                raise ValueError("Conflicting ATHENA nucleotide mapping.")
            cndo_to_json[cndo_node] = json_node
    if len(cndo_to_json) != len(cndo_topology):
        raise ValueError("ATHENA nucleotide mapping is incomplete.")
    return cndo_to_json, json_nodes


def _normalize(vector):
    length = math.sqrt(sum(value * value for value in vector))
    if length <= 1e-12:
        return (1.0, 0.0, 0.0)
    return tuple(value / length for value in vector)


def _lerp(first, second, fraction):
    return tuple(first[index] +
                 (second[index] - first[index]) * fraction
                 for index in range(3))


def _cndo_frames(topology, dnodes, triads, pairs):
    frames = {}
    for bp_id, scaffold_id, staple_id in pairs:
        if bp_id not in dnodes or bp_id not in triads:
            continue
        center = tuple(value / OXDNA_LENGTH_ANGSTROM
                       for value in dnodes[bp_id])
        triad = triads[bp_id]
        radial = _normalize(triad[3:6])
        tangent = _normalize(triad[6:9])
        for nt_id, sign in ((scaffold_id, 1.0), (staple_id, -1.0)):
            frames[nt_id] = {
                "pos": [center[index] + sign * 0.60 * radial[index]
                        for index in range(3)],
                "a1": [-sign * radial[index] for index in range(3)],
                "a3": [sign * tangent[index] for index in range(3)]}

    # ATHENA represents vertex poly-T segments as unpaired nucleotides.  CanDo
    # has no dNode for them, so interpolate between their paired neighbours.
    for order in _component_orders(topology):
        missing = [position for position, node in enumerate(order)
                   if node not in frames]
        while missing:
            start = missing[0]
            end = start
            while end + 1 < len(order) and order[end + 1] not in frames:
                end += 1
            before = start - 1 if start > 0 else None
            after = end + 1 if end + 1 < len(order) else None
            count = end - start + 1
            for offset, position in enumerate(range(start, end + 1), 1):
                if before is not None and after is not None:
                    fraction = float(offset) / float(count + 1)
                    first, second = frames[order[before]], frames[order[after]]
                    pos = _lerp(first["pos"], second["pos"], fraction)
                    a1 = _normalize(_lerp(first["a1"], second["a1"],
                                         fraction))
                    a3 = _normalize(_lerp(first["a3"], second["a3"],
                                         fraction))
                elif before is not None:
                    first = frames[order[before]]
                    pos = [first["pos"][axis] +
                           offset * OXDNA_AXIAL_STEP * first["a3"][axis]
                           for axis in range(3)]
                    a1, a3 = first["a1"], first["a3"]
                elif after is not None:
                    second = frames[order[after]]
                    distance = count - offset + 1
                    pos = [second["pos"][axis] -
                           distance * OXDNA_AXIAL_STEP * second["a3"][axis]
                           for axis in range(3)]
                    a1, a3 = second["a1"], second["a3"]
                else:
                    raise ValueError(
                        "ATHENA produced an entirely unpaired strand.")
                frames[order[position]] = {
                    "pos": list(pos), "a1": list(a1), "a3": list(a3)}
            missing = [position for position, node in enumerate(order)
                       if node not in frames]
    return frames


def _frame_key(json_node):
    strand_type = "scaffold" if json_node[0] == "scaf" else "staple"
    return "%s:%d:%d" % (strand_type, json_node[1], json_node[2])


def build_geometry_payload(json_obj, cndo_path):
    topology, dnodes, triads, pairs = _read_cndo(cndo_path)
    cndo_to_json, json_nodes = _topology_mapping(json_obj, topology)
    cndo_frames = _cndo_frames(topology, dnodes, triads, pairs)
    frames = {}
    for cndo_id, json_node in cndo_to_json.items():
        frames[_frame_key(json_node)] = cndo_frames[cndo_id]
    base_keys = sorted(frames)
    fingerprint = hashlib.sha256(
        "\n".join(base_keys).encode("utf-8")).hexdigest()
    return {
        "coordinate_units": "oxDNA",
        "source": "ATHENA CanDo pseudo-atomic target geometry",
        "unrelaxed": True,
        "base_count": len(frames),
        "base_set_fingerprint": fingerprint,
        "frames": frames}


def encode_geometry_payload(payload):
    raw = json.dumps(payload, separators=(",", ":"),
                     sort_keys=True).encode("utf-8")
    return base64.b64encode(gzip.compress(raw, 9)).decode("ascii")


def decode_geometry_payload(metadata):
    encoded = metadata.get("geometry_data")
    if not isinstance(encoded, str) or not encoded:
        raise ValueError("This ATHENA design has no embedded 3D mapping.")
    try:
        raw = gzip.decompress(base64.b64decode(encoded.encode("ascii")))
        return json.loads(raw.decode("utf-8"))
    except Exception as error:
        raise ValueError("The embedded ATHENA 3D mapping is damaged.") from error


def _find_numpy_python():
    candidates = [
        sys.executable,
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        shutil.which("python3")]
    seen = set()
    for candidate in candidates:
        if not candidate or candidate in seen or not os.path.isfile(candidate):
            continue
        seen.add(candidate)
        probe = subprocess.run(
            [candidate, "-c", "import numpy"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if probe.returncode == 0:
            return candidate
    return None


def _generate_original_pdb(cndo_path, output_dir):
    python = _find_numpy_python()
    if python is None:
        raise RuntimeError(
            "ATHENA PDB generation requires a Python interpreter with NumPy.")
    runner = os.path.join(
        athena_root(), "pdbgen", "pdbgen_runner.py")
    result = subprocess.run(
        [python, runner, cndo_path, output_dir],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "ATHENA PDB generation failed.\n\n%s" %
            "\n".join(result.stdout.splitlines()[-20:]))
    stem = os.path.splitext(os.path.basename(cndo_path))[0]
    path = os.path.join(output_dir, stem + ".pdb")
    if os.path.isfile(path):
        return path
    # PDBGen deliberately skips the conventional single-model PDB when a
    # structure has more than 63 strands (the one-character chain-ID limit).
    # Its segid output preserves every strand and is safe to convert to mmCIF.
    for suffix in ("-segid.pdb", "-multimodel.pdb"):
        fallback = os.path.join(output_dir, stem + suffix)
        if os.path.isfile(fallback):
            return fallback
    raise RuntimeError("ATHENA PDB generator produced no usable PDB file.")


def _pdb_to_mmcif(pdb_path, cif_path, name):
    rows = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as source:
        for line in source:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            rows.append({
                "group": line[0:6].strip(),
                "atom": line[12:16].strip() or "?",
                "residue": line[17:20].strip() or "DN",
                "chain": ("S" + line[72:76].strip()
                          if line[72:76].strip() else
                          (line[21:22].strip() or "A")),
                "residue_key": line[22:27],
                "x": line[30:38].strip(),
                "y": line[38:46].strip(),
                "z": line[46:54].strip(),
                "element": line[76:78].strip() or
                           (line[12:16].strip()[:1] or "C")})
    lines = [
        "data_%s" % safe_name(name),
        "_struct.title 'ATHENA pseudo-atomic DNA model'",
        "#", "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_PDB_model_num"]
    residue_numbers = {}
    previous_residue = {}
    entity_ids = {}
    for atom_id, row in enumerate(rows, 1):
        chain = row["chain"]
        if chain not in entity_ids:
            entity_ids[chain] = len(entity_ids) + 1
        residue_key = row["residue_key"]
        if previous_residue.get(chain) != residue_key:
            residue_numbers[chain] = residue_numbers.get(chain, 0) + 1
            previous_residue[chain] = residue_key
        lines.append(
            "%s %d %s %s %s %s %d %d %s %s %s 1.00 0.00 1" %
            (row["group"], atom_id, row["element"], row["atom"],
             row["residue"], chain, entity_ids[chain],
             residue_numbers[chain],
             row["x"], row["y"], row["z"]))
    lines.extend(("#", ""))
    with open(cif_path, "w", encoding="utf-8", newline="\n") as output:
        output.write("\n".join(lines))
    return len(rows)


def _copy_model_outputs(source_dir, models_dir, output_name):
    os.makedirs(models_dir, exist_ok=True)
    copied = []
    rules = (
        ("_cylinder_model.bild", output_name + "_cylindrical_model.bild"),
        ("_routing_multi.bild", output_name + "_routing_model_multi.bild"),
        ("_routing_two.bild", output_name + "_routing_model_two.bild"),
        ("_atomic_model_multi.bild",
         output_name + "_pseudo_atomic_model_multi.bild"),
        ("_atomic_model_two.bild",
         output_name + "_pseudo_atomic_model_two.bild"))
    for source_name in os.listdir(source_dir):
        for suffix, target_name in rules:
            if source_name.endswith(suffix):
                target = os.path.join(models_dir, target_name)
                shutil.copy2(os.path.join(source_dir, source_name), target)
                copied.append(target)
                break
    return copied


def create_project(spec, progress=None, cancelled=None):
    """Run ATHENA and create the clean project folder agreed by the user."""
    ply_path = os.path.abspath(spec["ply_path"])
    dimension = inspect_ply_dimension(ply_path)
    edge_type = str(spec["edge_type"]).upper()
    engine = str(spec.get("engine") or
                 recommended_engine(dimension, edge_type)).upper()
    if engine != recommended_engine(dimension, edge_type):
        # Manual override is allowed only when it is dimension/edge compatible.
        raise ValueError(
            "%s is not compatible with %s + %s." %
            (engine, dimension, edge_type))
    name = safe_name(spec.get("name") or
                     os.path.splitext(os.path.basename(ply_path))[0])
    output_name = wireframe_output_name(name, edge_type)
    project_root = os.path.abspath(spec["project_root"])
    input_dir = os.path.join(project_root, "input")
    models_dir = os.path.join(project_root, "models")
    internal_dir = os.path.join(project_root, ".athena")
    for folder in (project_root, input_dir, models_dir, internal_dir):
        os.makedirs(folder, exist_ok=True)

    project_ply = os.path.join(input_dir, output_name + ".ply")
    shutil.copy2(ply_path, project_ply)
    scaffold_sequence = "".join(
        char for char in str(spec["scaffold_sequence"]).upper()
        if char in "ACGT")
    if not scaffold_sequence:
        raise ValueError("The selected scaffold template is empty.")

    with tempfile.TemporaryDirectory(prefix="cadnano-athena-") as temp_dir:
        # ATHENA requires a sequence while routing, but the user explicitly
        # does not want a scaffold-sequence output file. Keep this input only
        # in the temporary backend directory and discard it on completion.
        scaffold_path = os.path.join(temp_dir, "scaffold_template.txt")
        with open(scaffold_path, "w", encoding="ascii",
                  newline="\n") as output:
            output.write(scaffold_sequence + "\n")
        command = backend_command(
            engine, temp_dir, project_ply, scaffold_path,
            spec["edge_length"])
        log = _run_backend(
            command, progress=progress, cancelled=cancelled)
        with open(os.path.join(
                internal_dir, output_name + "_backend.log"), "w",
                  encoding="utf-8", newline="\n") as output:
            output.write(log)
        json_files = sorted(Path(temp_dir).glob("*.json"))
        cndo_files = sorted(Path(temp_dir).glob("*.cndo"))
        if not json_files or not cndo_files:
            if "User-defined sequences are short" in log:
                required_values = re.findall(
                    r"# of nts in (?:the )?scaffold\s*:\s*(\d+)", log)
                required = (required_values[-1]
                            if required_values else "more")
                raise RuntimeError(
                    "The selected scaffold is too short. This routing "
                    "requires approximately %s nt. Choose a longer scaffold "
                    "or reduce the shortest edge length." % required)
            raise RuntimeError(
                "ATHENA completed without producing JSON/CanDo output.")
        with open(str(json_files[0]), "r", encoding="utf-8") as source:
            json_obj = json.load(source)
        source_cndo = os.path.join(
            internal_dir, output_name + "_source.cndo")
        shutil.copy2(str(cndo_files[0]), source_cndo)
        geometry = build_geometry_payload(json_obj, source_cndo)

        parameters = {
            "format": "cadnano-athena-project-v1",
            "name": name,
            "output_name": output_name,
            "dimension": dimension,
            "edge_type": edge_type,
            "engine": engine,
            "edge_length_bp": int(spec["edge_length"]),
            "scaffold_template": spec["scaffold_name"],
            "scaffold_template_length": len(scaffold_sequence),
            "unrelaxed": True}
        with open(os.path.join(
                input_dir, output_name + "_parameters.json"), "w",
                  encoding="utf-8", newline="\n") as output:
            json.dump(parameters, output, indent=2, sort_keys=True)
            output.write("\n")

        metadata = dict(parameters)
        metadata.update({
            "metadata_version": ATHENA_METADATA_VERSION,
            "project_root": project_root,
            "input_ply": os.path.relpath(project_ply, project_root),
            "source_cndo": os.path.relpath(source_cndo, project_root),
            "geometry_encoding": "gzip+base64+json",
            "geometry_data": encode_geometry_payload(geometry)})
        json_obj["athena_metadata"] = metadata
        json_path = os.path.join(project_root, output_name + ".json")
        with open(json_path, "w", encoding="utf-8", newline="\n") as output:
            json.dump(json_obj, output, separators=(",", ":"))

        _copy_model_outputs(temp_dir, models_dir, output_name)
        # Keep a readable copy of the geometry payload for diagnostics while
        # the embedded compressed copy makes the JSON self-contained.
        with open(os.path.join(
                internal_dir, output_name + "_geometry.json"), "w",
                  encoding="utf-8", newline="\n") as output:
            json.dump(geometry, output, separators=(",", ":"))

    return {
        "json_path": json_path,
        "project_root": project_root,
        "model_paths": sorted(
            str(path) for path in Path(models_dir).glob("*")),
        "metadata": metadata}

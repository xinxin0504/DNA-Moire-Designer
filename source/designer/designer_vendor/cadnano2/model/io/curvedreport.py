"""Data extraction for the Curved Design completion report."""

import json
import os
import subprocess
import sys

from .indelanalysis import (curved_pair_curvature_data,
                            frame_pair_curvature_data,
                            generated_single_helix_distribution_data)


def _domain_indels(design, rows, lattice, frame_only=False):
    """Return final signed indels for every complete native domain."""
    domain_size = 7 if lattice == "honeycomb" else 8
    ring_records = {
        int(record["helix"]): record
        for record in (design.get("curvature_indels", {}) or {}).get(
            "rings", [])
        if "helix" in record}
    frame_windows = []
    if frame_only:
        plan = dict((design.get("curved_metadata", {}) or {}).get(
            "frame_plan", {}) or {})
        frame_windows = [
            (float(centre)-float(length)/2.0,
             float(centre)+float(length)/2.0, vertex)
            for vertex, (centre, length) in enumerate(zip(
                plan.get("vertex_native_centres", ()),
                plan.get("bend_length_bp", ())))]
    result = []
    for helix, row in sorted(rows.items()):
        loops = list(row.get("loop", []))
        skips = list(row.get("skip", []))
        record = ring_records.get(helix, {})
        nominal_size = int(record.get(
            "nominal_bases", min(len(loops), len(skips))))
        nominal_size = max(0, min(nominal_size, len(loops), len(skips)))
        complete_domain_count = nominal_size // domain_size
        values = []
        domain_indices = []
        domain_vertices = []
        for domain in range(complete_domain_count):
            first = domain * domain_size
            last = first + domain_size
            centre = 0.5*(first+last)
            vertex = next((item[2] for item in frame_windows
                           if item[0] <= centre <= item[1]), None)
            if frame_only and vertex is None:
                continue
            insertions = sum(max(0, int(loops[index]))
                             for index in range(first, last))
            deletions = sum(max(0, -int(skips[index]))
                            for index in range(first, last))
            values.append(int(insertions - deletions))
            domain_indices.append(domain)
            domain_vertices.append(vertex)
        result.append({
            "helix": int(helix),
            "nominal_bases": nominal_size,
            "domain_size_bp": domain_size,
            "values": values, "domain_indices": domain_indices,
            "domain_vertices": domain_vertices,
            "scope": "bend-windows" if frame_only else "whole-loop"})
    return result


def _staple_centers(rows):
    nodes = set()
    five_prime = []
    for helix, row in rows.items():
        for index, connection in enumerate(row.get("stap", [])):
            if connection != [-1, -1, -1, -1]:
                nodes.add((helix, index))
                if int(connection[0]) < 0:
                    five_prime.append((helix, index))

    def follow(start):
        chain = []
        seen = set()
        node = start
        while node in nodes and node not in seen:
            seen.add(node)
            chain.append(node)
            connection = rows[node[0]]["stap"][node[1]]
            next_helix = int(connection[2])
            next_index = int(connection[3])
            if next_helix < 0:
                break
            node = (next_helix, next_index)
        return chain

    chains = []
    assigned = set()
    for start in sorted(five_prime):
        if start in assigned:
            continue
        chain = follow(start)
        if chain:
            chains.append(chain)
            assigned.update(chain)
    for start in sorted(nodes - assigned):
        chain = follow(start)
        if chain:
            chains.append(chain)
            assigned.update(chain)

    centers = []
    for chain in chains:
        weighted_nodes = []
        actual_length = 0
        for helix, index in chain:
            row = rows[helix]
            loops = row.get("loop", [])
            skips = row.get("skip", [])
            loop = int(loops[index]) if index < len(loops) else 0
            skip = int(skips[index]) if index < len(skips) else 0
            weight = max(0, 1 + max(0, loop) - max(0, -skip))
            weighted_nodes.append(((helix, index), weight))
            actual_length += weight
        if not weighted_nodes:
            continue
        midpoint = max(0.0, (actual_length - 1) / 2.0)
        cumulative = 0
        center = weighted_nodes[len(weighted_nodes) // 2][0]
        for node, weight in weighted_nodes:
            if weight > 0 and cumulative + weight > midpoint:
                center = node
                break
            cumulative += weight
        centers.append({
            "helix": int(center[0]),
            "base": int(center[1]),
            "length": int(actual_length)})
    return sorted(centers, key=lambda item: (
        item["helix"], item["base"], item["length"]))


def curved_report_data(json_path):
    """Return graph-ready final domain-indel and staple-center data."""
    with open(json_path, "r", encoding="utf-8") as source:
        design = json.load(source)
    rows = {int(row["num"]): row
            for row in design.get("vstrands", [])}
    lattice = str(
        design.get("lattice") or
        (design.get("curved_metadata", {}) or {}).get(
            "lattice", "honeycomb")).lower()
    metadata = dict(design.get("curved_metadata", {}) or {})
    is_frame = metadata.get("format") == "cadnano-frame-project-v1"
    pair_rows, pair_summary = (
        frame_pair_curvature_data(design) if is_frame else
        curved_pair_curvature_data(design))
    return {
        "domain_indels": _domain_indels(
            design, rows, lattice, frame_only=is_frame),
        "pair_curvature_rows": pair_rows,
        "pair_curvature_summary": pair_summary,
        "single_helix_distribution":
            generated_single_helix_distribution_data(
                design, frame_only=is_frame),
        "staple_nick_optimization": dict(
            metadata.get("frame_staple_nick_optimization", {}) or
            metadata.get("curved_staple_nick_optimization", {}) or {}),
        "frame_straight_common_mode_remove_twist": dict(
            metadata.get(
                "frame_straight_common_mode_remove_twist", {}) or {}),
        "report_mode": "frame" if is_frame else "curved",
        "analysis_scope": ("bend-windows-only" if is_frame else
                           "whole-closed-loop"),
        "staples": _staple_centers(rows),
        "helices": sorted(rows),
        "lattice": lattice,
        "base_count": max(
            [len(row.get("stap", [])) for row in rows.values()] or [0])}


def start_curved_report_export(json_path, summary_lines, output_path):
    """Start a detached report renderer that survives cadnano shutdown."""
    # A project can be regenerated into the same directory.  Do not let the
    # GUI watcher mistake the previous PNG for the report belonging to the
    # newly started worker.  The worker itself writes to ``.writing`` and
    # atomically renames only after painting has completed.
    for stale_path in (output_path, output_path + ".writing"):
        if os.path.exists(stale_path):
            try:
                os.unlink(stale_path)
            except OSError:
                pass
    environment = os.environ.copy()
    environment.setdefault("QT_QPA_PLATFORM", "offscreen")
    worker_path = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "..", "views",
        "curvedreportworker.py"))
    command = [
        sys.executable, worker_path,
        os.path.abspath(json_path), os.path.abspath(output_path),
        json.dumps(list(summary_lines), ensure_ascii=False)]
    kwargs = {
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
        "close_fds": True,
        "env": environment}
    if os.name == "nt":
        kwargs["creationflags"] = (
            getattr(subprocess, "DETACHED_PROCESS", 0x00000008) |
            getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x00000200))
    else:
        kwargs["start_new_session"] = True
    return subprocess.Popen(command, **kwargs)

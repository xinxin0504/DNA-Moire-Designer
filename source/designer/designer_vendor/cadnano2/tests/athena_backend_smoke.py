"""Smoke test for the bundled ATHENA backend and self-contained metadata."""

import json
import os
import tempfile

from cadnano2.model.io.athena import (create_project,
                                      decode_geometry_payload,
                                      preset_shapes)


assert len(preset_shapes()) >= 90
with tempfile.TemporaryDirectory() as directory:
    shape = next(item for item in preset_shapes()
                 if item["dimension"] == "3D" and
                 item["label"].lower() == "tetrahedron")
    summary = create_project({
        "ply_path": shape["path"],
        "edge_type": "DX",
        "engine": "DAEDALUS2",
        "edge_length": 42,
        "scaffold_name": "test scaffold",
        "scaffold_sequence": "ACGT" * 250,
        "name": "tetrahedron",
        "project_root": directory})
    assert os.path.isfile(summary["json_path"])
    assert os.path.isfile(summary["structure_path"])
    assert len(summary["model_paths"]) == 5
    encoded = json.load(open(summary["json_path"], encoding="utf-8"))
    metadata = encoded["athena_metadata"]
    geometry = decode_geometry_payload(metadata)
    assert geometry["base_count"] > 100
    spans = []
    positions = [frame["pos"] for frame in geometry["frames"].values()]
    for axis in range(3):
        spans.append(max(pos[axis] for pos in positions) -
                     min(pos[axis] for pos in positions))
    assert all(span > 1.0 for span in spans), spans
    all_files = [name for unused_root, unused_dirs, names
                 in os.walk(directory) for name in names]
    assert not any(name.lower().endswith(".csv") for name in all_files)
    assert not any("scaffold" in name.lower() and
                   name.lower().endswith(".txt") for name in all_files)

print("ATHENA backend smoke: OK")

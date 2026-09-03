#!/usr/bin/env python3
"""Verify exact 5'-anchored scaffold-sequence JSON round-tripping."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMPANION = ROOT / "source" / "cadnano_companion"
sys.path.insert(0, str(COMPANION))
# The upstream headless compatibility loader imports ``dummyqt`` as a
# top-level package. The release companion retains upstream's exact dummy
# files; source-only QA uses a test-scoped compatibility shim instead.
sys.path.insert(0, str(COMPANION / "cadnano2"))
sys.path.insert(0, str(ROOT / "tests" / "headless_qt"))

from cadnano2 import cadnano, util
# Force cadnano's model-only compatibility layer for this source-level test.
# The packaged companion itself is separately exercised with a real PyQt GUI
# in --self-test mode.
util.qtFrameworkList = ['Dummy']
from cadnano2.model.document import Document
from cadnano2.model.io.legacydecoder import import_legacy_dict
from cadnano2.model.io.legacyencoder import legacy_dict_from_doc


def canonical(records):
    return sorted(
        ({"start_vh": int(item["start_vh"]),
          "start_idx": int(item["start_idx"]),
          "sequence": str(item["sequence"])} for item in records),
        key=lambda item: (item["start_vh"], item["start_idx"]))


def digest(records):
    data = json.dumps(canonical(records), sort_keys=True,
                      separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def synthetic_fixture():
    path = (ROOT / "source" / "designer" / "moire_design_core" /
            "resources" / "Square_Seed_2L_newtemplate.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    first_document = Document()
    part = import_legacy_dict(first_document, payload)
    records = []
    alphabet = "ACGT"
    scaffolds = sorted(
        (oligo for oligo in part.oligos()
         if not oligo.isStaple() and oligo.strand5p() is not None),
        key=lambda oligo: (oligo.strand5p().virtualHelix().number(),
                           oligo.strand5p().idx5Prime()))
    for number, oligo in enumerate(scaffolds[:3]):
        strand = oligo.strand5p()
        sequence = "".join(alphabet[(index + number) % 4]
                           for index in range(oligo.length()))
        records.append({
            "start_vh": strand.virtualHelix().number(),
            "start_idx": strand.idx5Prime(),
            "sequence": sequence,
        })
    if not records:
        raise RuntimeError("The bundled test template has no scaffold oligo.")
    payload["scaffold_sequences"] = records
    payload["moire_structure_metadata"] = {
        "roundtrip_test": "preserve unknown Designer metadata"}
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", nargs="?", type=Path)
    args = parser.parse_args()
    application = cadnano.initAppWithoutGui([])
    # cadnano2's historical headless object omits this GUI-owned signal even
    # though Document emits it unconditionally.
    class Signal:
        @staticmethod
        def emit(*unused):
            return None
    application.documentWasCreatedSignal = Signal()
    application.prefs.squareSteps = 200
    application.prefs.honeycombSteps = 200
    application.prefs.honeycombRows = 50
    application.prefs.honeycombCols = 50
    payload = (json.loads(args.fixture.read_text(encoding="utf-8"))
               if args.fixture else synthetic_fixture())
    expected = canonical(payload.get("scaffold_sequences", []))
    if not expected:
        raise SystemExit("Fixture contains no scaffold_sequences records.")
    document = Document()
    part = import_legacy_dict(document, payload)
    coords = [(int(row["row"]), int(row["col"]))
              for row in payload["vstrands"]]
    encoded = legacy_dict_from_doc(document, "roundtrip.json", coords)
    actual = canonical(encoded.get("scaffold_sequences", []))
    report = {
        "input_records": len(expected),
        "output_records": len(actual),
        "input_nt": sum(len(item["sequence"]) for item in expected),
        "output_nt": sum(len(item["sequence"]) for item in actual),
        "input_sha256": digest(expected),
        "output_sha256": digest(actual),
        "exact_records": expected == actual,
        "metadata_preserved": (encoded.get("moire_structure_metadata") ==
                               payload.get("moire_structure_metadata")),
        "helices": part.numberOfVirtualHelices(),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["exact_records"] and report["metadata_preserved"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

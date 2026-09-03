#!/usr/bin/env python3
"""Render cadnano's Illustrator SVG offscreen for a supplied JSON state."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parents[2]
CADNANO_ROOT = PROJECT_ROOT / "work" / "cadnano2-modified"
if not (CADNANO_ROOT / "__init__.py").is_file():
    CADNANO_ROOT = HERE.parents[1] / "designer_vendor" / "cadnano2"
if not (CADNANO_ROOT / "__init__.py").is_file():
    raise RuntimeError("The packaged Designer cadnano engine was not found.")
sys.path.insert(0, str(CADNANO_ROOT.parent))

from cadnano2 import cadnano
from cadnano2.model.enum import LatticeType


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: vector_export_worker.py input.json output.svg")
    source = Path(sys.argv[1]).expanduser().resolve()
    target = Path(sys.argv[2]).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    application = cadnano.initAppWithGui(["moire-vector-export"])
    # Import the GUI-aware decoder only after QApplication exists; its legacy
    # styling imports may query the application font database.
    from cadnano2.model.io import decoder
    # DNA Moire Designer serializes every legacy view on the square caDNAno
    # canvas. Kagome and mixed designs use different strand routing, not a
    # Honeycomb caDNAno part. A 672-base canvas is divisible by both 21 and
    # 32, so the generic GUI decoder otherwise opens an invisible lattice
    # chooser and the offscreen export hangs indefinitely. Pin this dedicated
    # worker to the known square canvas.
    original_import = decoder.import_legacy_dict

    def import_square(document, payload, *unused_args, **unused_kwargs):
        return original_import(
            document, payload, LatticeType.Square,
            forceLatticeType=True)

    decoder.import_legacy_dict = import_square
    controller = next(iter(application.documentControllers))
    controller.openAfterMaybeSaveCallback(str(source))
    application.qApp.processEvents()
    controller._writeIllustratorSVG(str(target))
    application.qApp.processEvents()
    print(json.dumps({"svg": str(target)}, ensure_ascii=True))


if __name__ == "__main__":
    main()

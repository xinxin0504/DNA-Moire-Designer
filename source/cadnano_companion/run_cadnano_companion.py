#!/usr/bin/env python3
"""Launcher for the official-base cadnano companion application."""

from __future__ import annotations

import json
from importlib.metadata import PackageNotFoundError, version
import os
from pathlib import Path
import sys

from cadnano2 import cadnano


def self_test():
    # Importing encoder/decoder alone cannot catch missing frozen GUI modules.
    # Construct the same first document/window used by a normal launch while
    # forcing Qt's platform-independent backend for unattended builds.
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    application = cadnano.initAppWithGui([
        "cadnano2-self-test", "-platform", "offscreen"])
    try:
        release = version("cadnano2")
    except PackageNotFoundError:
        release = "2.4.13-source"
    from cadnano2.model.io import legacydecoder, legacyencoder
    report = {
        "application": "cadnano2 companion",
        "official_base_version": release,
        "package_root": str(Path(cadnano.__file__).resolve().parent),
        "legacy_decoder": str(Path(legacydecoder.__file__).resolve()),
        "legacy_encoder": str(Path(legacyencoder.__file__).resolve()),
        "designer_sequence_json_extension": True,
        "gui_initialized": bool(application.documentControllers),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    application.qApp.closeAllWindows()
    application.qApp.processEvents()
    application.qApp.quit()
    return 0


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    application = cadnano.initAppWithGui(argv)
    # The upstream launcher does not consume a positional JSON argument.
    # Supporting one is required for Designer's "Open cadnano" action and
    # does not change editing or file-format behavior.
    candidates = [Path(value).expanduser().resolve()
                  for value in argv[1:] if not value.startswith("-")]
    source = next((path for path in candidates if path.is_file()), None)
    if source is not None:
        controller = next(iter(application.documentControllers))
        controller.openAfterMaybeSaveCallback(str(source))
    application.exec_()
    return 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(self_test())
    raise SystemExit(main())

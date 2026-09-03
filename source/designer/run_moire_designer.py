#!/usr/bin/env python3
"""Development and frozen launcher for DNA Moiré Designer."""

import json
import os
from pathlib import Path
import re
import sys

from moire_runtime import (
    WORKER_MODULES, application_root, cadnano_executable,
    configure_designer_engine, dispatch_worker, source_root, tool_executable)


configure_designer_engine()


def self_test():
    """Run resource checks and construct the complete window off-screen."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import numpy
    import PyQt6
    import cadnano2
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import (
        QAbstractButton, QApplication, QComboBox, QGroupBox, QLabel,
        QTabWidget, QTableWidget, QTreeWidget, QWidget)
    from moire_design_core import models, structure, template
    from moire_designer import analysis_bulk
    from moire_designer.i18n import translate
    from moire_designer.mainwindow import MoireDesignerWindow

    app = QApplication.instance() or QApplication(
        ["DNA_Moire_Designer", "--self-test"])
    window = MoireDesignerWindow()
    window._localizer.retranslate()
    app.processEvents()
    window_initialized = (
        window.centralWidget() is not None
        and window.analysis_module_stack.count() == 1
    )

    # A source audit alone cannot detect a missing translations.json in a
    # frozen build.  Exercise the packaged catalog and scan the constructed
    # window so a Chinese presentation regression stops build.ps1.
    cjk = re.compile(r"[\u3400-\u9fff]")
    presentation_text = []

    def inspect(value, location):
        value = str(value or "")
        if cjk.search(value):
            presentation_text.append({"location": location, "text": value})

    objects = [window] + window.findChildren(QWidget) + \
        window.findChildren(QAction)
    for item in objects:
        name = item.objectName() or type(item).__name__
        if isinstance(item, QWidget) and item.isWindow():
            inspect(item.windowTitle(), name + ".windowTitle")
        if isinstance(item, (QLabel, QAbstractButton)):
            inspect(item.text(), name + ".text")
        if isinstance(item, QGroupBox):
            inspect(item.title(), name + ".title")
        if isinstance(item, QAction):
            inspect(item.text(), name + ".action")
        if isinstance(item, QComboBox):
            for index in range(item.count()):
                inspect(item.itemText(index), "%s.item[%d]" % (name, index))
        if isinstance(item, QTabWidget):
            for index in range(item.count()):
                inspect(item.tabText(index), "%s.tab[%d]" % (name, index))
        if isinstance(item, QTableWidget):
            for column in range(item.columnCount()):
                cell = item.horizontalHeaderItem(column)
                if cell is not None:
                    inspect(cell.text(), "%s.header[%d]" % (name, column))
        if isinstance(item, QTreeWidget):
            header = item.headerItem()
            if header is not None:
                for column in range(header.columnCount()):
                    inspect(header.text(column),
                            "%s.header[%d]" % (name, column))
    inspect(window.statusBar().currentMessage(), "statusBar")

    accepted_parameters = translate("已接受参数")
    catalog_ready = (
        translate("文件") == "File"
        and translate("序列分析") == "Sequence Analysis"
        and accepted_parameters != "已接受参数"
        and not cjk.search(accepted_parameters)
    )

    resources = source_root() / "moire_design_core" / "resources"
    required = (
        "Square_Seed_2L_newtemplate.json",
        "Square_SST_original_128.json",
        "Kagome_Seed_Ka-seed-pore_3L.json",
        "Kagome_SST_original_128.json",
    )
    missing = [name for name in required if not (resources / name).is_file()]
    companion = cadnano_executable()
    report = {
        "application_root": str(application_root()),
        "source_root": str(source_root()),
        "numpy_version": numpy.__version__,
        "pyqt6_imported": bool(PyQt6),
        "designer_engine": str(Path(cadnano2.__file__).resolve()),
        "designer_modules_imported": all((models, structure, template,
                                           analysis_bulk)),
        "worker_count": len(WORKER_MODULES),
        "missing_resources": missing,
        "cadnano_companion": str(companion),
        "cadnano_companion_present": companion.is_file(),
        "tesseract": tool_executable("tesseract"),
        "main_window_initialized": window_initialized,
        "analysis_module_count": window.analysis_module_stack.count(),
        "english_catalog_ready": catalog_ready,
        "visible_cjk_count": len(presentation_text),
        "visible_cjk": presentation_text[:20],
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    window.close()
    window.deleteLater()
    app.processEvents()
    return 0 if (not missing and window_initialized and catalog_ready and
                 not presentation_text) else 2


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--self-test":
        raise SystemExit(self_test())
    if len(sys.argv) >= 3 and sys.argv[1] == "--moire-worker":
        raise SystemExit(dispatch_worker(sys.argv[2], sys.argv[3:]))
    from moire_designer.app import main
    raise SystemExit(main())

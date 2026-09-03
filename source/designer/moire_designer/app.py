"""Standalone application entry point."""

from __future__ import annotations

import sys

from PyQt6.QtWidgets import QApplication

from .mainwindow import MoireDesignerWindow


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    application = QApplication.instance() or QApplication(argv)
    application.setApplicationName("DNA Moiré Designer")
    application.setOrganizationName("cadnano_xinxin")
    application.setQuitOnLastWindowClosed(True)
    window = MoireDesignerWindow()
    window.show()
    return application.exec()


if __name__ == "__main__":
    raise SystemExit(main())

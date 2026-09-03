"""Detached Curved Design report renderer."""

import json
import os
import sys


def main(arguments=None):
    arguments = list(sys.argv[1:] if arguments is None else arguments)
    if len(arguments) != 3:
        return 2
    source_path, output_path, summary_json = arguments
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    # Select the same PyQt6 binding as the desktop application before
    # cadnano's compatibility wrapper probes any legacy Qt installations.
    from PyQt6.QtWidgets import QApplication
    application = QApplication.instance() or QApplication([])
    from cadnano2.model.io.curvedreport import curved_report_data
    from cadnano2.model.io.indelanalysis import (
        write_generated_single_helix_distribution_csv,
        write_pair_curvature_csv, write_pair_curvature_svg)
    from cadnano2.views.curvedreport import create_curved_report_image

    report_data = curved_report_data(source_path)
    summary_lines = json.loads(summary_json)
    temporary_path = output_path + ".writing"
    try:
        create_curved_report_image(
            summary_lines, report_data, temporary_path)
        os.replace(temporary_path, output_path)
        stem, unused_extension = os.path.splitext(output_path)
        write_pair_curvature_csv(
            report_data.get("pair_curvature_rows", []),
            stem + "_pair_curvature.csv")
        write_pair_curvature_svg(
            report_data.get("pair_curvature_rows", []),
            report_data.get("pair_curvature_summary", {}),
            stem + "_pair_curvature.svg")
        write_generated_single_helix_distribution_csv(
            report_data.get("single_helix_distribution", []),
            stem + "_single_helix_distribution.csv")
    finally:
        if os.path.exists(temporary_path):
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
    # Keep the QApplication reference alive until all painting and PNG
    # compression has completed.
    del application
    return 0


if __name__ == "__main__":
    sys.exit(main())

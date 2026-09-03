"""Run the existing orthogonal-sequence designer inside DNA Moiré Designer.

The parameter dialog and generation/export flow intentionally mirror the
existing implementation.  This module owns the launch path so opening the
tool never requires a cadnano document controller or cadnano main window.
"""

from __future__ import annotations

import os
import json
import time
import traceback
from pathlib import Path

from PyQt6.QtCore import QDir, Qt
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QMessageBox,
    QProgressDialog,
)

from cadnano2.data.dnasequences import sequences as standard_scaffold_sequences
from cadnano2.model.io.orthogonalseq import (
    DEFAULT_SETTINGS as ORTHOGONAL_SEQUENCE_DEFAULTS,
    GenerationCancelled,
    generate_sequences,
    normalized_settings,
    read_sequence_text,
    write_orthogonal_workbook,
)
from cadnano2.views.orthogonalsequences import OrthogonalSequenceDialog
import cadnano2.views.orthogonalsequences as orthogonal_dialog_module
from cadnano2.views.primer3analysis import Primer3AnalysisDialog
from .i18n import localize_xlsx


def generate_orthogonal_sequences_automatic(length_counts, output_directory,
                                            parent=None):
    """Generate mutually orthogonal sequences for one or more exact lengths.

    Earlier length groups are supplied as background to every later group,
    so the complete returned set remains cross-length orthogonal.  A full
    diagnostic workbook is written for every distinct requested length.
    """
    output_root = Path(output_directory).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    generated = []
    by_length = {}
    reports = []
    background = []
    total = sum(int(count) for count in length_counts.values())
    completed_count = 0
    last_ui_update = [0.0]
    # Design the longest inputs first.  If hundreds of short inputs are fixed
    # before a small number of longer inputs, their combined 8-mer coverage
    # can exhaust the long-sequence search space even though a valid joint set
    # exists.  Longest-first preserves every hard orthogonality threshold and
    # leaves the much larger short-sequence candidate space for the final
    # groups.  Single-length workflows are unchanged.
    for length, count in sorted(
            ((int(length), int(count)) for length, count in
             length_counts.items()), key=lambda item: item[0], reverse=True):
        settings = dict(ORTHOGONAL_SEQUENCE_DEFAULTS)
        settings.update({
            "length": length,
            "count": count,
            "gc_min": 0.40,
            "gc_max": 0.60,
            "max_same_substring": 7,
            "max_cross_complement": 7,
        })
        settings = normalized_settings(settings)

        def progress(accepted, attempts):
            if parent is None:
                return
            now = time.monotonic()
            if accepted < count and now - last_ui_update[0] < 0.10:
                return
            last_ui_update[0] = now
            try:
                parent.statusBar().showMessage(
                    "Designing SST sublattice input sequences: %d of %d; "
                    "current length: %d nt; candidates evaluated: %d" %
                    (completed_count + accepted, total, length, attempts))
                # Processing every candidate made the click callback re-enter
                # thousands of times. A throttled pass keeps macOS responsive
                # without exposing the slot to repeated button activation.
                QApplication.processEvents()
            except RuntimeError:
                # The window may be closing while generation unwinds. This is
                # not a sequence-generation failure and must not abort PyQt.
                return

        result = generate_sequences(
            settings,
            background_sequences=list(background),
            scaffold_sequences=(), progress=progress,
            cancelled=(lambda: False))
        result["input_file"] = ""
        if not result.get("complete") or len(result.get("sequences", ())) != count:
            raise RuntimeError(
                "Only %d of %d requested %d-nt sequences were generated. "
                "Use Expert Mode to adjust the design constraints." %
                (len(result.get("sequences", ())), count, length))
        workbook = output_root / (
            "orthogonal_input_%dnt_analysis.xlsx" % length)
        write_orthogonal_workbook(str(workbook), result)
        localize_xlsx(workbook)
        values = [str(item).upper() for item in result["sequences"]]
        generated.extend(values)
        by_length[length] = values
        background.extend(values)
        completed_count += len(values)
        reports.append({
            "length": length,
            "count": count,
            "workbook": str(workbook),
            "attempts": int(result.get("attempts", 0)),
            "rejections": dict(result.get("rejections", {})),
        })
    summary = {"sequences": generated, "by_length": by_length,
               "reports": reports,
               "settings": {"gc_min": 0.40, "gc_max": 0.60,
                            "max_same_substring": 7,
                            "max_cross_complement": 7}}
    manifest = output_root / "automatic_input_sequence_manifest.json"
    manifest.write_text(json.dumps(
        summary, ensure_ascii=False, indent=2), encoding="utf-8")
    summary["manifest"] = str(manifest)
    return summary


def rejection_details(result):
    """Format rejection counters exactly as the existing designer does."""
    labels = {
        "global_gc": "Global GC content",
        "local_gc": "Local GC content",
        "homopolymer": "Homopolymer run",
        "entropy": "Low sequence complexity",
        "self_complement": "Self-complementarity",
        "hairpin": "Hairpin formation",
        "forbidden_motif": "Forbidden motif",
        "same_substring": "Same-orientation exact match",
        "cross_complement": "Interstrand complementarity",
        "hamming": "Hamming distance",
    }
    lines = [
        "Candidates evaluated: %d" % result.get("attempts", 0),
        "Background sequences screened: %d"
        % result.get("background_count", 0),
        "",
        "Candidates rejected by criterion:",
    ]
    for key, count in sorted(
        result.get("rejections", {}).items(),
        key=lambda item: (-item[1], item[0]),
    ):
        lines.append("%s: %d" % (labels.get(key, key), count))
    return "\n".join(lines)


def _show_message(parent, icon, title, text, details=None):
    message = QMessageBox(
        icon,
        title,
        text,
        QMessageBox.StandardButton.Ok,
        parent,
    )
    if details:
        message.setDetailedText(details)
    message.exec()


def run_orthogonal_sequence_designer(
    parent,
    project_filename=None,
    primer3_entries=(),
    suggested_directory=None,
):
    """Open and run the unchanged designer flow as a child of ``parent``.

    Returns a dictionary containing the exported path and Primer3 entries on
    success.  Cancellation and validation failures return ``None``.
    """
    icon_directory = (
        Path(orthogonal_dialog_module.__file__).resolve().parents[1]
        / "ui"
        / "mainwindow"
        / "images"
    )
    QDir.addSearchPath("icons", str(icon_directory))

    directory = (str(Path(suggested_directory).expanduser().resolve())
                 if suggested_directory else
                 (os.path.expanduser("~/Desktop")
                  if not project_filename else
                  str(Path(project_filename).resolve().parent)))
    Path(directory).mkdir(parents=True, exist_ok=True)
    dialog = OrthogonalSequenceDialog(
        ORTHOGONAL_SEQUENCE_DEFAULTS,
        parent,
        primer3_entries=primer3_entries,
        suggested_directory=directory,
    )
    if dialog.exec() != QDialog.DialogCode.Accepted:
        return None

    requested = dialog.settings()
    input_filename = requested.pop("input_file", "")
    scaffold_name = requested.pop("scaffold_name", "")
    try:
        settings = normalized_settings(requested)
    except ValueError as error:
        _show_message(
            parent, QMessageBox.Icon.Warning, "Orthogonal Sequence Design",
            str(error)
        )
        return None

    input_sequences = []
    input_errors = []
    if input_filename:
        try:
            input_sequences, input_errors = read_sequence_text(input_filename)
        except (IOError, OSError, UnicodeError) as error:
            _show_message(
                parent,
                QMessageBox.Icon.Critical,
                "Orthogonal Sequence Input Error",
                "The selected TXT file could not be read.",
                str(error),
            )
            return None
    if input_errors:
        _show_message(
            parent,
            QMessageBox.Icon.Critical,
            "Orthogonal Sequence Input Error",
            "%d invalid entries were found in the TXT file. Sequence "
            "generation was not started." % len(input_errors),
            "\n".join(input_errors),
        )
        return None

    scaffold_background = []
    if scaffold_name:
        scaffold_sequence = standard_scaffold_sequences.get(scaffold_name)
        if not scaffold_sequence:
            _show_message(
                parent,
                QMessageBox.Icon.Critical,
                "Orthogonal Sequence Design",
                "The selected scaffold sequence could not be found: %s"
                % scaffold_name,
            )
            return None
        scaffold_background.append((scaffold_name, scaffold_sequence))

    design_name = (
        "cadnano"
        if not project_filename
        else os.path.splitext(os.path.basename(str(project_filename)))[0]
    )
    suggested = os.path.join(
        directory, design_name + "_orthogonal_sequences.xlsx"
    )
    filename = QFileDialog.getSaveFileName(
        parent,
        "Orthogonal Sequence Design — Save Results",
        suggested,
        "Orthogonal Sequence Workbook (*.xlsx)",
    )
    if isinstance(filename, (tuple, list)):
        filename = filename[0]
    if not filename or os.path.isdir(filename):
        return None
    if not filename.lower().endswith(".xlsx"):
        filename += ".xlsx"

    progress_dialog = QProgressDialog(
        "Generating and evaluating candidate sequences…", "Cancel", 0, 0,
        parent
    )
    progress_dialog.setWindowTitle("Orthogonal Sequence Design")
    progress_dialog.setWindowModality(Qt.WindowModality.WindowModal)
    progress_dialog.setMinimumDuration(0)

    def update_progress(accepted, attempts):
        progress_dialog.setLabelText(
            "Generated %d of %d sequences; evaluated %d candidates…"
            % (accepted, settings["count"], attempts)
        )
        QApplication.processEvents()

    try:
        result = generate_sequences(
            settings,
            background_sequences=input_sequences,
            scaffold_sequences=scaffold_background,
            progress=update_progress,
            cancelled=progress_dialog.wasCanceled,
        )
        result["input_file"] = input_filename
    except GenerationCancelled:
        progress_dialog.close()
        parent.statusBar().showMessage(
            "Orthogonal-sequence generation was canceled.", 5000)
        return None
    except Exception as error:
        progress_dialog.close()
        traceback.print_exc()
        _show_message(
            parent,
            QMessageBox.Icon.Critical,
            "Orthogonal Sequence Design Failed",
            "An error occurred while generating the sequences.",
            str(error),
        )
        return None
    finally:
        progress_dialog.close()

    if not result["sequences"]:
        _show_message(
            parent,
            QMessageBox.Icon.Warning,
            "Orthogonal Sequence Design",
            "No sequence satisfying the current constraints was found "
            "within the maximum number of attempts. Relax the constraints, "
            "particularly the GC-content, same-orientation exact-match, or "
            "interstrand-complementarity limits. If advanced criteria are "
            "enabled, disable or relax them individually.",
            rejection_details(result),
        )
        return None
    if not result["complete"]:
        choice = QMessageBox.question(
            parent,
            "Orthogonal Sequence Design Incomplete",
            "%d sequences were requested, but only %d passed all criteria "
            "within the permitted number of attempts. Save the sequences "
            "that passed?"
            % (settings["count"], len(result["sequences"])),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if choice != QMessageBox.StandardButton.Yes:
            return None
    try:
        write_orthogonal_workbook(filename, result)
        localize_xlsx(filename)
    except (IOError, OSError, ValueError) as error:
        _show_message(
            parent,
            QMessageBox.Icon.Critical,
            "Orthogonal Sequence Export Failed",
            "The XLSX workbook could not be written.",
            str(error),
        )
        return None

    next_primer3_entries = [
        ("Input", "Input-%03d" % index, sequence)
        for index, sequence in enumerate(result.get("input_sequences", ()), 1)
    ]
    next_primer3_entries.extend(
        ("Newly generated", "New sequence-%03d" % index, sequence)
        for index, sequence in enumerate(result.get("sequences", ()), 1)
    )
    completion = QMessageBox(
        QMessageBox.Icon.Information,
        "Orthogonal Sequence Design Complete",
        "Read %d input sequences, used %d scaffold sequences, and generated "
        "%d new sequences. The workbook was saved to:\n%s\n\n"
        "The Sequence Analysis worksheet distinguishes input, scaffold, and "
        "newly generated sequences. Scaffold sequences are used for "
        "screening only and are excluded from Pairwise Analysis."
        % (
            len(result.get("input_sequences", ())),
            len(result.get("scaffold_sequences", ())),
            len(result["sequences"]),
            filename,
        ),
        QMessageBox.StandardButton.Ok,
        parent,
    )
    completion.setDetailedText(rejection_details(result))
    primer3_button = completion.addButton(
        "Primer3 Thermodynamic Analysis…", QMessageBox.ButtonRole.ActionRole
    )
    completion.exec()
    if completion.clickedButton() is primer3_button:
        Primer3AnalysisDialog(next_primer3_entries, directory, parent).exec()

    return {
        "filename": filename,
        "primer3_entries": next_primer3_entries,
        "result": result,
    }

"""Subprocess-safe API for the staged sequence assignment workflow."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable

from moire_runtime import worker_command


def _run(mode: str, *arguments: str) -> Dict[str, Any]:
    values = [str(item) for item in arguments]

    # CreateProcess receives one flat command-line string on Windows.  A full
    # set of SST assignments can readily exceed its limit when JSON is passed
    # as a positional argument (WinError 206: filename or extension too long).
    # Use a short response file for large requests on every platform so the
    # packaged and development runtimes exercise exactly the same path.
    command_size = sum(len(item) + 3 for item in values)
    with tempfile.TemporaryDirectory(prefix="dmd_seq_") as directory:
        if command_size > 8000:
            arguments_file = Path(directory) / "args.json"
            arguments_file.write_text(
                json.dumps(values, ensure_ascii=False), encoding="utf-8")
            command_values = ["@arguments-file", str(arguments_file)]
        else:
            command_values = values
        command = worker_command(
            "sequence-workflow", mode, *command_values)
        completed = subprocess.run(
            command, check=False, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(detail or "The sequence workflow did not finish.")
    try:
        return json.loads(completed.stdout)
    except Exception as error:
        raise RuntimeError(
            "The sequence workflow did not return a valid report.") from error


def analyze_sequence_design(filename: str) -> Dict[str, Any]:
    return _run("analyze", filename)


def extract_scaffold_sequence(filename: str,
                              target: Dict[str, Any]) -> Dict[str, Any]:
    return _run("extract-scaffold", filename,
                json.dumps(target, ensure_ascii=False))


def list_standard_scaffolds(target_length: int, multiple: bool = False,
                            used_names: Iterable[str] = ()) \
        -> Dict[str, Any]:
    return _run("list-scaffolds", str(int(target_length)),
                "1" if multiple else "0",
                json.dumps(list(used_names), ensure_ascii=False))


def assign_standard_scaffold_sequence(
        target: Dict[str, Any], name: str, multiple: bool = False,
        used_names: Iterable[str] = ()) -> Dict[str, Any]:
    return _run("assign-scaffold",
                json.dumps(target, ensure_ascii=False), name,
                "1" if multiple else "0",
                json.dumps(list(used_names), ensure_ascii=False))


def export_sst_input_template(design: str, filename: str,
                              layers_identical: bool) -> Dict[str, Any]:
    return _run("export-input-template", design, filename,
                "1" if layers_identical else "0")


def import_sst_input_template(design: str, filename: str,
                              layers_identical: bool) -> Dict[str, Any]:
    return _run("import-input-template", design, filename,
                "1" if layers_identical else "0")


def build_sequenced_design(design: str, filename: str,
                           assignments: Iterable[Dict[str, Any]]) \
        -> Dict[str, Any]:
    return _run("build-sequenced", design, filename,
                json.dumps(list(assignments), ensure_ascii=False))


def export_final_package(project_file: str, sequenced_design: str,
                         output_directory: str,
                         workflow: Dict[str, Any]) -> Dict[str, Any]:
    return _run("final-export", project_file or "", sequenced_design,
                output_directory, json.dumps(workflow, ensure_ascii=False))

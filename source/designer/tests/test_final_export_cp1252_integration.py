"""Real final-export regression for legacy Windows stdout encodings."""

from __future__ import annotations

import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from moire_design_core import (
    finalize_structure,
    generate_scaffold_review,
    write_shifted_sst,
)
from moire_design_core import sequence_workflow_worker
from moire_design_core.sequence_workflow_worker import analyze, build_sequenced


class FinalExportCp1252IntegrationTests(unittest.TestCase):
    def test_real_zero_spacing_final_export_survives_cp1252_stdout(self):
        """A real 0-bp auxiliary design must export and report on Windows."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sst = Path(write_shifted_sst(
                str(root / "sst.json"), 64, 0, 64, 128, 128,
                lattice_type="square_kagome"))
            scaffold = root / "scaffold.json"
            final = root / "final.json"
            sequenced = root / "with_sequences.json"
            generate_scaffold_review(str(scaffold), str(sst))
            finalize_structure(str(scaffold), str(final))

            report = analyze(str(final))
            assignments = []
            for group in report["targets"].values():
                for target in group:
                    sequence = ("ACGT" * (
                        (int(target["length"]) + 3) // 4))[
                            :int(target["length"])]
                    assignments.append({
                        "target_id": target["id"],
                        "sequence": sequence,
                    })
            build_report = build_sequenced(
                str(final), str(sequenced), assignments)
            self.assertEqual(build_report["unresolved_output_bases"], 0)

            workflow = {
                "automatic_design_exports": {"sst": str(sst)},
                "scaffold_accepted": str(scaffold),
                "structure_complete": str(final),
                "sequence_assignments": assignments,
            }
            argv = [
                "sequence_workflow_worker.py", "final-export", "",
                str(sequenced), str(root / "export"),
                json.dumps(workflow, ensure_ascii=False),
            ]
            raw = io.BytesIO()
            stdout = io.TextIOWrapper(
                raw, encoding="cp1252", errors="strict", newline="\n")
            try:
                with patch.object(sys, "argv", argv), \
                        patch.object(sys, "stdout", stdout):
                    sequence_workflow_worker.main()
                    stdout.flush()
                transported = raw.getvalue().decode("cp1252")
            finally:
                stdout.detach()

            result = json.loads(transported)
            export_root = Path(result["root"])
            self.assertTrue(export_root.is_dir())
            self.assertTrue((export_root / "PDB∕oxView files").is_dir())
            self.assertTrue((export_root / "Oligonucleotide sequences" /
                             "all_sequences.xlsx").is_file())
            self.assertTrue((export_root / "caDNAno design files" /
                             "sst_scaffold_staple_capture_no_sequence.json").is_file())
            self.assertIn(r"\u2215", transported)


if __name__ == "__main__":
    unittest.main()

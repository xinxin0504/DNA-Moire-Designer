from contextlib import redirect_stdout
import io
import json
import math
import sys
import unittest
from unittest.mock import patch

try:
    import numpy as np
except ImportError:  # The production Qt venv delegates this worker externally.
    np = None

if np is not None:
    from moire_design_core import image_analysis_worker as image_worker
    from moire_design_core.image_analysis_worker import (
        _axis_fingerprint_alignment, _welch_lattice_fft, analyze_fft,
        analyze_tem, detect_scale_bar)
else:
    image_worker = None
    _axis_fingerprint_alignment = _welch_lattice_fft = analyze_fft = None
    analyze_tem = detect_scale_bar = None


@unittest.skipIf(np is None, "NumPy is provided by the external image worker")
class ImageAnalysisTests(unittest.TestCase):
    def test_cli_accepts_complete_moire_analysis_contract(self):
        mocked_result = {
            "analysis_kind": "single",
            "period_candidates_px": [],
        }
        argv = [
            "image_analysis_worker.py", "tem", "synthetic.pgm",
            "--original", "synthetic.tif",
            "--analysis-kind", "single",
            "--output-dir", "analysis-output",
            "--theoretical-symmetry", "Square",
            "--theoretical-a-nm", "2.0",
            "--pixel-size-nm", "0.32786885245901637",
        ]
        with patch.object(sys, "argv", argv), \
                patch.object(image_worker, "read_pgm",
                             return_value=np.zeros((8, 8), dtype=float)), \
                patch.object(image_worker, "detect_scale_bar",
                             return_value=None), \
                patch.object(image_worker, "ocr_scale",
                             return_value=(None, "")), \
                patch.object(image_worker, "analyze_tem",
                             return_value=mocked_result) as analyze:
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                image_worker.main()

        payload = json.loads(stdout.getvalue())
        self.assertEqual(payload["analysis_kind"], "single")
        analyze.assert_called_once()
        arguments = analyze.call_args.args
        self.assertEqual(arguments[3], "single")
        self.assertEqual(arguments[4], "analysis-output")
        self.assertEqual(arguments[6], [2.0])
        self.assertEqual(arguments[7], ["Square"])

    def test_scale_bar_detection(self):
        image = np.ones((400, 600), dtype=float) * 128
        image[350:356, 350:510] = 255
        bar = detect_scale_bar(image)
        self.assertIsNotNone(bar)
        self.assertAlmostEqual(bar["pixel_length"], 160, delta=2)

    def test_fft_two_square_orientations(self):
        size = 512
        image = np.zeros((size, size), dtype=float)
        center = size / 2
        for offset in (0, 8):
            for index in range(4):
                angle = math.radians(offset + 90 * index)
                x = int(round(center + 120 * math.cos(angle)))
                y = int(round(center - 120 * math.sin(angle)))
                image[y - 2:y + 3, x - 2:x + 3] = 255
        result = analyze_fft(image)
        self.assertAlmostEqual(result["twist_angle_deg"], 8, delta=1.5)
        self.assertGreaterEqual(result["peak_count"], 4)

    def test_kagome_fingerprint_alignment_preserves_non_60_degree_angles(self):
        first = [10.0, 66.0, 130.0]  # gaps 56, 64, 60 degrees
        second = [18.0, 74.0, 138.0]
        aligned = _axis_fingerprint_alignment(first, second, 60.0)
        self.assertIsNotNone(aligned)
        self.assertAlmostEqual(aligned["twist_angle_deg"], 8.0, delta=.01)
        self.assertAlmostEqual(aligned["maximum_residual_deg"], 0.0,
                               delta=.01)

    def test_tem_period_candidates(self):
        yy, xx = np.mgrid[:512, :512]
        image = (128 + 35 * np.cos(2 * np.pi * xx / 64) +
                 25 * np.cos(2 * np.pi * yy / 64) +
                 15 * np.cos(2 * np.pi * (xx + yy) / 16))
        result = analyze_tem(image)
        self.assertAlmostEqual(result["period_candidates_px"][0], 64,
                               delta=4)
        # A coarse low-frequency candidate is not promoted to a reported
        # bilayer period without a validated two-layer first-order FFT pair.
        self.assertIsNone(result["moire_period_px"])
        self.assertFalse(result["fft_twist_reliable"])
        self.assertIsNotNone(result["lattice_constant_px"])

    def test_tem_fft_resolves_two_square_first_order_families(self):
        size = 512
        yy, xx = np.mgrid[:size, :size]
        lattice = 10.0
        twist = math.radians(4.0)
        rotated_x = xx * math.cos(twist) + yy * math.sin(twist)
        rotated_y = -xx * math.sin(twist) + yy * math.cos(twist)
        image = (128 + 18 * np.cos(2 * np.pi * xx / lattice) +
                 18 * np.cos(2 * np.pi * yy / lattice) +
                 18 * np.cos(2 * np.pi * rotated_x / lattice) +
                 18 * np.cos(2 * np.pi * rotated_y / lattice))
        result = _welch_lattice_fft(image)
        self.assertTrue(result["valid"])
        self.assertAlmostEqual(result["lattice_constant_px"], lattice,
                               delta=0.25)
        self.assertAlmostEqual(result["twist_angle_deg"], 4.0, delta=0.25)

if __name__ == "__main__":
    unittest.main()

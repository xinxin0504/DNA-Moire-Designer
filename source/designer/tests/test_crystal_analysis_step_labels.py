import ast
import unittest
from pathlib import Path


SOURCE_PATH = Path(__file__).parents[1] / "moire_designer" / "analysis_bulk.py"


class CrystalAnalysisStepLabelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        module = ast.parse(SOURCE_PATH.read_text(encoding="utf-8"))
        cls.functions = {
            node.name: ast.unparse(node)
            for node in ast.walk(module)
            if isinstance(node, ast.FunctionDef)
        }

    def test_analysis_button_keeps_step_three_in_single_and_bulk_modes(self):
        source = self.functions["_bulk_analysis_changed"]
        self.assertIn("3. Automatically Identify and Analyze", source)
        self.assertIn("3. Run Batch Detection and Analysis", source)
        self.assertIn("Run Batch Analysis", source)
        self.assertIn("Run Analysis", source)

    def test_crystal_workflow_moves_scale_to_upload_and_keeps_results_visible(self):
        source = self.functions["_build_crystal_analysis_panel"]
        self.assertIn("self.crystal_section_group = QButtonGroup(self)", source)
        self.assertIn("self.crystal_section_group.setExclusive(False)", source)
        self.assertIn("button.setObjectName('parameterStepButton')", source)
        self.assertIn("button.toggled.connect(content.setVisible)", source)
        self.assertIn("button.setChecked(False)", source)
        self.assertIn("content.hide()", source)
        for title in (
            "1. 上传 TEM 图像",
            "2. 选择单层或双层分析",
        ):
            self.assertIn(repr(title), source)
        self.assertIn(
            "self.analysis_calibration_box.setTitle(translate('4. 结果摘要'))",
            source,
        )
        self.assertIn("self.analysis_calibration_box.hide()", source)
        self.assertIn("Scale bar pixel length", source)
        self.assertIn("Scale bar label value", source)
        self.assertIn("self.crystal_run_section = QWidget()", source)
        self.assertIn(
            "self.crystal_run_step_button = crystal_step_button('3. Automatically Identify and Analyze', self.crystal_run_section)",
            source,
        )
        self.assertIn(
            "self.run_image_analysis_button = QPushButton('Run Analysis')",
            source)
        self.assertNotIn("crystal_step_button('3. 自动识别并分析'", source)
        self.assertIn("left.addWidget(calibration_box)", source)
        self.assertNotIn("self.crystal_result_step_button", source)

    def test_step_four_title_switches_with_bulk_mode(self):
        source = self.functions["_bulk_analysis_changed"]
        self.assertIn("self.analysis_calibration_box.setTitle", source)
        self.assertIn(
            "self.analysis_calibration_box.setVisible(bool(enabled))", source)
        self.assertIn("self.single_result_summary.hide()", source)
        self.assertIn("4. 批量运行状态", source)
        self.assertIn("4. 结果摘要", source)

    def test_single_scale_preflight_keeps_analysis_section_openable(self):
        source = self.functions["_update_single_scale_gate"]
        self.assertGreaterEqual(
            source.count("self.crystal_mode_step_button.setEnabled(True)"), 2)
        self.assertNotIn("self.crystal_mode_step_button.setEnabled(ready)", source)
        self.assertIn("self.crystal_run_step_button.setEnabled(True)", source)
        self.assertNotIn("self.crystal_run_step_button.setEnabled(ready)", source)
        self.assertIn("self.run_image_analysis_button.setEnabled(ready)", source)
        self.assertIn("has_image and pixels_valid and value_valid", source)
        self.assertIn("#ffe6e6", source)

    def test_analysis_page_matches_design_inset_and_preview_frames(self):
        tab_source = self.functions["_build_analysis_tab"]
        panel_source = self.functions["_build_crystal_analysis_panel"]
        self.assertIn("outer.setContentsMargins(0, 8, 0, 8)", tab_source)
        self.assertIn("panel = QGroupBox(title)", panel_source)
        self.assertIn("panel.setObjectName('previewGroupBox')", panel_source)
        self.assertNotIn("panel.setObjectName('analysisFigurePanel')",
                         panel_source)
        for title in ("Original TEM", "Reconstructed lattice", "Inverse FFT"):
            self.assertIn(repr(title), panel_source)

    def test_fft_references_are_independent_and_english(self):
        build = self.functions["_build_theoretical_reference_widget"]
        values = self.functions["_theoretical_reference_values"]
        worker = self.functions["_run_image_worker"]
        self.assertIn("FFT Recognition References (Optional)", build)
        self.assertIn("Expected symmetries (maximum 2)", build)
        self.assertIn("Expected lattice constants a (multiple values)", build)
        self.assertNotIn("do not correspond row by row", build)
        self.assertNotIn("对称性", build)
        self.assertIn("'symmetries'", values)
        self.assertIn("'a_nm'", values)
        self.assertIn("for symmetry in", worker)
        self.assertIn("for a_nm in", worker)

    def test_upload_runs_scale_preflight_and_immediately_sets_preview(self):
        select_source = self.functions["select_tem_image"]
        detect_source = self.functions["_detect_single_image_scale"]
        preview_source = self.functions["_show_scale_preflight_preview"]
        self.assertIn("self.tem_analysis_image.set_image(source_image)",
                      select_source)
        self.assertIn("self._detect_single_image_scale()", select_source)
        self.assertIn("self._run_scale_detection(source)", detect_source)
        self.assertIn("Scale bar detected · value not recognized",
                      preview_source)
        self.assertIn("font.setPixelSize(14)", preview_source)
        self.assertNotIn("int(15 * scale)", preview_source)
        self.assertIn("self._normalize_annotation_dpi(image)",
                      preview_source)
        self.assertIn("rect.right() + gap + text_width", preview_source)
        self.assertNotIn("rect.top() - text_height", preview_source)

    def test_analysis_standardizer_preserves_step_button_style(self):
        source = self.functions["_standardize_analysis_button_heights"]
        self.assertIn("'parameterStepButton'", source)

    def test_bulk_to_single_intro_does_not_read_analysis_metrics(self):
        source = self.functions["_update_crystal_intro"]
        self.assertNotIn("metrics", source)
        self.assertIn("if self._bulk_enabled()", source)
        self.assertIn("else:", source)

    def test_bulk_toggle_preserves_independent_single_image_state(self):
        changed = self.functions["_bulk_analysis_changed"]
        capture = self.functions["_capture_analysis_mode_state"]
        restore = self.functions["_restore_analysis_mode_state"]
        self.assertIn("self._single_analysis_state = current_state", changed)
        self.assertIn("self._bulk_analysis_state = current_state", changed)
        self.assertIn("self._restore_analysis_mode_state", changed)
        self.assertNotIn("self._tem_image_paths = []", changed)
        for key in ("paths", "path", "records", "preflight",
                    "scale_pixels", "scale_nm"):
            self.assertIn(repr(key), capture)
        self.assertIn("self._show_scale_preflight_preview()", restore)

    def test_crystal_sidebar_uses_content_fitted_compact_width(self):
        source = self.functions["_build_crystal_analysis_panel"]
        self.assertIn("self._analysis_sidebar_width(", source)
        self.assertIn("minimum=390", source)
        self.assertNotIn("splitter.setSizes([430, 1080])", source)

    def test_crystal_sidebar_width_helper_is_packaged_with_analysis(self):
        self.assertIn("_analysis_sidebar_width", self.functions)
        source = self.functions["_analysis_sidebar_width"]
        self.assertIn("widget.sizeHint().width()", source)
        self.assertIn("max(int(minimum), min(int(maximum)", source)

    def test_bulk_scale_help_has_room_for_wrapped_text(self):
        build_source = self.functions["_build_crystal_analysis_panel"]
        mode_source = self.functions["_bulk_scale_mode_changed"]
        self.assertIn("self.bulk_scale_help.setMinimumHeight(44)", build_source)
        self.assertIn(
            "QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum",
            build_source,
        )
        self.assertIn("Use raw TIFF files with embedded physical pixel size", mode_source)

    def test_bilayer_original_restores_moire_cells_and_period_annotation(self):
        paint = self.functions["_paint_annotations"]
        cells = self.functions["_draw_moire_cells"]
        self.assertIn(
            "self._draw_moire_cells(painter, record, draw_cells=True)",
            paint,
        )
        self.assertIn("TEM-derived moiré period = %.3f nm", cells)
        self.assertIn("basis_vectors_px", cells)
        self.assertIn("representative_pair", cells)

    def test_fft_first_order_frames_are_cyclic_closed_dashed_and_unfilled(self):
        source = self.functions["_draw_fft_inset"]
        self.assertIn("ordered = sorted(polygon", source)
        self.assertIn("math.atan2", source)
        self.assertIn("Qt.PenStyle.DashLine", source)
        self.assertIn("outline.closeSubpath()", source)
        self.assertIn("painter.setBrush(Qt.BrushStyle.NoBrush)", source)
        self.assertIn("painter.drawPath(outline)", source)
        self.assertNotIn("label_box", source)

    def test_reconstructed_same_symmetry_bilayer_reports_one_lattice_a(self):
        source = self.functions["_draw_reconstruction_annotations"]
        self.assertIn("shared_bilayer_a", source)
        self.assertIn("len(entry_symmetries) == 1", source)
        self.assertIn("len(entries) > 1 and (not shared_bilayer_a)", source)
        self.assertIn("'a = %s'", source)

    def test_raster_annotation_dpi_is_normalized_before_painting(self):
        normalize = self.functions["_normalize_annotation_dpi"]
        font = self.functions["_annotation_font"]
        render = self.functions["_render_analysis"]
        self.assertIn("setDotsPerMeterX", normalize)
        self.assertIn("setDotsPerMeterY", normalize)
        self.assertIn("font.setPixelSize(int(pixel_size))", font)
        self.assertIn("self._normalize_annotation_dpi(raster)", render)


if __name__ == "__main__":
    unittest.main()

import ast
import os
import unittest
from pathlib import Path
from unittest.mock import patch


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QApplication, QFrame, QLabel, QMessageBox, QSizePolicy)
    from moire_designer.mainwindow import MoireDesignerWindow
except Exception:  # pragma: no cover - core-only test environments
    QApplication = None
    QMessageBox = None
    QLabel = None
    QFrame = None
    QSizePolicy = None
    Qt = None
    MoireDesignerWindow = None


class MainWindowWorkflowSourceTests(unittest.TestCase):
    def test_parameter_changes_reset_all_downstream_design_and_sequence_state(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        module = ast.parse(source_path.read_text(encoding="utf-8"))
        window_class = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef)
            and node.name == "MoireDesignerWindow")
        functions = {
            node.name: ast.unparse(node)
            for node in window_class.body
            if isinstance(node, ast.FunctionDef)
        }

        reset = functions["_invalidate_downstream_design_state"]
        for key in (
                "sst_two_layer", "scaffold_review", "scaffold_accepted",
                "structure_complete", "structure_accepted",
                "automatic_design_exports", "sequence_analysis",
                "sequence_assignments", "sequence_scaffold_accepted",
                "sequence_sst_detected", "sequence_sst_accepted",
                "sequence_source", "sequence_exports"):
            self.assertIn(repr(key), reset)
        self.assertIn("self._sequence_analysis = None", reset)
        self.assertIn("self._sequence_assignments = {}", reset)

        for name in ("_invalidate_design_basis_acceptance",
                     "_invalidate_parameter_acceptance"):
            self.assertIn(
                "self._invalidate_downstream_design_state(workflow)",
                functions[name])

        restore = functions["_restore_structure_workflow"]
        self.assertLess(
            restore.index("workflow.get('settings_signature')"),
            restore.index("basis_ok = bool"))
        self.assertIn(
            "self._invalidate_downstream_design_state(workflow)", restore)
        self.assertIn(
            "self._clear_sequence_cards(self.scaffold_cards_layout)",
            restore)
        self.assertIn(
            "self._clear_sequence_cards(self.sst_cards_layout)", restore)

    def test_mixed_lattice_consistency_does_not_copy_input_sequences(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        module = ast.parse(source_path.read_text(encoding="utf-8"))
        window_class = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef)
            and node.name == "MoireDesignerWindow")
        function = next(
            node for node in window_class.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_sst_layers_identical")
        source = ast.unparse(function)
        self.assertIn(
            "self.project.settings.lattice_symmetry != 'square_kagome'",
            source)

    def test_moire_title_and_red_blue_parameter_display_order(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn('"1 Moiré 参数输入"', source)
        self.assertLess(
            source.index("sst_form.addRow(self.sst_z3_label"),
            source.index("sst_form.addRow(self.sst_spacing_label"))
        self.assertLess(
            source.index("sst_form.addRow(self.sst_spacing_label"),
            source.index("sst_form.addRow(self.sst_z1_label"))
        self.assertLess(
            source.index("seed_form.addWidget(self.seed_z3_label, 0, 0)"),
            source.index("seed_form.addWidget(self.seed_z2_label, 1, 0)"))
        self.assertLess(
            source.index("seed_form.addWidget(self.seed_z2_label, 1, 0)"),
            source.index("seed_form.addWidget(self.seed_z1_label, 2, 0)"))

    def test_startup_is_non_modal_and_mode_switch_is_in_workflow_bar(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertNotIn(
            "QTimer.singleShot(0, self._prompt_startup)", source)
        self.assertIn(
            'self.mode_switch_button = QPushButton('
            '"Switch to Analysis Mode")', source)
        self.assertIn(
            'self.mode_switch_button.setText("Switch to Design Mode")',
            source)
        self.assertIn(
            '"Analysis Mode. No project file is required or saved."',
            source)

        module = ast.parse(source)
        window_class = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef)
            and node.name == "MoireDesignerWindow")
        prompt = next(
            node for node in window_class.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_prompt_startup")
        self.assertNotIn("QMessageBox", ast.unparse(prompt))

    def test_analysis_workspace_contains_only_moire_analysis(self):
        mainwindow_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        mainwindow_source = mainwindow_path.read_text(encoding="utf-8")
        self.assertNotIn('QAction("Gel analysis"', mainwindow_source)
        self.assertNotIn('QAction("Particle analysis"', mainwindow_source)
        self.assertIn('"Moiré analysis", self)', mainwindow_source)

        analysis_path = Path(__file__).parents[1] / "moire_designer" / \
            "analysis_bulk.py"
        analysis_source = analysis_path.read_text(encoding="utf-8")
        module = ast.parse(analysis_source)
        mixin = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef)
            and node.name == "AnalysisBulkMixin")
        build_tab = next(
            node for node in mixin.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_build_analysis_tab")
        build_source = ast.unparse(build_tab)
        self.assertNotIn("_build_gel_analysis_panel", build_source)
        self.assertNotIn("_build_particle_analysis_panel", build_source)
        self.assertIn("_build_crystal_analysis_panel", build_source)

    def test_analysis_buttons_share_the_analyze_button_height(self):
        root = Path(__file__).parents[1] / "moire_designer"
        mainwindow_source = (root / "mainwindow.py").read_text(
            encoding="utf-8")
        analysis_source = (root / "analysis_bulk.py").read_text(
            encoding="utf-8")
        self.assertIn("ANALYSIS_BUTTON_HEIGHT = 34", analysis_source)
        self.assertIn(
            "button.setFixedHeight(self.ANALYSIS_BUTTON_HEIGHT)",
            analysis_source)
        self.assertIn(
            'button.setObjectName("primaryButton")', analysis_source)
        self.assertIn(
            "self._standardize_analysis_button_heights(tab)",
            analysis_source)
        self.assertEqual(
            mainwindow_source.count(
                "setFixedHeight(self.ANALYSIS_BUTTON_HEIGHT)"), 2)
        self.assertIn(
            '"workflowButton", "workflowButton")', mainwindow_source)
        self.assertIn(
            'self.mode_switch_button.setObjectName("primaryButton")',
            mainwindow_source)

    def test_sequence_workflow_uses_collapsible_step_buttons(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        for label in (
                "3.1 Assign scaffold sequences",
                "3.2 Assign SST sublattice input sequences",
                "3.3 Final export"):
            self.assertIn('"%s"' % label, source)
        self.assertIn(
            'button.setObjectName("parameterStepButton")', source)
        self.assertIn(
            'button.toggled.connect(content.setVisible)', source)
        self.assertIn('button.setChecked(False)', source)
        self.assertIn('content.hide()', source)

    def test_project_actions_are_embedded_above_design_step_1_1(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn(
            'self.new_project_button = QPushButton("New Project")', source)
        self.assertIn(
            'self.open_project_button = QPushButton("Open Project")', source)
        self.assertIn(
            'button.setObjectName("projectActionButton")', source)
        self.assertIn('button.setFixedHeight(34)', source)
        self.assertIn(
            'self.new_project_button.clicked.connect(self.new_project)',
            source)
        self.assertIn(
            'self.open_project_button.clicked.connect(self.open_project)',
            source)
        self.assertIn(
            'if not self._ensure_project_for_parameter_acceptance():',
            source)
        self.assertIn(
            '"Create DNA Moiré Project", "Create Project"', source)

    def test_accept_without_project_uses_non_resetting_project_prompt(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        module = ast.parse(source_path.read_text(encoding="utf-8"))
        window_class = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef)
            and node.name == "MoireDesignerWindow")
        ensure = next(
            node for node in window_class.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_ensure_project_for_parameter_acceptance")
        ensure_source = ast.unparse(ensure)
        self.assertIn("self._project_setup_selection", ensure_source)
        self.assertIn("self.recalculate()", ensure_source)
        self.assertIn("self._save_current_project", ensure_source)
        self.assertNotIn("self.apply_paper_preset", ensure_source)
        self.assertNotIn("self.seed_cross_section_picker.reset_default",
                         ensure_source)

    def test_history_buttons_show_success_feedback(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn('button.setText("✓")', source)
        self.assertIn('QTimer.singleShot(900, restore_button)', source)
        self.assertIn('action="undo"', source)
        self.assertIn('action="redo"', source)
        self.assertIn('撤销成功：已恢复到“%s”。', source)
        self.assertIn('重做成功：已恢复到“%s”。', source)
        self.assertIn('self._history_restore_token += 1', source)
        self.assertIn('self._finish_history_restore(token)', source)

    def test_seed_support_is_read_only_and_spacing_is_0_to_160(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn("self.seed_z1.setRange(128, 128)", source)
        self.assertIn("self.seed_z3.setRange(128, 128)", source)
        self.assertIn("self.seed_z1.hide()", source)
        self.assertIn("self.seed_z3.hide()", source)
        self.assertIn("self.seed_z1_overlap_readout", source)
        self.assertIn("self.seed_z3_overlap_readout", source)
        self.assertIn("return list(range(0, 161, 8))", source)
        for retired in ("s8x7", "8×7", "8x7", "r4x3", "seed87"):
            self.assertNotIn(retired, source)

    def test_kagome_kagome_stage2_and_stage3_use_kagome_generator(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn(
            '"square_square_c4", "kagome_kagome", "square_kagome")',
            source)
        self.assertIn(
            'if symmetry == "kagome_kagome":\n        return "kagome"',
            source)
        self.assertGreaterEqual(
            source.count("_sst_lattice_for_symmetry("), 3)
        self.assertGreaterEqual(
            source.count("layers_design_sequence_identical=bool("), 2)

    def test_preview_values_are_in_full_width_rows_above_both_canvases(self):
        preview_path = Path(__file__).parents[1] / "moire_designer" / \
            "preview.py"
        source = preview_path.read_text(encoding="utf-8")
        self.assertIn(
            "visual_twist = _right_handed_preview_rotation(self._angle)",
            source)
        self.assertIn(
            "_right_handed_preview_rotation(self._twist_angle)", source)
        self.assertNotIn("_draw_layer_legend", source)
        self.assertNotIn("include_spacing=False", source)
        self.assertIn("self._draw_twist_annotation(painter, scale)", source)
        self.assertIn('"Twist: %+.1f° (%s)"', source)
        self.assertIn('handedness = "right-handed"', source)
        self.assertIn('handedness = "left-handed"', source)
        self.assertIn("self._draw_period_annotation(painter, scale)", source)
        self.assertIn("self._draw_axial_dimensions(", source)
        self.assertIn("screen space, so they stay", source)
        self.assertIn('"1st layer\\n%s"', source)
        self.assertIn('"Spacing\\n%s"', source)
        self.assertIn('"Z1\\n%s"', source)
        self.assertIn("fully outside the Seed footprint", source)
        self.assertIn("first_step*first_vector.x()", source)
        self.assertIn("included_angle = math.radians(60.0)", source)
        self.assertIn("included_angle = math.radians(90.0)", source)
        module = ast.parse(source)
        classes = {
            node.name: node for node in module.body
            if isinstance(node, ast.ClassDef)
        }
        for name in ("BilayerPreview", "MoireTopViewPreview"):
            paint = next(
                node for node in classes[name].body
                if isinstance(node, ast.FunctionDef)
                and node.name == "paintEvent")
            self.assertNotIn("_draw_parameter_banner", ast.unparse(paint))
        mainwindow_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        mainwindow_source = mainwindow_path.read_text(encoding="utf-8")
        self.assertIn(
            '("Seed Z1", "#2a78d1")',
            mainwindow_source)
        self.assertIn(
            '("Seed Z2", "#d9dee3")',
            mainwindow_source)
        self.assertIn(
            '("Seed Z3", "#d65b74")',
            mainwindow_source)
        self.assertLess(
            mainwindow_source.index(
                "right_layout.addWidget(self.setup_preview_parameters)"),
            mainwindow_source.index(
                "right_layout.addWidget(self.setup_preview, 1)"))
        self.assertLess(
            mainwindow_source.index(
                "right_layout.addWidget(self.side_preview_parameters)"),
            mainwindow_source.index("right_layout.addWidget(self.preview, 1)"))
        self.assertIn("background: #070b10", mainwindow_source)
        self.assertNotIn("Current Twist Angle", mainwindow_source)
        self.assertNotIn("Current Moiré Period", mainwindow_source)
        self.assertIn("Twist angle: ", mainwindow_source)
        self.assertIn("Moiré period: ", mainwindow_source)
        self.assertIn(
            'width="100%%" cellspacing="0" cellpadding="0"',
            mainwindow_source)
        self.assertIn(
            "QLabel#previewCommonTitle { color: #f3f6f9; "
            "background: #070b10;",
            mainwindow_source)
        self.assertIn(
            "QLabel#previewParameterBanner { color: #f3f6f9; "
            "background: #070b10; border: 1px dashed #52616f; "
            "border-radius: 0; "
            "padding: 2px 6px; margin: 0; }",
            mainwindow_source)
        self.assertNotIn('bgcolor="#34424f"', mainwindow_source)
        self.assertIn(
            "combined_preview_layout.setSpacing(0)", mainwindow_source)
        self.assertIn(
            "combined_preview_layout.setContentsMargins(10, 10, 10, 8)",
            mainwindow_source)
        self.assertIn(
            '<table width="90%%" align="center" cellspacing="0" cellpadding="0"',
            mainwindow_source)
        self.assertIn(
            'style="color:#f3f6f9;font-size:11px;font-weight:650;"',
            mainwindow_source)
        self.assertIn('1st layer:&nbsp;%s', mainwindow_source)
        self.assertIn('Spacing:&nbsp;%s', mainwindow_source)
        self.assertIn('2nd layer:&nbsp;%s', mainwindow_source)
        self.assertIn(
            '<td width="20%%" align="left">SST sublattice</td>',
                      mainwindow_source)
        self.assertIn('<td width="20%%" align="left">Seed</td>',
                      mainwindow_source)
        self.assertIn('Z1:&nbsp;%s', mainwindow_source)
        self.assertIn('<tr><td height="11" colspan="4"></td></tr>',
                      mainwindow_source)
        self.assertIn(
            "font-size: 19px; font-weight: 750;", mainwindow_source)
        self.assertIn(
            "self.design_preview_title.setMinimumHeight(44)",
            mainwindow_source)
        self.assertIn(
            "self.design_preview_title.setContentsMargins(0, 4, 0, 0)",
            mainwindow_source)
        self.assertEqual(
            mainwindow_source.count("DNA Moiré Superlattice"), 2)
        self.assertIn("matrix.setContentsMargins(0, 4, 0, 0)",
                      mainwindow_source)
        self.assertIn("self.setup_preview_parameters.hide()",
                      mainwindow_source)
        self.assertIn("self.side_preview_parameters.hide()",
                      mainwindow_source)
        self.assertIn("self.design_preview_top_heading.hide()",
                      mainwindow_source)
        self.assertIn("self.design_preview_side_heading.hide()",
                      mainwindow_source)
        self.assertIn("self.design_preview_vertical_separator.hide()",
                      mainwindow_source)

        self.assertIn('setObjectName("previewSurface")', mainwindow_source)
        self.assertIn(
            "preview_surface_layout.setContentsMargins(0, 0, 0, 0)",
            mainwindow_source)
        self.assertIn(
            "preview_surface_layout.setSpacing(0)", mainwindow_source)
        self.assertNotIn("previewSectionDivider", mainwindow_source)
        self.assertIn(
            "QFrame#previewSurface, QFrame#previewMatrix",
            mainwindow_source)
        self.assertIn(
            '("SST sublattice 1st layer", "#2a78d1")',
            mainwindow_source)
        self.assertIn(
            '("SST sublattice 2nd layer", "#d65b74")',
            mainwindow_source)
        self.assertLess(
            mainwindow_source.index(
                '(\"SST sublattice 1st layer\", \"#2a78d1\")'),
            mainwindow_source.index('(\"Seed\", \"#ffffff\")'))
        self.assertLess(
            mainwindow_source.index('(\"Seed\", \"#ffffff\")'),
            mainwindow_source.index(
                '(\"SST sublattice 2nd layer\", \"#d65b74\")'))
        self.assertNotIn("SST superlattice helix", mainwindow_source)
        self.assertIn("PREVIEW_CONTENT_TOP = 12.0", source)
        self.assertIn("PREVIEW_CONTENT_BOTTOM = 42.0", source)
        self.assertIn("PREVIEW_MODEL_HEIGHT_FACTOR = 1.48", source)
        self.assertIn(
            "target_height = usable_height/PREVIEW_MODEL_HEIGHT_FACTOR",
            source)
        self.assertIn(
            "height*PREVIEW_MODEL_HEIGHT_FACTOR", source)
        self.assertIn("content_top+content_height*.50", source)
        self.assertEqual(mainwindow_source.count(
            'color:#ffffff;font-size:11px;font-weight:700;'), 2)
        self.assertNotIn(
            "matrix.addWidget(self.setup_preview_parameters, 0, 0)",
            mainwindow_source)
        self.assertIn(
            "matrix.addWidget(self.setup_preview_box, 0, 0)",
            mainwindow_source)
        self.assertEqual(source.count("_draw_model_coverage_frame("), 3)
        self.assertIn('QColor("#52616f")', source)
        self.assertNotIn("design_preview_channel", mainwindow_source)
        self.assertNotIn(
            "Left drag: rotate · Right drag: pan", source)
        self.assertNotIn(
            "左拖旋转 · 右拖平移", source)
        self.assertEqual(
            mainwindow_source.count(
                'setObjectName("previewGroupBox")'), 5)
        self.assertEqual(
            mainwindow_source.count("setSizes([460, 1000])"), 5)
        self.assertGreaterEqual(
            mainwindow_source.count("setMinimumWidth(455)"), 4)
        self.assertGreaterEqual(
            mainwindow_source.count("setMaximumWidth(535)"), 4)
        self.assertGreaterEqual(
            mainwindow_source.count("setMinimumWidth(360)"), 2)
        self.assertNotIn("right_panel.setMinimumWidth(640)",
                         mainwindow_source)
        self.assertNotIn("preview_panel.setMinimumWidth(640)",
                         mainwindow_source)
        self.assertEqual(
            mainwindow_source.count(
                "layout.setContentsMargins(0, 8, 0, 8)"), 3)
        self.assertIn("right.setContentsMargins(0, 0, 0, 0)",
                      mainwindow_source)
        self.assertIn(
            "root_layout.setContentsMargins(10, 8, 10, 0)",
            mainwindow_source)
        self.assertIn(
            "self.messageChanged.connect(self._sync_visibility)",
            mainwindow_source)
        self.assertIn(
            "self.setVisible(bool(str(message).strip()))",
            mainwindow_source)
        self.assertIn('"2 Automated DNA Design", "3  序列导出"',
                      mainwindow_source)
        self.assertIn('"Next: Automated DNA Design"', mainwindow_source)
        self.assertIn('"3.3 Final export")', mainwindow_source)
        self.assertNotIn('QGroupBox("4. Final Export")', mainwindow_source)
        self.assertIn(
            "QGroupBox#previewGroupBox::title {", mainwindow_source)
        self.assertIn("top: 2px;", mainwindow_source)
        self.assertIn("font-size: 17px;", mainwindow_source)
        self.assertIn("font-weight: 750;", mainwindow_source)

    def test_structure_and_sequence_preview_body_text_is_uniform(self):
        preview_path = Path(__file__).parents[1] / "moire_designer" / \
            "preview.py"
        preview_source = preview_path.read_text(encoding="utf-8")
        self.assertIn("def _preview_text_font", preview_source)
        self.assertIn("font.setPixelSize(11)", preview_source)
        self.assertIn(
            "painter.setFont(_preview_text_font(QFont.Weight.DemiBold))",
            preview_source)
        self.assertIn("Qt.TextFlag.TextWordWrap", preview_source)
        self.assertNotIn(
            "Scaffold routing Path · Drag/scroll to zoom", preview_source)
        self.assertNotIn(
            "Drag dividers to resize panels", preview_source)
        self.assertIn(
            'painter, rect, "Face %d" % (face_index + 1)',
            preview_source)
        self.assertNotIn('"base %d" % low', preview_source)
        self.assertNotIn('"base %d" % high', preview_source)

    def test_top_and_side_headings_have_no_flanking_dashed_lines(self):
        mainwindow_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = mainwindow_path.read_text(encoding="utf-8")
        self.assertIn("heading_layout.addWidget(subtitle, 1)", source)
        self.assertNotIn('setObjectName("previewDashedLine")', source)
        self.assertNotIn("QFrame#previewDashedLine", source)
        self.assertIn("(self.setup_preview_box, self.setup_preview)", source)
        self.assertIn("(self.side_preview_box, self.preview)", source)
        self.assertIn("preview_canvas.setMinimumHeight(260)", source)

    def test_compact_menu_and_no_calibration_restore_control(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn('menu_bar.addMenu("文件")', source)
        for removed_menu in ('编辑', '视图', '工具', '语言'):
            self.assertNotIn(
                'menu_bar.addMenu("%s")' % removed_menu, source)
        self.assertNotIn("language_actions", source)
        self.assertNotIn("_change_interface_language", source)
        self.assertNotIn('QPushButton("恢复校准基准参数")', source)
        self.assertNotIn("paper_preset_action_status", source)

    def test_live_preview_status_has_reviewed_english_translation(self):
        i18n_path = Path(__file__).parents[1] / "moire_designer" / "i18n.py"
        source = i18n_path.read_text(encoding="utf-8")
        self.assertIn(
            '"暂无设计参数": "No accepted Moiré parameters"', source)
        self.assertIn(
            '"尚未接受设计图": "No accepted DNA design"', source)
        self.assertIn(
            '"Live update: nominal Z2=%d bp, actual Z2/spacing=%.1f bp, "',
            source)

    def test_twist_label_omits_the_input_range(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn('twist_label = QLabel("Twist")', source)
        self.assertNotIn('QLabel("Twist（−45° 至 45°）")', source)
        self.assertIn(
            'twist_label.setStyleSheet("font-weight:700")', source)
        self.assertIn(
            'self.period_label.setStyleSheet("font-weight:700")', source)
        self.assertNotIn(
            'twist_label.setStyleSheet("color:', source)
        self.assertNotIn(
            'self.period_label.setStyleSheet("color:', source)
        self.assertIn(
            "self.angle.setRange(-1.0e9, 1.0e9)", source)
        self.assertNotIn(
            "self.angle.setRange(-45.0, 45.0)", source)
        self.assertIn("class SignedAngleSpinBox", source)
        self.assertIn('return "right-handed"', source)
        self.assertIn('return "left-handed"', source)

    def test_current_project_badge_is_not_in_the_workflow_bar(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertNotIn("self.current_project_label", source)
        self.assertNotIn('setObjectName("currentProject")', source)

    def test_capture_column_readouts_show_the_minimum(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertGreaterEqual(
            source.count("capture columns (minimum 4)"), 4)

    def test_indel_readout_states_the_structural_limit(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        source = source_path.read_text(encoding="utf-8")
        self.assertIn(
            'self.mean_indel.setSuffix(\n'
            '            " / helix (minimum -12, maximum +10)")',
            source)
        self.assertIn('"minimum_seed_deletion_per_helix"', source)
        self.assertIn('"seed_indel_limit_exceeded"', source)
        self.assertIn(
            'mean_indel_label = QLabel("Mean insertion/deletion")', source)

    def test_one_click_structure_generation_does_not_accept_intermediate(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / "mainwindow.py"
        module = ast.parse(source_path.read_text(encoding="utf-8"))
        function = next(
            node for node in ast.walk(module)
            if isinstance(node, ast.FunctionDef)
            and node.name == "generate_simple_structure_design")
        calls = [node for node in ast.walk(function) if isinstance(node, ast.Call)]
        method_names = {
            node.func.attr for node in calls if isinstance(node.func, ast.Attribute)}
        self.assertIn("generate_scaffold_design", method_names)
        self.assertIn("generate_complete_structure", method_names)
        self.assertNotIn("accept_scaffold", method_names)
        complete_call = next(
            node for node in calls
            if isinstance(node.func, ast.Attribute)
            and node.func.attr == "generate_complete_structure")
        self.assertIn("scaffold_source",
                      [keyword.arg for keyword in complete_call.keywords])


@unittest.skipIf(QApplication is None, "PyQt6 UI runtime is unavailable")
class MainWindowWorkflowUiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        with patch.object(MoireDesignerWindow, "_prompt_startup",
                          lambda unused_self: None):
            self.window = MoireDesignerWindow()
        self.window._localizer.retranslate()
        self.window.statusBar().retranslate()

    def tearDown(self):
        self.window.close()

    def test_design_and_project_independent_analysis_modes_round_trip(self):
        self.assertEqual(self.window._app_mode, "design")
        self.assertEqual(
            self.window.mode_switch_button.text(),
            "Switch to Analysis Mode")
        self.window._go_to_step(2)
        design_project = self.window.project
        self.window.project_path = "/tmp/design-state.moire.json"

        self.window._switch_app_mode()
        self.assertEqual(self.window._app_mode, "analysis")
        self.assertEqual(self.window.tabs.currentIndex(), 3)
        self.assertEqual(
            self.window.mode_switch_button.text(),
            "Switch to Design Mode")
        self.assertIs(self.window.project, design_project)
        self.assertEqual(
            self.window.project_path, "/tmp/design-state.moire.json")
        self.assertIsNone(self.window._project_output_dir("analysis/test"))
        with patch.object(self.window, "_save_current_project") as save:
            self.window._autosave_project()
        save.assert_not_called()

        self.window._open_analysis_module(1)
        self.window._switch_app_mode()
        self.assertEqual(self.window._app_mode, "design")
        self.assertEqual(self.window.tabs.currentIndex(), 1)
        self.assertIs(self.window.project, design_project)
        self.assertEqual(
            self.window.mode_switch_button.text(),
            "Switch to Analysis Mode")

        self.window._switch_app_mode()
        self.assertEqual(self.window.analysis_module_stack.currentIndex(), 1)

    def test_mixed_lattice_consistent_lengths_use_independent_inputs(self):
        self.window.project.settings.layers_design_sequence_identical = True
        self.window.project.settings.lattice_symmetry = "square_kagome"
        self.assertFalse(self.window._sst_layers_identical())

        self.window.project.settings.lattice_symmetry = "square_square"
        self.assertTrue(self.window._sst_layers_identical())

        self.window.project.settings.lattice_symmetry = "kagome_kagome"
        self.assertTrue(self.window._sst_layers_identical())

    def test_preview_parameters_and_simplified_structure_controls(self):
        self.assertEqual(
            self.window.workflow_buttons[0].text(),
            "1 Moiré Parameter Input")
        self.assertEqual(self.window.new_project_button.text(), "New Project")
        self.assertEqual(self.window.open_project_button.text(), "Open Project")
        self.assertEqual(self.window.new_project_button.height(), 34)
        self.assertEqual(self.window.open_project_button.height(), 34)
        self.assertLess(
            self.window.design_basis_left_layout.indexOf(
                self.window.project_action_bar),
            self.window.design_basis_left_layout.indexOf(
                self.window.design_parameter_sections))
        self.assertTrue(self.window.workflow_buttons[1].isHidden())
        for index in (0, 2, 3):
            self.assertEqual(self.window.workflow_buttons[index].width(), 190)
            self.assertEqual(self.window.workflow_buttons[index].height(), 32)
        self.assertTrue(self.window.history_back_button.isHidden())
        self.assertTrue(self.window.history_forward_button.isHidden())
        self.assertEqual(
            self.window.symmetry_step_button.text(),
            "1.1 Select Bilayer Symmetry")
        self.assertEqual(
            self.window.twist_period_step_button.text(),
            "1.2 Enter Twist or Moiré Period")
        self.assertEqual(
            self.window.sst_parameter_step_button.text(),
            "1.3 Enter SST Superlattice Parameters")
        self.assertEqual(
            self.window.accept_parameters_button.text(),
            "Accept Current Moiré Parameters")
        self.assertEqual(
            self.window.accept_parameters_button.objectName(),
            "primaryButton")
        self.assertTrue(self.window.confirm_design_basis_button.isHidden())
        self.assertTrue(self.window.design_basis_next_button.isHidden())
        self.assertEqual(
            self.window.sst_box.title(),
            "Layer Length and Interlayer Spacing")
        self.assertTrue(self.window.symmetry_box.isHidden())
        self.assertTrue(self.window.target_box.isHidden())
        self.assertTrue(self.window.sst_seed_parameter_container.isHidden())
        self.assertTrue(self.window.seed_cross_section_box.isHidden())
        self.window.symmetry_step_button.setChecked(True)
        self.window.twist_period_step_button.setChecked(True)
        self.app.processEvents()
        self.assertFalse(self.window.symmetry_box.isHidden())
        self.assertFalse(self.window.target_box.isHidden())
        self.assertTrue(self.window.seed_cross_section_box.isHidden())
        sections = self.window.design_parameter_section_widgets
        self.assertEqual(len(sections), 3)
        expected_pairs = (
            (self.window.symmetry_step_button, self.window.symmetry_box),
            (self.window.twist_period_step_button, self.window.target_box),
            (self.window.sst_parameter_step_button,
             self.window.sst_seed_parameter_container),
        )
        for section, expected in zip(sections, expected_pairs):
            section_layout = section.layout()
            self.assertEqual(section_layout.spacing(), 7)
            self.assertIs(section_layout.itemAt(0).widget(), expected[0])
            self.assertIs(section_layout.itemAt(1).widget(), expected[1])
        self.assertEqual(
            self.window.design_parameter_sections.layout().spacing(), 12)
        self.assertEqual(self.window.design_basis_left_layout.spacing(), 12)
        for button in (
                self.window.symmetry_step_button,
                self.window.twist_period_step_button,
                self.window.sst_parameter_step_button,
                self.window.accept_parameters_button,
                self.window.parameters_next_button):
            self.assertEqual(button.minimumHeight(), 34)
            self.assertEqual(button.maximumHeight(), 34)
        self.assertEqual(
            self.window.design_parameter_sections.sizePolicy().verticalPolicy(),
            QSizePolicy.Policy.Maximum)
        for section in sections:
            self.assertEqual(
                section.sizePolicy().verticalPolicy(),
                QSizePolicy.Policy.Maximum)
        navigation = self.window.design_parameter_navigation_layout
        self.assertIs(
            navigation.itemAt(0).widget(),
            self.window.accept_parameters_button)
        self.assertIs(
            navigation.itemAt(1).widget(),
            self.window.parameters_action_status)
        self.assertIs(
            navigation.itemAt(2).widget(),
            self.window.parameters_next_button)
        self.assertTrue(self.window.design_basis_action_status.isHidden())
        self.assertEqual(
            self.window.design_basis_left_layout.alignment(),
            Qt.AlignmentFlag.AlignTop)
        self.assertFalse(any(
            self.window.design_basis_left_layout.itemAt(index).layout() is
            self.window.design_basis_legacy_buttons_layout
            for index in range(
                self.window.design_basis_left_layout.count())))
        self.assertEqual(
            self.window.side_preview_box.title(),
            "")
        self.assertEqual(
            self.window.setup_preview_box.title(),
            "")
        self.assertEqual(
            self.window.design_combined_preview_box.title(),
            "2D Lattice and Moiré Preview")
        self.assertEqual(
            self.window.design_preview_top_heading.findChild(QLabel).text(),
            "Top view")
        self.assertEqual(
            self.window.design_preview_side_heading.findChild(QLabel).text(),
            "Side view")
        self.assertTrue(self.window.design_preview_top_heading.isHidden())
        self.assertTrue(self.window.design_preview_side_heading.isHidden())
        self.assertIsNone(
            self.window.design_preview_top_heading.findChild(
                QFrame, "previewDashedLine"))
        self.assertIsNone(
            self.window.design_preview_side_heading.findChild(
                QFrame, "previewDashedLine"))
        self.assertFalse(self.window.setup_preview_box.isHidden())
        self.assertFalse(self.window.side_preview_box.isHidden())
        self.assertEqual(
            self.window.capture_preview_box.title(),
            "Embedded cadnano design vision")
        self.assertTrue(self.window.structure_preview_channel.isHidden())
        self.assertEqual(
            self.window.sequence_preview_box.title(),
            "Sequence position and structure preview")
        self.assertTrue(self.window.sequence_preview_status.isHidden())
        self.assertIn("Twist angle:",
                      self.window.setup_preview_parameters.text())
        side_text = self.window.side_preview_parameters.text()
        self.assertIn("SST", side_text)
        self.assertIn("1st layer", side_text)
        self.assertIn("Seed", side_text)
        self.assertIn("Z1:", side_text)
        self.assertIn("color:#2a78d1", side_text)
        self.assertIn("color:#8a61bb", side_text)
        self.assertIn("color:#d65b74", side_text)
        self.assertEqual(
            self.window.side_preview_parameters.alignment(),
            Qt.AlignmentFlag.AlignCenter)
        self.assertLessEqual(self.window.angle.minimum(), -1.0e9)
        self.assertGreaterEqual(self.window.angle.maximum(), 1.0e9)
        self.assertTrue(self.window.setup_preview_parameters.isHidden())
        self.assertTrue(self.window.side_preview_parameters.isHidden())
        self.assertTrue(self.window.design_preview_summary.isHidden())
        self.assertAlmostEqual(self.window.preview._yaw, 0.0)
        self.assertAlmostEqual(self.window.preview._pitch, 0.0)
        self.assertTrue(self.window.seed_cross_section_preset.isHidden())
        self.assertFalse(
            self.window.seed_cross_section_preset_display.isHidden())
        self.assertEqual(
            self.window.seed_cross_section_preset_display.text(),
            "8×8 + 4×4 pore")
        self.assertTrue(self.window.structure_expert_button.isHidden())
        self.assertTrue(self.window.structure_next_button.isHidden())
        self.assertFalse(self.window.inspect_final_design_button.isHidden())
        self.assertFalse(self.window.inspect_final_design_button.isEnabled())
        self.assertIn("Optional",
                      self.window.inspect_final_design_button.text())

    def test_twist_above_45_is_preserved_for_later_feasibility_check(self):
        self.window.angle.setValue(90.0)
        self.app.processEvents()
        self.assertAlmostEqual(self.window.angle.value(), 90.0, places=6)
        self.assertEqual(
            self.window.angle.text(), "+90.0° (right-handed)")
        self.assertIn(
            "+90.0° (right-handed)",
            self.window.setup_preview_parameters.text())
        self.assertGreater(self.window.mean_indel.value(), 10.0)
        self.assertTrue(self.window.project.prediction[
            "seed_insertion_limit_exceeded"])
        self.assertIn("color:#c62828", self.window.mean_indel.styleSheet())

        self.window.angle.setValue(3.3)
        self.app.processEvents()
        self.assertLessEqual(self.window.mean_indel.value(), 10.0)
        self.assertNotIn("color:#c62828", self.window.mean_indel.styleSheet())

        self.window.angle.setValue(-20.0)
        self.app.processEvents()
        self.assertEqual(
            self.window.angle.text(), "-20.0° (left-handed)")
        self.assertIn(
            "-20.0° (left-handed)",
            self.window.setup_preview_parameters.text())

    def test_undo_then_redo_restores_parameter_without_losing_future(self):
        self.window._app_mode = "design"
        self.window._history = []
        self.window._history_index = -1
        initial_angle = self.window.angle.value()
        self.window._record_history("Initial state")

        self.window.angle.setValue(initial_angle + 1.0)
        self.app.processEvents()
        changed_angle = self.window.angle.value()
        self.assertGreaterEqual(len(self.window._history), 2)

        self.window._history_back()
        self.app.processEvents()
        self.assertAlmostEqual(self.window.angle.value(), initial_angle)
        self.assertTrue(self.window.history_forward_button.isEnabled())

        self.window._history_forward()
        self.app.processEvents()
        self.assertAlmostEqual(self.window.angle.value(), changed_angle)
        self.assertEqual(
            self.window._history_index, len(self.window._history) - 1)

    def test_seed_indel_and_actual_z2_follow_twist_without_clipping(self):
        self.window.sst_spacing.setCurrentIndex(
            self.window.sst_spacing.findData(32))
        self.window.angle.setValue(20.0)
        self.app.processEvents()

        expected_indel = self.window.project.settings.mean_indel_per_helix
        expected_z2 = self.window.project.prediction[
            "actual_z2_spacing_bp"]
        self.assertGreater(expected_indel, 12.0)
        self.assertAlmostEqual(
            self.window.mean_indel.value(), expected_indel, places=1)
        self.assertAlmostEqual(
            self.window.actual_z2_spacing.value(), expected_z2, places=1)
        self.assertIn(
            "Z2", self.window.side_preview_parameters.text())
        self.assertIn(
            "%.1f bp" % expected_z2,
            self.window.side_preview_parameters.text())
        self.assertGreaterEqual(
            self.window.side_preview_parameters.text().count(
                "%.1f bp" % expected_z2),
            2,
        )
        self.assertIn(
            "Seed Z1/Z2/Z3", self.window.design_preview_length_summary.text())
        self.assertIn(
            "/%s/" % ("%.1f" % expected_z2),
            self.window.design_preview_length_summary.text())
        self.assertFalse(hasattr(self.window, "seed_indel_limit_status"))

    def test_insertion_limit_is_reported_only_when_accepting_1b(self):
        self.window.angle.setValue(20.0)
        self.app.processEvents()
        workflow = self.window._workflow()
        workflow["design_basis_accepted"] = True
        with patch.object(QMessageBox, "warning") as warning:
            self.window.accept_parameters()
        warning.assert_called_once()
        self.assertIn("Reduce the Twist magnitude or increase the spacing",
                      warning.call_args.args[2])
        self.assertFalse(workflow.get("parameters_accepted", False))

    def test_deletion_limit_tracks_spacing_and_blocks_acceptance(self):
        self.window.sst_spacing.setCurrentIndex(
            self.window.sst_spacing.findData(32))
        self.window.angle.setValue(-10.0)
        self.app.processEvents()
        self.assertTrue(self.window.project.prediction[
            "seed_deletion_limit_exceeded"])
        self.assertIn("minimum -12", self.window.mean_indel.suffix())
        self.assertIn("color:#c62828", self.window.mean_indel.styleSheet())

        workflow = self.window._workflow()
        workflow["design_basis_accepted"] = True
        with patch.object(QMessageBox, "warning") as warning:
            self.window.accept_parameters()
        warning.assert_called_once()
        self.assertIn("Each 8-bp domain permits at most 3",
                      warning.call_args.args[2])
        self.assertFalse(workflow.get("parameters_accepted", False))

        self.window.sst_spacing.setCurrentIndex(
            self.window.sst_spacing.findData(64))
        self.window.angle.setValue(-20.0)
        self.app.processEvents()
        self.assertFalse(self.window.project.prediction[
            "seed_deletion_limit_exceeded"])
        self.assertIn("minimum -24", self.window.mean_indel.suffix())
        self.assertNotIn("color:#c62828", self.window.mean_indel.styleSheet())

    def test_sequence_expert_actions_are_ordered_on_right(self):
        layout = self.window.sst_template_actions.layout()
        buttons = [layout.itemAt(index).widget().text()
                   for index in range(layout.count())
                   if layout.itemAt(index).widget() is not None
                   and hasattr(layout.itemAt(index).widget(), "text")]
        self.assertEqual(
            buttons[-3:],
            ["Export Input Template", "正交序列设计",
             "Import Input Sequences"])
        self.assertIs(self.window.sst_template_actions.parentWidget(),
                      self.window.sequence_results_box)
        self.assertIs(self.window.accept_added_sst_button.parentWidget(),
                      self.window.sequence_results_box)

    def test_sequence_steps_expand_from_compact_buttons(self):
        buttons = (
            self.window.sequence_scaffold_step_button,
            self.window.sequence_sst_step_button,
            self.window.sequence_export_step_button,
        )
        contents = (
            self.window.sequence_scaffold_section_content,
            self.window.sequence_sst_section_content,
            self.window.sequence_export_section_content,
        )
        self.assertEqual(
            [button.text() for button in buttons],
            ["3.1 Assign scaffold sequences",
             "3.2 Assign SST sublattice input sequences",
             "3.3 Final export"])
        self.assertTrue(all(button.isCheckable() for button in buttons))
        self.assertTrue(all(content.isHidden() for content in contents))
        for button, content in zip(buttons, contents):
            button.click()
            self.app.processEvents()
            self.assertFalse(content.isHidden())
            button.click()
            self.app.processEvents()
            self.assertTrue(content.isHidden())

    def test_sequence_expert_mode_is_independently_available(self):
        self.assertTrue(self.window.sequence_expert_button.isEnabled())
        self.window.sequence_expert_button.setChecked(True)
        self.app.processEvents()
        self.assertFalse(self.window.sst_template_actions.isHidden())
        self.assertTrue(
            self.window.open_orthogonal_sequences_button.isEnabled())
        self.assertFalse(
            self.window.export_sst_input_template_button.isEnabled())
        self.assertFalse(
            self.window.import_sst_input_template_button.isEnabled())

    def test_accepted_parameter_summary_is_english_after_restore(self):
        workflow = self.window._workflow()
        workflow["design_basis_accepted"] = True
        workflow["parameters_accepted"] = True

        self.window._restore_structure_workflow()

        summary = self.window.accepted_parameters_summary.text()
        self.assertTrue(summary.startswith("Accepted parameters · Twist "))
        self.assertNotIn("已接受参数", summary)

    def test_one_click_structure_flow_keeps_only_final_acceptance(self):
        workflow = self.window._workflow()
        workflow["parameters_accepted"] = True

        def scaffold_stage():
            workflow["sst_two_layer"] = "/tmp/example_sst.json"
            workflow["scaffold_review"] = "/tmp/example_sst_scaffold.json"

        received = []

        def final_stage(scaffold_source=None):
            received.append(scaffold_source)
            workflow["structure_complete"] = "/tmp/example_final.json"

        with patch.object(self.window, "_structure_template_supported",
                          return_value=True), \
                patch.object(self.window, "generate_scaffold_design",
                             side_effect=scaffold_stage), \
                patch.object(self.window, "generate_complete_structure",
                             side_effect=final_stage), \
                patch.object(self.window, "_refresh_structure_preview"), \
                patch.object(self.window, "_record_history"):
            self.window.generate_simple_structure_design()

        self.assertEqual(received, ["/tmp/example_sst_scaffold.json"])
        self.assertNotIn("scaffold_accepted", workflow)
        self.assertFalse(workflow.get("structure_accepted", False))
        self.assertTrue(self.window.inspect_final_design_button.isEnabled())
        self.assertEqual(
            workflow["automatic_design_exports"], {
                "sst": "/tmp/example_sst.json",
                "sst_scaffold_routing": "/tmp/example_sst_scaffold.json",
                "sst_scaffold_routing_staple_capture":
                    "/tmp/example_final.json",
            })

    def test_optional_inspection_is_a_chooser_but_acceptance_uses_final(self):
        source_path = Path(__file__).parents[1] / "moire_designer" / \
            "mainwindow.py"
        module = ast.parse(source_path.read_text(encoding="utf-8"))
        functions = {
            node.name: node for node in ast.walk(module)
            if isinstance(node, ast.FunctionDef)
        }
        accept_source = ast.unparse(functions["accept_complete_structure"])
        connect_source = ast.unparse(functions["_connect_signals"])
        chooser_source = ast.unparse(
            functions["inspect_selected_design_in_cadnano"])
        newest_source = ast.unparse(
            functions["_latest_complete_structure_file"])
        self.assertIn("_latest_complete_structure_file", accept_source)
        self.assertIn("inspect_selected_design_in_cadnano", connect_source)
        self.assertIn("QFileDialog.getOpenFileName", chooser_source)
        self.assertIn("cadnano_inspection_files", chooser_source)
        self.assertIn("require_staples=True", newest_source)
        self.assertIn('workflow.get("structure_complete")', newest_source)
        self.assertIn("st_mtime_ns", newest_source)
        self.assertIn("capture_bridge_component_count", newest_source)

    def test_auto_input_slot_contains_unexpected_errors(self):
        def fail():
            raise KeyError("synthetic failure")

        self.window._auto_design_and_add_sst_inputs_impl = fail
        with patch.object(QMessageBox, "critical", return_value=None) as critical:
            self.window.auto_design_and_add_sst_inputs()
        critical.assert_called_once()
        self.assertFalse(self.window._auto_input_design_running)
        self.assertIn("synthetic failure",
                      self.window.sst_sequence_status.text())


if __name__ == "__main__":
    unittest.main()

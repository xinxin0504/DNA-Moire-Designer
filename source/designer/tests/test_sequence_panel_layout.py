import ast
import unittest
from pathlib import Path


SOURCE_PATH = Path(__file__).parents[1] / "moire_designer" / "mainwindow.py"


class SequencePanelLayoutSourceTests(unittest.TestCase):

    def test_automatic_multilength_inputs_generate_longest_first(self):
        source = (Path(__file__).parents[1] / "moire_designer" /
                  "orthogonal_sequence_tool.py").read_text(encoding="utf-8")
        self.assertIn(
            "key=lambda item: item[0], reverse=True", source)

    def test_vector_export_forces_the_square_cadnano_canvas(self):
        source = (Path(__file__).parents[1] / "moire_design_core" /
                  "vector_export_worker.py").read_text(encoding="utf-8")
        self.assertIn("LatticeType.Square", source)
        self.assertIn("forceLatticeType=True", source)

    @classmethod
    def setUpClass(cls):
        cls.module = ast.parse(SOURCE_PATH.read_text(encoding="utf-8"))

    @classmethod
    def function_source(cls, name):
        function = next(
            node for node in ast.walk(cls.module)
            if isinstance(node, ast.FunctionDef) and node.name == name)
        return ast.unparse(function)

    def test_sequence_results_use_a_vertical_splitter(self):
        source = self.function_source("_build_sequence_tab")
        self.assertIn(
            "self.sequence_results_splitter = QSplitter(Qt.Orientation.Vertical)",
            source)
        for widget in (
                "self.sequence_preview_box",
                "self.sst_sequence_table_box",
                "self.sst_template_actions"):
            self.assertIn(
                "self.sequence_results_splitter.addWidget(%s)" % widget,
                source)
        self.assertIn(
            "self.sequence_results_splitter.setChildrenCollapsible(True)",
            source)
        self.assertIn(
            "self.sequence_results_splitter.setOpaqueResize(True)", source)
        self.assertIn(
            "self.sequence_results_splitter.setCollapsible(index, True)",
            source)

    def test_result_panels_are_not_locked_by_minimum_heights(self):
        build_source = self.function_source("_build_sequence_tab")
        toggle_source = self.function_source("_toggle_sequence_expert")
        balance_source = self.function_source(
            "_balance_sequence_result_panels")
        for widget in (
                "self.sequence_preview_box",
                "self.sequence_preview",
                "self.sst_sequence_table",
                "self.sst_sequence_table_box",
                "self.sst_template_actions"):
            self.assertIn("%s.setMinimumHeight(0)" % widget, build_source)
        self.assertNotIn("setMinimumHeight(235)", build_source)
        self.assertNotIn("setMinimumHeight(205)", build_source)
        self.assertNotIn("setMinimumHeight(190)", build_source)
        self.assertIn("button.setMinimumHeight(40)", build_source)
        self.assertIn("QTimer.singleShot", toggle_source)
        self.assertIn("splitter.setSizes", balance_source)
        self.assertIn("splitter.sizes()", balance_source)
        self.assertIn("all((sizes[index] > 0", balance_source)
        self.assertIn("_sequence_result_visibility_state", balance_source)
        self.assertNotIn("max(splitter.height(), 720)", balance_source)
        self.assertIn(
            "QTimer.singleShot(0, self._balance_sequence_result_panels)",
            self.function_source("_populate_sst_sequence_table"))

    def test_sst_statuses_follow_their_actions_and_restore(self):
        source = self.function_source("_build_sequence_tab")
        ordered_tokens = [
            "sst_layout.addWidget(self.auto_design_sst_inputs_button)",
            "sst_layout.addWidget(self.sst_auto_import_status)",
            "sst_layout.addWidget(self.sequence_expert_button)",
            "sst_layout.addWidget(self.sst_expert_import_status)",
            "sst_layout.addWidget(self.accept_added_sst_button)",
            "sst_layout.addWidget(self.sst_acceptance_status)",
        ]
        positions = [source.index(token) for token in ordered_tokens]
        self.assertEqual(positions, sorted(positions))

        restore = self.function_source("_restore_structure_workflow")
        for token in (
                "sequence_sst_detected", "sequence_sst_import_method",
                "sequence_sst_import_status",
                "sequence_sst_acceptance_status",
                "analyze_sequence_design"):
            self.assertIn(token, restore)
        self.assertIn("automatic_orthogonal_input", restore)

        automatic = self.function_source(
            "_auto_design_and_add_sst_inputs_impl")
        expert = self.function_source("import_sequence_sst_template")
        accepted = self.function_source("accept_added_sst_inputs")
        self.assertIn('_record_sst_import_status(\'automatic\'', automatic)
        self.assertIn('_record_sst_import_status(\'expert\'', expert)
        self.assertIn("sst_acceptance_status", accepted)

    def test_expert_status_is_hidden_when_expert_mode_is_closed(self):
        toggle = self.function_source("_toggle_sequence_expert")
        self.assertIn("self.sst_expert_import_status.show()", toggle)
        self.assertIn("self.sst_expert_import_status.hide()", toggle)
        restore = self.function_source("_restore_structure_workflow")
        self.assertIn("self.sequence_expert_button.isChecked()", restore)
        self.assertIn("self.sst_expert_import_status.hide()", restore)

    def test_wrapped_feedback_reserves_its_complete_text_height(self):
        feedback = self.function_source("_show_action_feedback")
        self.assertIn("label.heightForWidth(width)", feedback)
        self.assertIn("label.setMinimumHeight(required)", feedback)
        self.assertIn("QTimer.singleShot(0, fit_height)", feedback)

        export = self.function_source("export_sequence_final_package")
        self.assertIn("Final export completed successfully.", export)
        self.assertIn("Output folder:", export)
        self.assertIn("Next: open the output folder", export)


if __name__ == "__main__":
    unittest.main()

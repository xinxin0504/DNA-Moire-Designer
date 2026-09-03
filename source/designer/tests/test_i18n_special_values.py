import ast
from pathlib import Path
import unittest

from moire_designer.i18n import set_language, translate


I18N_PATH = Path(__file__).parents[1] / "moire_designer" / "i18n.py"


class SpecialValueLocalizationTests(unittest.TestCase):
    def tearDown(self):
        set_language("en")

    def test_academic_english_special_value_labels(self):
        set_language("en")
        self.assertEqual(translate("自动建议"), "Auto (recommended)")
        self.assertEqual(translate("请输入"), "Enter a value")

    def test_particle_statistics_headers_use_academic_english(self):
        set_language("en")
        expected = {
            "Group归属": "Group ID",
            "颗粒归属": "Particle Classification",
            "面积 nm²": "Area (nm²)",
            "圆度": "Circularity",
            "边缘": "Edge Contact",
            "主轴长度 / 次轴宽度":
                "Major-axis Length / Minor-axis Width",
            "中心路径长度": "Centerline Path Length",
        }
        for source, target in expected.items():
            self.assertEqual(translate(source), target)

    def test_gel_workflow_uses_two_numbered_steps(self):
        set_language("en")
        self.assertEqual(translate("2. 分析"), "2. Analyze")
        self.assertEqual(translate("2. 开始分析"), "2. Start analysis")

    def test_crystal_analysis_button_keeps_step_three(self):
        set_language("en")
        self.assertEqual(
            translate("3. 自动识别并分析"),
            "3. Automatically identify and analyze",
        )
        self.assertEqual(
            translate("3. 批量自动识别并分析"),
            "3. Run batch detection and analysis",
        )

    def test_design_buttons_use_consistent_title_case(self):
        set_language("en")
        expected = {
            "恢复当前预设": "Restore Current Preset",
            "重新选择点阵 / Seed 截面":
                "Reselect Lattice / Seed Cross-section",
            "导入 Moiré 工程 (.moire.json)":
                "Import Moiré Project (.moire.json)",
            "在 cadnano 内专家编辑完成 Scaffold routing":
                "Edit Scaffold Routing in caDNAno",
            "载入 cadnano 专家编辑后的 JSON":
                "Load Expert-Edited caDNAno JSON",
            "接受 cadnano 当前保存的 Scaffold routing":
                "Accept Current Scaffold Routing Saved in caDNAno",
            "生成 Staple / Capture 设计":
                "Generate Staple/Capture Design",
            "在 cadnano 内专家编辑完成结构":
                "Edit Structure in caDNAno",
            "接受当前 Added Scaffold":
                "Accept assigned scaffold sequences",
        }
        for source, target in expected.items():
            self.assertEqual(translate(source), target)

    def test_ui_localizer_handles_spinbox_special_value_text(self):
        module = ast.parse(I18N_PATH.read_text(encoding="utf-8"))
        localizer = next(
            node for node in module.body
            if isinstance(node, ast.ClassDef) and node.name == "UiLocalizer")
        retranslate = next(
            node for node in localizer.body
            if isinstance(node, ast.FunctionDef) and
            node.name == "retranslate")
        source = ast.unparse(retranslate)
        self.assertIn("QSpinBox", source)
        self.assertIn("QDoubleSpinBox", source)
        self.assertIn("item.specialValueText", source)
        self.assertIn("item.setSpecialValueText", source)


if __name__ == "__main__":
    unittest.main()

import ast
import re
from pathlib import Path
import unittest


TOOL_PATH = (
    Path(__file__).parents[1]
    / "moire_designer"
    / "orthogonal_sequence_tool.py"
)


class OrthogonalSequenceReportEnglishTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = TOOL_PATH.read_text(encoding="utf-8")
        cls.module = ast.parse(cls.source)

    def test_tool_contains_no_chinese_user_facing_text(self):
        string_literals = [
            node.value
            for node in ast.walk(self.module)
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
        ]
        chinese_literals = [
            value for value in string_literals
            if re.search(r"[\u3400-\u9fff]", value)
        ]
        self.assertEqual(chinese_literals, [])

    def test_completion_report_uses_academic_english_labels(self):
        required = (
            "Orthogonal Sequence Design Complete",
            "Candidates evaluated:",
            "Background sequences screened:",
            "Candidates rejected by criterion:",
            "Same-orientation exact match",
            "Interstrand complementarity",
            "Primer3 Thermodynamic Analysis…",
        )
        for text in required:
            self.assertIn(text, self.source)


if __name__ == "__main__":
    unittest.main()

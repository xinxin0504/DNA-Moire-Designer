import re
import tempfile
import unittest
import zipfile
from pathlib import Path

from moire_designer.i18n import localize_xlsx


WORKBOOK_XML = b'''<?xml version="1.0" encoding="UTF-8"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheets><sheet name="\xe5\xba\x8f\xe5\x88\x97\xe5\x88\x86\xe6\x9e\x90" sheetId="1"/></sheets>
</workbook>'''

SHEET_XML = b'''<?xml version="1.0" encoding="UTF-8"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheetData><row r="1">
    <c r="A1" t="inlineStr"><is><t>\xe6\x9d\xa5\xe6\xba\x90</t></is></c>
    <c r="B1" t="inlineStr"><is><t>\xe4\xba\x92\xe8\xa1\xa5\xe5\xba\x8f\xe5\x88\x97\xef\xbc\x885\xe2\x80\xb2\xe2\x86\x923\xe2\x80\xb2\xef\xbc\x89</t></is></c>
  </row><row r="2">
    <c r="A2" t="inlineStr"><is><t>\xe6\x96\xb0\xe7\x94\x9f\xe6\x88\x90</t></is></c>
    <c r="B2" t="inlineStr"><is><t>\xe6\x96\xb0\xe5\xba\x8f\xe5\x88\x97-001</t></is></c>
    <c r="C2" t="inlineStr"><is><r><rPr><color rgb="FF000000"/></rPr><t>ACGT</t></r><r><rPr><color rgb="FFFF0000"/></rPr><t>TGCA</t></r></is></c>
    <c r="D2" t="inlineStr"><is><t>\xe9\x80\x9a\xe8\xbf\x87</t></is></c>
  </row></sheetData>
</worksheet>'''


class OrthogonalWorkbookLocalizationTests(unittest.TestCase):
    def test_without_openpyxl_preserves_rich_sequence_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            workbook_path = Path(directory) / "orthogonal-analysis.xlsx"
            styles = b"<styleSheet><marker>unchanged</marker></styleSheet>"
            with zipfile.ZipFile(workbook_path, "w") as workbook:
                workbook.writestr("xl/workbook.xml", WORKBOOK_XML)
                workbook.writestr("xl/worksheets/sheet1.xml", SHEET_XML)
                workbook.writestr("xl/styles.xml", styles)

            localize_xlsx(workbook_path, "en")

            with zipfile.ZipFile(workbook_path) as workbook:
                workbook_xml = workbook.read(
                    "xl/workbook.xml").decode("utf-8")
                sheet_xml = workbook.read(
                    "xl/worksheets/sheet1.xml").decode("utf-8")
                self.assertEqual(workbook.read("xl/styles.xml"), styles)

        combined_xml = workbook_xml + sheet_xml
        self.assertIsNone(re.search(r"[\u3400-\u9fff]", combined_xml))
        self.assertIn('name="Sequence Analysis"', workbook_xml)
        self.assertIn("Source", sheet_xml)
        self.assertIn("Reverse complement (5′→3′)", sheet_xml)
        self.assertIn("Newly generated", sheet_xml)
        self.assertIn("New sequence-001", sheet_xml)
        self.assertIn("Pass", sheet_xml)
        self.assertIn('<color rgb="FF000000"', sheet_xml)
        self.assertIn('<color rgb="FFFF0000"', sheet_xml)
        self.assertIn("ACGT", sheet_xml)
        self.assertIn("TGCA", sheet_xml)


if __name__ == "__main__":
    unittest.main()

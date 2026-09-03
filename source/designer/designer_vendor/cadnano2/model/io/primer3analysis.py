"""Primer3 thermodynamic analysis for short orthogonal DNA sequences."""

import math
import os
import posixpath
from itertools import combinations
from xml.etree import ElementTree
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, BadZipFile, ZipFile

from .orthogonalseq import sanitize_sequence
from .sequencexlsx import (_MAIN_NS, _OFFICE_REL_NS, _PACKAGE_REL_NS,
                           _cell_value, _column_index, _qname,
                           _shared_strings)


PRIMER3_SETTINGS = {
    "mv_conc": 50.0,
    "dv_conc": 10.0,
    "dntp_conc": 0.0,
    "dna_conc": 100.0,
    "temp_c": 37.0,
    "max_loop": 30,
}

ANALYSIS_LABELS = {
    "hairpin": "Hairpin",
    "homodimer": "Homodimer",
    "heterodimer": "Heterodimer",
}

RESULT_HEADERS = (
    "Analysis Type", "Sequence 1 Source", "Sequence 1 Name",
    "Sequence 1 (5′→3′)", "Sequence 2 Source", "Sequence 2 Name",
    "Sequence 2 (5′→3′)", "Tm (°C)", "ΔG (kcal/mol)",
    "ΔH (kcal/mol)", "ΔS (cal/(K·mol))", "Structure Status",
    "Predicted Pairing Structure", "Notes",
)

SOURCE_LABELS = {
    "输入": "Input",
    "Input": "Input",
    "新生成": "Newly generated",
    "Newly generated": "Newly generated",
    "导入": "Imported",
    "Imported": "Imported",
    "骨架链": "Scaffold",
    "Scaffold": "Scaffold",
}


class Primer3Unavailable(RuntimeError):
    """Raised when primer3-py is not installed or cannot be imported."""


class Primer3Cancelled(Exception):
    """Raised when the user cancels a batch analysis."""


def normalized_entries(entries):
    """Validate, de-duplicate and normalize imported sequence records."""
    normalized = []
    seen = set()
    for index, entry in enumerate(entries, 1):
        if isinstance(entry, dict):
            source = str(entry.get("source", "Imported") or
                         "Imported").strip()
            name = str(entry.get("name", "Sequence-%03d" % index) or
                       "Sequence-%03d" % index).strip()
            raw_sequence = entry.get("sequence", "")
        else:
            values = tuple(entry)
            if len(values) == 3:
                source, name, raw_sequence = values
            elif len(values) == 2:
                source, name, raw_sequence = "Imported", values[0], values[1]
            else:
                source, name, raw_sequence = \
                    "Imported", "Sequence-%03d" % index, values[0]
            source, name = str(source).strip(), str(name).strip()
        source = SOURCE_LABELS.get(source, source)
        if source == "Scaffold":
            continue
        sequence = sanitize_sequence(raw_sequence)
        if sequence is None:
            continue
        key = (source, name, sequence)
        if key in seen:
            continue
        seen.add(key)
        normalized.append({"source": source or "Imported",
                           "name": name or "Sequence-%03d" % index,
                           "sequence": sequence})
    return normalized


def read_primer3_text(filename):
    """Read ``sequence`` or ``name<TAB>sequence`` records from a TXT file."""
    entries = []
    errors = []
    with open(filename, "r", encoding="utf-8-sig") as input_file:
        for line_number, raw_line in enumerate(input_file, 1):
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            if "\t" in value:
                name, sequence = value.rsplit("\t", 1)
            elif "," in value:
                name, sequence = value.rsplit(",", 1)
            else:
                name, sequence = \
                    "Imported-%03d" % (len(entries) + 1), value
            sequence = sanitize_sequence(sequence)
            if sequence is None:
                errors.append(
                    "Line %d is not a valid A/C/G/T sequence." % line_number)
                continue
            entries.append(("Imported", name.strip(), sequence))
    return normalized_entries(entries), errors


def _sheet_path(workbook_file, requested_name):
    requested_names = ((requested_name,) if isinstance(requested_name, str)
                       else tuple(requested_name))
    workbook_root = ElementTree.fromstring(
        workbook_file.read("xl/workbook.xml"))
    relationship_id = None
    available = []
    for sheet in workbook_root.iter(_qname(_MAIN_NS, "sheet")):
        name = sheet.attrib.get("name", "").strip()
        available.append(name)
        if name in requested_names:
            relationship_id = sheet.attrib.get(_qname(_OFFICE_REL_NS, "id"))
            break
    if relationship_id is None:
        raise ValueError(
            "The XLSX workbook does not contain the required worksheet "
            "(%s). Available worksheets: %s" %
            (" or ".join(requested_names), ", ".join(available) or "none"))
    rels_root = ElementTree.fromstring(
        workbook_file.read("xl/_rels/workbook.xml.rels"))
    for relationship in rels_root.iter(
            _qname(_PACKAGE_REL_NS, "Relationship")):
        if relationship.attrib.get("Id") == relationship_id:
            target = relationship.attrib.get("Target", "")
            if target.startswith("/"):
                return target.lstrip("/")
            return posixpath.normpath(posixpath.join("xl", target))
    raise ValueError(
        "The worksheet data could not be located for: %s" %
        " or ".join(requested_names))


def read_primer3_workbook(filename):
    """Read input/generated strands from an orthogonal-sequence workbook."""
    try:
        with ZipFile(filename, "r") as workbook_file:
            shared_strings = _shared_strings(workbook_file)
            path = _sheet_path(
                workbook_file, ("Sequence Analysis", "序列分析"))
            sheet_root = ElementTree.fromstring(workbook_file.read(path))
    except (BadZipFile, KeyError, ElementTree.ParseError) as error:
        raise ValueError("The selected file is not a readable XLSX: %s" % error)

    table = []
    sheet_data = sheet_root.find(_qname(_MAIN_NS, "sheetData"))
    if sheet_data is None:
        raise ValueError(
            "The Sequence Analysis worksheet contains no tabular data.")
    for row in sheet_data.findall(_qname(_MAIN_NS, "row")):
        values = {}
        for cell in row.findall(_qname(_MAIN_NS, "c")):
            column = _column_index(cell.attrib.get("r"))
            if column is not None:
                values[column] = _cell_value(cell, shared_strings)
        if values:
            table.append(values)
    if not table:
        return []
    headers = dict((str(value).strip(), column)
                   for column, value in table[0].items())
    sequence_column = next(
        (headers[name] for name in
         ("Sequence (5′→3′)", "Sequence (5'→3')", "序列（5′→3′）")
         if name in headers), None)
    if sequence_column is None:
        raise ValueError(
            "The Sequence Analysis worksheet lacks a Sequence (5′→3′) "
            "column.")
    source_column = next(
        (headers[name] for name in ("Source", "来源") if name in headers),
        None)
    name_column = next(
        (headers[name] for name in ("Name", "名称") if name in headers),
        None)
    entries = []
    for index, row in enumerate(table[1:], 1):
        source = (str(row.get(source_column, "Imported")).strip()
                  if source_column is not None else "Imported")
        source = SOURCE_LABELS.get(source, source)
        if source == "Scaffold":
            continue
        # Orthogonal workbooks use these two classes.  Unknown non-empty
        # sources are retained as imported short strands rather than treated
        # as scaffold data.
        name = (str(row.get(name_column, "Imported-%03d" % index)).strip()
                if name_column is not None else "Imported-%03d" % index)
        entries.append((source or "Imported", name,
                        row.get(sequence_column, "")))
    return normalized_entries(entries)


def read_primer3_sequences(filename):
    extension = os.path.splitext(filename)[1].lower()
    if extension == ".xlsx":
        return read_primer3_workbook(filename), []
    return read_primer3_text(filename)


def create_primer3_analyzer(settings=None):
    """Create the reusable primer3-py thermodynamic analysis object."""
    try:
        from primer3.thermoanalysis import ThermoAnalysis
    except (ImportError, OSError) as error:
        raise Primer3Unavailable(
            "primer3-py could not be loaded. Install primer3-py in the "
            "Python environment used by cadnano: %s" % error)
    values = dict(PRIMER3_SETTINGS)
    if settings:
        values.update(settings)
    try:
        return ThermoAnalysis(**values)
    except Exception as error:
        raise Primer3Unavailable("Primer3 initialization failed: %s" % error)


def _structure_text(result):
    lines = getattr(result, "ascii_structure_lines", None)
    if lines:
        return "\n".join(str(line) for line in lines)
    value = getattr(result, "ascii_structure", "")
    return str(value or "")


def _structure_payload(line):
    return line.split("\t", 1)[1] if "\t" in line else line


def parse_primer3_structure(kind, structure):
    """Return base and pair data from Primer3's ASCII structure output."""
    lines = str(structure or "").splitlines()
    if kind == "hairpin":
        if len(lines) < 2:
            return {"type": "hairpin", "sequence": "", "pairs": ()}
        notation = _structure_payload(lines[0])
        sequence = "".join(base for base in _structure_payload(lines[1])
                           if base in "ACGT")
        notation = notation[:len(sequence)].ljust(len(sequence), "-")
        left = [index for index, marker in enumerate(notation)
                if marker == "/"]
        right = [index for index, marker in enumerate(notation)
                 if marker == "\\"]
        pair_count = min(len(left), len(right))
        pairs = tuple(zip(left[-pair_count:],
                          reversed(right[:pair_count]))) \
            if pair_count else ()
        return {"type": "hairpin", "sequence": sequence,
                "notation": notation, "pairs": pairs}

    if len(lines) < 4:
        return {"type": "dimer", "columns": (),
                "top_sequence": "", "bottom_sequence": ""}
    top_outer, top_inner, bottom_inner, bottom_outer = [
        _structure_payload(line) for line in lines[:4]]
    width = max(map(len, (top_outer, top_inner,
                          bottom_inner, bottom_outer)))
    values = [value.ljust(width) for value in
              (top_outer, top_inner, bottom_inner, bottom_outer)]
    columns = []
    top_sequence = []
    bottom_display = []
    top_index = 0
    bottom_display_index = 0
    for column in range(width):
        top = (values[1][column] if values[1][column] in "ACGT" else
               values[0][column] if values[0][column] in "ACGT" else "")
        bottom = (values[2][column] if values[2][column] in "ACGT" else
                  values[3][column] if values[3][column] in "ACGT" else "")
        paired = (values[1][column] in "ACGT" and
                  values[2][column] in "ACGT")
        record = {"column": column, "top": top, "bottom": bottom,
                  "paired": paired, "top_index": None,
                  "bottom_display_index": None}
        if top:
            record["top_index"] = top_index
            top_index += 1
            top_sequence.append(top)
        if bottom:
            record["bottom_display_index"] = bottom_display_index
            bottom_display_index += 1
            bottom_display.append(bottom)
        columns.append(record)
    return {"type": "dimer", "columns": tuple(columns),
            "top_sequence": "".join(top_sequence),
            # Primer3 displays strand 2 from 3′ to 5′, hence reverse it to
            # recover the originally supplied 5′→3′ sequence.
            "bottom_sequence": "".join(reversed(bottom_display))}


def _finite_number(value, scale=1.0, digits=3):
    try:
        value = float(value) / scale
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(value):
        return ""
    return round(value, digits)


def _result_row(kind, first, second, thermo_result, error=""):
    if thermo_result is None:
        return {
            "kind": kind, "kind_label": ANALYSIS_LABELS[kind],
            "first": first, "second": second,
            "tm": "", "dg": "", "dh": "", "ds": "",
            "structure_found": False, "structure": "", "error": error,
        }
    return {
        "kind": kind, "kind_label": ANALYSIS_LABELS[kind],
        "first": first, "second": second,
        "tm": _finite_number(getattr(thermo_result, "tm", None), digits=2),
        "dg": _finite_number(getattr(thermo_result, "dg", None),
                             scale=1000.0),
        "dh": _finite_number(getattr(thermo_result, "dh", None),
                             scale=1000.0),
        "ds": _finite_number(getattr(thermo_result, "ds", None)),
        "structure_found": bool(
            getattr(thermo_result, "structure_found", False)),
        "structure": _structure_text(thermo_result), "error": error,
    }


def run_primer3_analysis(entries, modes=("hairpin", "homodimer",
                                         "heterodimer"),
                         analyzer=None, progress=None, cancelled=None):
    """Analyze selected strands and return rows sorted by most negative ΔG."""
    entries = normalized_entries(entries)
    modes = tuple(mode for mode in modes if mode in ANALYSIS_LABELS)
    if not entries or not modes:
        return []
    analyzer = analyzer or create_primer3_analyzer()
    jobs = []
    if "hairpin" in modes:
        jobs.extend(("hairpin", entry, None) for entry in entries)
    if "homodimer" in modes:
        jobs.extend(("homodimer", entry, None) for entry in entries)
    if "heterodimer" in modes:
        jobs.extend(("heterodimer", first, second)
                    for first, second in combinations(entries, 2))

    rows = []
    total = len(jobs)
    for completed, (kind, first, second) in enumerate(jobs, 1):
        if cancelled and cancelled():
            raise Primer3Cancelled()
        first_sequence = first["sequence"]
        second_sequence = second["sequence"] if second else None
        error = ""
        thermo_result = None
        try:
            if kind in ("hairpin", "homodimer") and len(first_sequence) > 60:
                raise ValueError(
                    "Primer3 requires sequences of 60 nt or fewer for this "
                    "analysis.")
            if kind == "heterodimer" and len(first_sequence) > 60 and \
                    len(second_sequence) > 60:
                raise ValueError(
                    "Heterodimer analysis requires at least one sequence to "
                    "be 60 nt or fewer.")
            if kind == "hairpin":
                thermo_result = analyzer.calc_hairpin(
                    first_sequence, output_structure=True)
            elif kind == "homodimer":
                thermo_result = analyzer.calc_homodimer(
                    first_sequence, output_structure=True)
            else:
                thermo_result = analyzer.calc_heterodimer(
                    first_sequence, second_sequence, output_structure=True)
        except Exception as caught:
            error = str(caught)
        rows.append(_result_row(kind, first, second, thermo_result, error))
        if progress:
            progress(completed, total, "%s: %s%s" % (
                ANALYSIS_LABELS[kind], first["name"],
                (" + " + second["name"]) if second else ""))

    def sort_key(row):
        value = row["dg"]
        return (value == "", float(value) if value != "" else math.inf,
                row["kind_label"], row["first"]["name"],
                row["second"]["name"] if row["second"] else "")
    rows.sort(key=sort_key)
    return rows


def result_values(row):
    first = row["first"]
    second = row["second"] or {"source": "", "name": "", "sequence": ""}
    return (
        row["kind_label"], first["source"], first["name"],
        first["sequence"], second["source"], second["name"],
        second["sequence"], row["tm"], row["dg"], row["dh"], row["ds"],
        ("Structure found" if row["structure_found"] else
         "No structure found"),
        row["structure"], row["error"],
    )


def _column_name(index):
    name = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _cell(reference, value, header=False):
    style = ' s="1"' if header else ""
    if isinstance(value, (int, float)):
        return '<c r="%s"%s><v>%s</v></c>' % (reference, style, value)
    text = escape(str(value or ""))
    preserve = ' xml:space="preserve"' if "\n" in text or \
        text != text.strip() else ""
    return ('<c r="%s" t="inlineStr"%s><is><t%s>%s</t></is></c>' %
            (reference, style, preserve, text))


def _worksheet(headers, rows, widths):
    all_rows = [headers] + list(rows)
    xml_rows = []
    for row_number, row in enumerate(all_rows, 1):
        cells = [_cell("%s%d" % (_column_name(column), row_number), value,
                       row_number == 1)
                 for column, value in enumerate(row, 1)]
        xml_rows.append('<row r="%d">%s</row>' %
                        (row_number, "".join(cells)))
    columns = "".join(
        '<col min="%d" max="%d" width="%s" customWidth="1"/>' %
        (index, index, width) for index, width in enumerate(widths, 1))
    last_column = _column_name(len(headers))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/'
            'spreadsheetml/2006/main"><dimension ref="A1:%s%d"/>'
            '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" '
            'topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
            '</sheetView></sheetViews><cols>%s</cols><sheetData>%s</sheetData>'
            '<autoFilter ref="A1:%s%d"/></worksheet>' %
            (last_column, len(all_rows), columns, "".join(xml_rows),
             last_column, len(all_rows)))


def write_primer3_workbook(filename, rows):
    """Write sorted Primer3 results and fixed calculation settings to XLSX."""
    result_rows = [result_values(row) for row in rows]
    settings_rows = (
        ("Parameter", "Value"),
        ("Monovalent cation concentration", "50 mM"),
        ("Divalent cation concentration", "10 mM"),
        ("dNTP concentration", "0 mM"),
        ("DNA concentration", "100 nM"),
        ("ΔG calculation temperature", "37 °C"),
        ("Sorting rule", "Ascending ΔG (most negative/highest risk first)"),
        ("Calculation engine", "Primer3 / primer3-py"),
    )
    sheets = (
        ("Primer3 Analysis", RESULT_HEADERS, result_rows,
         (15, 13, 18, 48, 13, 18, 48, 12, 18, 18, 20, 14, 72, 42)),
        ("Settings", settings_rows[0], settings_rows[1:], (32, 62)),
    )
    overrides = "".join(
        '<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % index
        for index in range(1, len(sheets) + 1))
    content_types = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '%s<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        '</Types>') % overrides
    package_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        '</Relationships>')
    workbook_sheets = "".join(
        '<sheet name="%s" sheetId="%d" r:id="rId%d"/>' %
        (escape(sheet[0]), index, index)
        for index, sheet in enumerate(sheets, 1))
    workbook = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets>%s</sheets></workbook>') % workbook_sheets
    relationships = "".join(
        '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>' %
        (index, index) for index in range(1, len(sheets) + 1))
    workbook_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">%s'
        '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
        '</Relationships>') % (relationships, len(sheets) + 1)
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<fonts count="2"><font/><font><b/><color rgb="FFFFFFFF"/></font></fonts>'
        '<fills count="3"><fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FF385D8A"/></patternFill></fill></fills>'
        '<borders count="1"><border/></borders><cellStyleXfs count="1">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
        '</cellXfs></styleSheet>')
    with ZipFile(filename, "w", ZIP_DEFLATED) as workbook_file:
        workbook_file.writestr("[Content_Types].xml", content_types)
        workbook_file.writestr("_rels/.rels", package_rels)
        workbook_file.writestr("xl/workbook.xml", workbook)
        workbook_file.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        workbook_file.writestr("xl/styles.xml", styles)
        for index, unused_sheet in enumerate(sheets, 1):
            workbook_file.writestr(
                "xl/worksheets/sheet%d.xml" % index,
                _worksheet(sheets[index - 1][1], sheets[index - 1][2],
                           sheets[index - 1][3]))

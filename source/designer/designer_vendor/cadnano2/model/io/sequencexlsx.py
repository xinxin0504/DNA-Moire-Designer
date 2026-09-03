"""Dependency-free XLSX input/output for caDNAno sequence workflows."""

import posixpath
import re
from xml.etree import ElementTree
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, BadZipFile, ZipFile


HEADERS = ('Start', 'End', 'Sequence', 'Length', 'Color')
_MAIN_NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
_OFFICE_REL_NS = ('http://schemas.openxmlformats.org/officeDocument/'
                  '2006/relationships')
_PACKAGE_REL_NS = ('http://schemas.openxmlformats.org/package/2006/'
                   'relationships')


def _qname(namespace, name):
    return '{%s}%s' % (namespace, name)


def _column_name(index):
    name = ''
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _column_index(reference):
    match = re.match(r'([A-Za-z]+)', reference or '')
    if not match:
        return None
    index = 0
    for character in match.group(1).upper():
        index = index * 26 + ord(character) - 64
    return index - 1


def _cell_xml(reference, value, is_header=False, style_index=None):
    if is_header:
        style = ' s="1"'
    elif style_index is not None:
        style = ' s="%d"' % style_index
    else:
        style = ''
    if isinstance(value, int):
        return '<c r="%s"%s><v>%d</v></c>' % (reference, style, value)
    text = escape(str(value))
    preserve = ' xml:space="preserve"' if text != text.strip() else ''
    return ('<c r="%s" t="inlineStr"%s><is><t%s>%s</t></is></c>' %
            (reference, style, preserve, text))


def _normalized_color(value):
    color = str(value).strip().lower()
    if re.match(r'^#[0-9a-f]{6}$', color):
        return color
    return None


def _worksheet_xml(rows, color_styles=None):
    all_rows = [HEADERS] + list(rows)
    xml_rows = []
    for row_index, row in enumerate(all_rows, 1):
        row_style = None
        if row_index > 1 and color_styles and len(row) > 4:
            row_style = color_styles.get(_normalized_color(row[4]))
        cells = []
        for column_index, value in enumerate(row, 1):
            reference = '%s%d' % (_column_name(column_index), row_index)
            cells.append(_cell_xml(reference, value, row_index == 1,
                                   row_style))
        xml_rows.append('<row r="%d">%s</row>' %
                        (row_index, ''.join(cells)))
    last_row = len(all_rows)
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/'
            'spreadsheetml/2006/main">'
            '<dimension ref="A1:E%d"/>'
            '<sheetViews><sheetView workbookViewId="0">'
            '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" '
            'state="frozen"/></sheetView></sheetViews>'
            '<cols><col min="1" max="2" width="14" customWidth="1"/>'
            '<col min="3" max="3" width="60" customWidth="1"/>'
            '<col min="4" max="4" width="10" customWidth="1"/>'
            '<col min="5" max="5" width="12" customWidth="1"/></cols>'
            '<sheetData>%s</sheetData><autoFilter ref="A1:E%d"/>'
            '</worksheet>' % (last_row, ''.join(xml_rows), last_row))


def _write_workbook(filename, sheets, use_row_colors=False):
    sheets = [(name, list(rows)) for name, rows in sheets]
    row_colors = []
    if use_row_colors:
        row_colors = sorted(set(
            color
            for unused_name, rows in sheets
            for row in rows
            for color in [_normalized_color(row[4]) if len(row) > 4
                          else None]
            if color is not None))
    color_styles = dict((color, index + 2)
                        for index, color in enumerate(row_colors))
    sheet_overrides = ''.join(
        '<Override PartName="/xl/worksheets/sheet%d.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.'
        'spreadsheetml.worksheet+xml"/>' % index
        for index in range(1, len(sheets) + 1))
    content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
%s
<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>''' % sheet_overrides
    package_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>'''
    workbook_sheets = ''.join(
        '<sheet name="%s" sheetId="%d" r:id="rId%d"/>' %
        (escape(name), index, index)
        for index, (name, rows) in enumerate(sheets, 1))
    workbook = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets>%s</sheets>
</workbook>''' % workbook_sheets
    sheet_relationships = ''.join(
        '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/'
        'officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet%d.xml"/>' % (index, index)
        for index in range(1, len(sheets) + 1))
    workbook_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
%s
<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>''' % (sheet_relationships, len(sheets) + 1)
    color_fonts = ''.join(
        '<font><color rgb="FF%s"/></font>' % color[1:].upper()
        for color in row_colors)
    color_cell_xfs = ''.join(
        '<xf numFmtId="0" fontId="%d" fillId="0" borderId="0" '
        'xfId="0" applyFont="1"/>' % (index + 2)
        for index in range(len(row_colors)))
    styles = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="%d"><font/><font><b/><color rgb="FFFFFFFF"/></font>%s</fonts>
<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF0F766E"/><bgColor indexed="64"/></patternFill></fill></fills>
<borders count="1"><border/></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="%d"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>%s</cellXfs>
</styleSheet>''' % (len(row_colors) + 2, color_fonts,
                     len(row_colors) + 2, color_cell_xfs)

    with ZipFile(filename, 'w', ZIP_DEFLATED) as workbook_file:
        workbook_file.writestr('[Content_Types].xml', content_types)
        workbook_file.writestr('_rels/.rels', package_rels)
        workbook_file.writestr('xl/workbook.xml', workbook)
        workbook_file.writestr('xl/_rels/workbook.xml.rels', workbook_rels)
        workbook_file.writestr('xl/styles.xml', styles)
        for index, (name, rows) in enumerate(sheets, 1):
            workbook_file.writestr('xl/worksheets/sheet%d.xml' % index,
                                   _worksheet_xml(rows, color_styles))


def write_sequence_workbook(filename, input_rows, output_rows):
    """Write input and output sequence rows to a two-sheet XLSX file."""
    _write_workbook(filename, (('input', input_rows),
                               ('output', output_rows)),
                    use_row_colors=True)


def write_sequence_template(filename, scaffold_rows):
    """Write an editable one-sheet template for bulk sequence import."""
    _write_workbook(filename, (('input', scaffold_rows),))


def _shared_strings(workbook_file):
    if 'xl/sharedStrings.xml' not in workbook_file.namelist():
        return []
    root = ElementTree.fromstring(workbook_file.read('xl/sharedStrings.xml'))
    strings = []
    for item in root.findall(_qname(_MAIN_NS, 'si')):
        strings.append(''.join(
            element.text or ''
            for element in item.iter(_qname(_MAIN_NS, 't'))))
    return strings


def _cell_value(cell, shared_strings):
    cell_type = cell.attrib.get('t')
    if cell_type == 'inlineStr':
        return ''.join(
            element.text or ''
            for element in cell.iter(_qname(_MAIN_NS, 't')))
    value_node = cell.find(_qname(_MAIN_NS, 'v'))
    if value_node is None or value_node.text is None:
        return ''
    value = value_node.text
    if cell_type == 's':
        try:
            return shared_strings[int(value)]
        except (IndexError, ValueError):
            raise ValueError('Workbook contains an invalid shared string')
    if cell_type in ('str', 'e', 'd'):
        return value
    try:
        numeric = float(value)
        return int(numeric) if numeric.is_integer() else numeric
    except ValueError:
        return value


def _input_sheet_path(workbook_file):
    workbook_root = ElementTree.fromstring(
        workbook_file.read('xl/workbook.xml'))
    relationship_id = None
    for sheet in workbook_root.iter(_qname(_MAIN_NS, 'sheet')):
        if sheet.attrib.get('name', '').strip().lower() == 'input':
            relationship_id = sheet.attrib.get(
                _qname(_OFFICE_REL_NS, 'id'))
            break
    if relationship_id is None:
        raise ValueError("Workbook must contain a sheet named 'input'")

    rels_root = ElementTree.fromstring(
        workbook_file.read('xl/_rels/workbook.xml.rels'))
    for relationship in rels_root.iter(
            _qname(_PACKAGE_REL_NS, 'Relationship')):
        if relationship.attrib.get('Id') == relationship_id:
            target = relationship.attrib.get('Target', '')
            if target.startswith('/'):
                return target.lstrip('/')
            return posixpath.normpath(posixpath.join('xl', target))
    raise ValueError("Could not locate the 'input' worksheet data")


def read_sequence_template(filename):
    """Return ``(headers, rows)`` from the workbook's input sheet.

    Each data row is returned as ``(xlsx_row_number, five_cell_values)``.
    Both shared-string workbooks written by Excel and inline-string workbooks
    written by caDNAno are supported.
    """
    try:
        with ZipFile(filename, 'r') as workbook_file:
            shared_strings = _shared_strings(workbook_file)
            sheet_path = _input_sheet_path(workbook_file)
            sheet_root = ElementTree.fromstring(
                workbook_file.read(sheet_path))
    except (BadZipFile, KeyError, ElementTree.ParseError) as error:
        raise ValueError('The selected file is not a readable XLSX workbook: '
                         '%s' % error)

    parsed_rows = []
    sheet_data = sheet_root.find(_qname(_MAIN_NS, 'sheetData'))
    if sheet_data is None:
        raise ValueError("The 'input' sheet does not contain tabular data")
    for fallback_row_number, row in enumerate(
            sheet_data.findall(_qname(_MAIN_NS, 'row')), 1):
        try:
            row_number = int(row.attrib.get('r', fallback_row_number))
        except ValueError:
            row_number = fallback_row_number
        values = [''] * len(HEADERS)
        for cell in row.findall(_qname(_MAIN_NS, 'c')):
            column_index = _column_index(cell.attrib.get('r'))
            if column_index is not None and column_index < len(values):
                values[column_index] = _cell_value(cell, shared_strings)
        if any(value not in ('', None) for value in values):
            parsed_rows.append((row_number, values))
    if not parsed_rows:
        return (), []
    headers = tuple(str(value).strip() for value in parsed_rows[0][1])
    return headers, parsed_rows[1:]

"""Strict validation for XLSX-based scaffold sequence imports."""

import re

from .sequencexlsx import HEADERS


def _text(value):
    return '' if value is None else str(value).strip()


def _integer(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value) if value.is_integer() else None
    try:
        number = float(_text(value))
        return int(number) if number.is_integer() else None
    except ValueError:
        return None


def _import_targets(part):
    targets = {}
    for oligo in part.oligos():
        strand5p = oligo.strand5p()
        if (strand5p is None or oligo.isLoop() or
                not strand5p.strandSet().isScaffold()):
            continue
        start, end = oligo.sequenceEndpoints()
        targets[start] = {
            'oligo': oligo,
            'end': end,
            'length': oligo.actualLength(),
            'color': str(oligo.color()).lower()
        }
    return targets


def validate_sequence_import(part, headers, workbook_rows):
    """Validate a complete template before returning any apply operations.

    Returns ``(operations, errors)``. Each operation is ``(oligo, sequence,
    start, row_number)``. If any error exists, operations is always empty.
    """
    if tuple(headers) != HEADERS:
        return [], [
            "Header error: expected %s, found %s." %
            (' | '.join(HEADERS), ' | '.join(headers) if headers else
             '<no header row>')]

    targets = _import_targets(part)
    seen_starts = {}
    operations = []
    errors = []

    for row_number, values in workbook_rows:
        padded = list(values[:len(HEADERS)])
        padded.extend([''] * (len(HEADERS) - len(padded)))
        start, end, raw_sequence, raw_length, color = padded
        sequence = re.sub(r'\s+', '', _text(raw_sequence)).upper()
        # Partial imports are supported: rows without a sequence are ignored
        # and do not participate in design or duplicate checks.
        if not sequence:
            continue
        start = _text(start)
        row_label = "Row %d" % row_number
        chain_label = "chain %s" % (start or '<blank Start>')

        if not start:
            errors.append("%s: Start is blank." % row_label)
            continue
        if start in seen_starts:
            errors.append(
                "%s, %s, Start: duplicate chain; it already appears in "
                "row %d." % (row_label, chain_label, seen_starts[start]))
            continue
        seen_starts[start] = row_number

        target = targets.get(start)
        if target is None:
            errors.append(
                "%s, %s, Start: no matching scaffold chain exists in the "
                "current design." % (row_label, chain_label))
            continue

        row_has_error = False
        actual_end = _text(end)
        if actual_end != target['end']:
            errors.append(
                "%s, %s, End: expected %s, found %s." %
                (row_label, chain_label, target['end'],
                 actual_end or '<blank>'))
            row_has_error = True

        actual_length = _integer(raw_length)
        if actual_length != target['length']:
            errors.append(
                "%s, %s, Length: expected %d, found %s." %
                (row_label, chain_label, target['length'],
                 _text(raw_length) or '<blank>'))
            row_has_error = True

        actual_color = _text(color).lower()
        if actual_color != target['color']:
            errors.append(
                "%s, %s, Color: expected %s, found %s." %
                (row_label, chain_label, target['color'],
                 actual_color or '<blank>'))
            row_has_error = True

        invalid = [(index + 1, character)
                   for index, character in enumerate(sequence)
                   if character not in 'ACGT']
        if invalid:
            locations = ', '.join(
                "%r at base %d" % (character, index)
                for index, character in invalid[:20])
            if len(invalid) > 20:
                locations += ", and %d more" % (len(invalid) - 20)
            errors.append(
                "%s, %s, Sequence: invalid character(s): %s." %
                (row_label, chain_label, locations))
            row_has_error = True
        if len(sequence) != target['length']:
            errors.append(
                "%s, %s, Sequence: expected %d bases, found %d." %
                (row_label, chain_label, target['length'],
                 len(sequence)))
            row_has_error = True
        if not row_has_error:
            operations.append((target['oligo'], sequence,
                               start, row_number))

    return ([], errors) if errors else (operations, [])

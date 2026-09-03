"""Graphical, NUPACK-like rendering of Primer3 secondary structures."""

import math

from cadnano2.model.io.primer3analysis import parse_primer3_structure
import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), ['QLineF', 'QRectF', 'Qt'])
util.qtWrapImport('QtGui', globals(), [
    'QBrush', 'QColor', 'QFont', 'QPainter', 'QPen'])
util.qtWrapImport('QtWidgets', globals(), [
    'QGraphicsScene', 'QGraphicsView'])


BASE_COLORS = {
    'A': QColor('#2e9d58'),
    'T': QColor('#d65353'),
    'C': QColor('#347bc1'),
    'G': QColor('#e29732'),
}


class Primer3StructureView(QGraphicsView):
    """Draw paired bases, backbone, loops and bulges from Primer3 output."""

    def __init__(self, parent=None):
        super(Primer3StructureView, self).__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(
            QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setMinimumHeight(230)
        self._baseFont = QFont('Arial', 10)
        self._baseFont.setBold(True)
        self._labelFont = QFont('Arial', 9)
        self._titleFont = QFont('Arial', 11)
        self._titleFont.setBold(True)
        self._autoFit = True

    def wheelEvent(self, event):
        self._autoFit = False
        factor = 1.18 if event.angleDelta().y() > 0 else 1.0 / 1.18
        self.scale(factor, factor)

    def showEvent(self, event):
        super(Primer3StructureView, self).showEvent(event)
        if self._autoFit:
            self._fitStructure()

    def resizeEvent(self, event):
        super(Primer3StructureView, self).resizeEvent(event)
        if self._autoFit:
            self._fitStructure()

    def _fitStructure(self):
        bounds = self.scene().sceneRect()
        if not bounds.isEmpty() and self.viewport().width() > 40 and \
                self.viewport().height() > 40:
            self.fitInView(bounds, Qt.AspectRatioMode.KeepAspectRatio)

    def showResult(self, result):
        scene = self.scene()
        scene.clear()
        self.resetTransform()
        self._autoFit = True
        if not result or not result.get('structure'):
            message = result.get('error', '') if result else ''
            self._message(
                message or 'Primer3 did not return a drawable structure.')
            return
        parsed = parse_primer3_structure(
            result.get('kind', ''), result.get('structure', ''))
        if parsed['type'] == 'hairpin':
            self._drawHairpin(parsed, result)
        else:
            self._drawDimer(parsed, result)
        bounds = scene.itemsBoundingRect().adjusted(-35, -35, 35, 35)
        scene.setSceneRect(bounds)
        self._fitStructure()

    def _message(self, text):
        item = self.scene().addText(
            text or 'No drawable structure was found.', self._titleFont)
        item.setDefaultTextColor(QColor('#5b6470'))

    def _line(self, first, second, color, width=2.0,
              style=Qt.PenStyle.SolidLine, z=0):
        pen = QPen(color, width, style)
        item = self.scene().addLine(
            QLineF(first[0], first[1], second[0], second[1]), pen)
        item.setZValue(z)
        return item

    def _text(self, text, x, y, color='#374151', font=None, centered=False):
        item = self.scene().addText(str(text), font or self._labelFont)
        item.setDefaultTextColor(QColor(color))
        if centered:
            rectangle = item.boundingRect()
            item.setPos(x - rectangle.width() / 2.0,
                        y - rectangle.height() / 2.0)
        else:
            item.setPos(x, y)
        item.setZValue(4)
        return item

    def _base(self, position, base, index, strand):
        x, y = position
        radius = 13.0
        color = BASE_COLORS.get(base, QColor('#7b8794'))
        circle = self.scene().addEllipse(
            QRectF(x - radius, y - radius, radius * 2, radius * 2),
            QPen(QColor('#ffffff'), 1.2), QBrush(color))
        circle.setZValue(2)
        circle.setToolTip('%s: %s, position %d (from the 5′ end)' %
                          (strand, base, index + 1))
        letter = self.scene().addText(base, self._baseFont)
        letter.setDefaultTextColor(QColor('#ffffff'))
        rectangle = letter.boundingRect()
        letter.setPos(x - rectangle.width() / 2.0,
                      y - rectangle.height() / 2.0 - 1)
        letter.setZValue(3)

    def _legend(self, x, y):
        self._line((x, y), (x + 34, y), QColor('#6f9fc8'), 2.2,
                   Qt.PenStyle.DashLine)
        self._text('Base pairing', x + 42, y - 12)
        self._line((x + 135, y), (x + 169, y), QColor('#59636e'), 2.2)
        self._text('Backbone direction', x + 177, y - 12)
        offset = x + 285
        for base in 'ATCG':
            self._base((offset, y), base, 0, 'Legend')
            offset += 38

    @staticmethod
    def _fillArmBulges(coords, paired_indices, offset):
        bulges = []
        for first, second in zip(paired_indices, paired_indices[1:]):
            count = second - first - 1
            if count <= 0:
                continue
            bulges.append(tuple(range(first + 1, second)))
            first_point, second_point = coords[first], coords[second]
            for rank, index in enumerate(range(first + 1, second), 1):
                ratio = rank / float(count + 1)
                coords[index] = (
                    first_point[0] * (1.0 - ratio) + second_point[0] * ratio +
                    offset * math.sin(math.pi * ratio),
                    first_point[1] * (1.0 - ratio) + second_point[1] * ratio)
        return bulges

    def _drawHairpin(self, parsed, result):
        sequence = parsed['sequence']
        pairs = list(parsed['pairs'])
        title = '%s — %s' % (result['kind_label'], result['first']['name'])
        self._text(title, 15, 5, font=self._titleFont)
        if not sequence or not pairs:
            self._text(
                'No stable pairing was found; the sequence is shown linearly.',
                15, 35,
                       color='#6b7280')
            coords = dict((index, (50 + index * 34, 100))
                          for index in range(len(sequence)))
            for index in range(len(sequence) - 1):
                self._line(coords[index], coords[index + 1],
                           QColor('#59636e'), 2.2)
            for index, base in enumerate(sequence):
                self._base(coords[index], base, index, 'Single strand')
            return

        left = [pair[0] for pair in pairs]
        right_for_left = [pair[1] for pair in pairs]
        right = sorted(right_for_left)
        pair_count = len(pairs)
        coords = {}
        for rank, (left_index, right_index) in enumerate(pairs):
            y = 125 + (pair_count - 1 - rank) * 39
            coords[left_index] = (275, y)
            coords[right_index] = (505, y)

        left_bulges = self._fillArmBulges(coords, left, -58)
        right_bulges = self._fillArmBulges(coords, right, 58)

        inner_left = left[-1]
        inner_right = right[0]
        loop_indices = list(range(inner_left + 1, inner_right))
        for rank, index in enumerate(loop_indices, 1):
            ratio = rank / float(len(loop_indices) + 1)
            angle = math.pi * (1.0 - ratio)
            coords[index] = (390 + 115 * math.cos(angle),
                             125 - 95 * math.sin(angle))

        outer_y = 125 + (pair_count - 1) * 39
        for index in range(left[0] - 1, -1, -1):
            distance = left[0] - index
            coords[index] = (255 - min(distance, 3) * 8,
                             outer_y + distance * 34)
        for index in range(right[-1] + 1, len(sequence)):
            distance = index - right[-1]
            coords[index] = (525 + min(distance, 3) * 8,
                             outer_y + distance * 34)

        # Defensive fallback for unusual Primer3 notations with an unmarked
        # position outside the recognized stem/loop ranges.
        for index in range(len(sequence)):
            if index not in coords:
                previous = coords.get(index - 1, (390, outer_y + 40))
                coords[index] = (previous[0] + 34, previous[1])

        for index in range(len(sequence) - 1):
            self._line(coords[index], coords[index + 1],
                       QColor('#59636e'), 2.2, z=0)
        for left_index, right_index in pairs:
            self._line(coords[left_index], coords[right_index],
                       QColor('#6f9fc8'), 2.2,
                       Qt.PenStyle.DashLine, z=1)
        for index, base in enumerate(sequence):
            self._base(coords[index], base, index, 'Single strand')

        if loop_indices:
            self._text('Loop: %d nt' % len(loop_indices), 390, -8,
                       color='#785e9b', centered=True)
        for group in left_bulges + right_bulges:
            middle = coords[group[len(group) // 2]]
            direction = -1 if group in left_bulges else 1
            self._text('bulge %d nt' % len(group),
                       middle[0] + direction * 18, middle[1] - 12,
                       color='#a25b38', centered=(direction < 0))
        self._text('5′', coords[0][0] - 34, coords[0][1] - 13,
                   color='#1f4f7a', font=self._titleFont)
        self._text('3′', coords[len(sequence) - 1][0] + 18,
                   coords[len(sequence) - 1][1] - 13,
                   color='#1f4f7a', font=self._titleFont)
        self._legend(25, max(point[1] for point in coords.values()) + 55)

    def _drawDimer(self, parsed, result):
        columns = [column for column in parsed['columns']
                   if column['top'] or column['bottom']]
        title = '%s — %s + %s' % (
            result['kind_label'], result['first']['name'],
            result['second']['name'] if result['second'] else
            result['first']['name'])
        self._text(title, 15, 5, font=self._titleFont)
        if not columns:
            self._message(
                'Primer3 did not return a parseable dimer structure.')
            return
        minimum_column = min(column['column'] for column in columns)
        top_points = []
        bottom_points = []
        pair_columns = []
        base_spacing = 34
        for column in columns:
            x = 65 + (column['column'] - minimum_column) * base_spacing
            if column['top']:
                y = 145 if column['paired'] else 76
                top_points.append((column, (x, y)))
            if column['bottom']:
                y = 215 if column['paired'] else 284
                bottom_points.append((column, (x, y)))
            if column['paired'] and column['top'] and column['bottom']:
                pair_columns.append(column['column'])
                self._line((x, 145), (x, 215), QColor('#6f9fc8'), 2.2,
                           Qt.PenStyle.DashLine, z=1)

        for points in (top_points, bottom_points):
            for (unused_first, first), (unused_second, second) in zip(
                    points, points[1:]):
                self._line(first, second, QColor('#59636e'), 2.2)

        bottom_length = len(parsed['bottom_sequence'])
        for column, position in top_points:
            self._base(
                position, column['top'], column['top_index'], 'Sequence 1')
        for column, position in bottom_points:
            original_index = bottom_length - 1 - column['bottom_display_index']
            self._base(position, column['bottom'], original_index,
                       'Sequence 2')

        if top_points:
            self._text('5′', top_points[0][1][0] - 38,
                       top_points[0][1][1] - 13, '#1f4f7a', self._titleFont)
            self._text('3′', top_points[-1][1][0] + 18,
                       top_points[-1][1][1] - 13, '#1f4f7a', self._titleFont)
        if bottom_points:
            self._text('3′', bottom_points[0][1][0] - 38,
                       bottom_points[0][1][1] - 13, '#1f4f7a', self._titleFont)
            self._text('5′', bottom_points[-1][1][0] + 18,
                       bottom_points[-1][1][1] - 13, '#1f4f7a', self._titleFont)

        # Regions between two paired columns are internal loops or bulges.
        for first, second in zip(pair_columns, pair_columns[1:]):
            if second <= first + 1:
                continue
            between = [column for column in columns
                       if first < column['column'] < second]
            top_count = sum(bool(column['top']) for column in between)
            bottom_count = sum(bool(column['bottom']) for column in between)
            if not top_count and not bottom_count:
                continue
            x = 65 + (((first + second) / 2.0) - minimum_column) * base_spacing
            self._text('Loop/bulge: top %d nt, bottom %d nt' %
                       (top_count, bottom_count), x, 180,
                       color='#785e9b', centered=True)
        right_edge = max(point[0] for unused, point in top_points + bottom_points)
        self._legend(25, 350)
        self._text(
            'Outward-projecting bases are unpaired; dashed lines indicate '
            'predicted base pairs.', min(right_edge + 35, 620), 335,
            color='#6b7280')

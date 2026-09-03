"""Interactive Twist and Bend editor for ordinary single-lattice designs."""

from __future__ import division

import math
from copy import deepcopy

import cadnano2.util as util
from ..model.io.twistbend import (SNUPI_CALIBRATION_VERSION, TARGET_PITCH,
                                  TwistBendError,
                                  calibrate_saved_plan,
                                  estimate_global_twist, merge_plans,
                                  plan_add_twist, plan_bend,
                                  plan_remove_twist, validate_regions)
from ..model.io.simulationtwist import (SimulationAnalysisError,
                                        analyze_simulation_file)

util.qtWrapImport('QtCore', globals(), ['QPoint', 'QRectF', 'Qt', 'pyqtSignal'])
util.qtWrapImport('QtGui', globals(), ['QBrush', 'QColor', 'QFont', 'QPainter',
                                       'QPainterPath', 'QPen'])
util.qtWrapImport('QtWidgets', globals(), ['QAbstractItemView', 'QComboBox',
                                           'QCheckBox',
                                           'QDialog', 'QDialogButtonBox',
                                           'QDoubleSpinBox', 'QFormLayout',
                                           'QFileDialog',
                                           'QFrame', 'QGridLayout',
                                           'QGroupBox', 'QHBoxLayout', 'QLabel',
                                           'QListWidget', 'QMessageBox',
                                           'QPushButton', 'QSpinBox',
                                           'QScrollArea',
                                           'QSlider',
                                           'QStackedWidget', 'QTabWidget',
                                           'QVBoxLayout', 'QWidget'])


def _helix_data(part):
    """Describe real helix geometry, coverage, existing indels and unsafe bases."""
    result = {}
    for vh in part.getVirtualHelices():
        strands = [strand for strandSet in vh.getStrandSets()
                   for strand in strandSet]
        scaffold_intervals = [(strand.lowIdx(), strand.highIdx())
                              for strand in vh.scaffoldStrandSet()]
        staple_intervals = [(strand.lowIdx(), strand.highIdx())
                            for strand in vh.stapleStrandSet()]
        covered = set()
        forbidden = set()
        crossovers = set()
        deletion_protected = set()
        for strand in strands:
            low, high = strand.idxs()
            covered.update(range(low, high + 1))
            # Every endpoint is a nick, physical end, or crossover endpoint.
            forbidden.update((low, high))
            for connected, idx in ((strand.connection5p(), strand.idx5Prime()),
                                   (strand.connection3p(), strand.idx3Prime())):
                if (connected is not None and
                        connected.virtualHelix() is not vh):
                    crossovers.add((int(idx),
                                    int(connected.virtualHelix().number())))
        for strand in vh.stapleStrandSet():
            try:
                oligo = strand.oligo()
                length = (oligo.actualLength()
                          if hasattr(oligo, 'actualLength') else oligo.length())
            except (AttributeError, TypeError):
                continue
            if int(length) <= 21:
                deletion_protected.update(
                    range(strand.lowIdx(), strand.highIdx()+1))
        if strands:
            low = min(strand.lowIdx() for strand in strands)
            high = max(strand.highIdx() for strand in strands)
        else:
            low, high = 0, part.maxBaseIdx()
        forbidden.update(idx for idx in range(low, high + 1)
                         if idx not in covered)
        coord = vh.coord()
        x, y = part.latticeCoordToPositionXY(*coord)
        # cadnano's drawing radius is 1.125 nm (2.25 nm center spacing).
        # Twist/Bend mechanics use the user-confirmed 2.8 nm center spacing
        # while retaining the exact Square/Honeycomb staggered geometry.
        mechanical_scale = 2.8 / (2.0 * float(part.radius()))
        x, y = x * mechanical_scale, y * mechanical_scale
        result[vh.number()] = {
            'number': vh.number(),
            'coord': (float(x), float(y)), 'lattice_coord': coord,
            'low': low, 'high': high, 'forbidden': forbidden,
            'deletion_protected': deletion_protected,
            'scaffold_intervals': scaffold_intervals,
            'staple_intervals': staple_intervals,
            'crossovers': sorted(crossovers),
            'twist_per_base': float(part._twistPerBase),
            'twist_offset': float(part._twistOffset),
            'insertions': dict((int(idx), int(item.length()))
                               for idx, item in
                               part.insertions().get(coord, {}).items())}
    return result


def _preview_data_before_recorded_tasks(data, metadata):
    """Recover the pre-task indel baseline for exact metadata replay.

    Saved task transforms already represent the macroscopic deformation.
    Subtracting their edits from preview-only data prevents reopening a JSON
    from counting the same designed indels once as existing axial length and
    again as the restored Twist/Bend transform.  Planning still uses ``data``.
    """
    baseline = deepcopy(data)
    for plan in (metadata or {}).get('last_plans', ()):
        for edit in plan.get('edits', ()):
            number, idx = int(edit['helix']), int(edit['idx'])
            if number not in baseline:
                continue
            insertions = baseline[number].setdefault('insertions', {})
            value = int(insertions.get(idx, 0)) - int(edit['length'])
            if value:
                insertions[idx] = value
            else:
                insertions.pop(idx, None)
    return baseline


class HelixPicker(QWidget):
    selectionChanged = pyqtSignal()

    def __init__(self, part, data, parent=None):
        super(HelixPicker, self).__init__(parent)
        self._part, self._data = part, data
        self._selected = set()
        self._screen = {}
        self._lasso = None
        self._lassoOrigin = None
        self._lassoBase = set()
        self.setMinimumHeight(205)
        self.setToolTip('单击选择；在空白处拖框可批量选择；Command/Ctrl/Shift 可追加')

    def selected(self):
        return sorted(self._selected)

    def setSelected(self, values):
        self._selected = set(int(value) for value in values
                             if int(value) in self._data)
        self.selectionChanged.emit()
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#f6f8fb'))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        if not self._data:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             '当前设计没有 helix')
            return
        xs = [item['coord'][0] for item in self._data.values()]
        ys = [item['coord'][1] for item in self._data.values()]
        min_x, max_x, min_y, max_y = min(xs), max(xs), min(ys), max(ys)
        scale = min((self.width() - 36) / max(1.0, max_x - min_x + 2.3),
                    (self.height() - 36) / max(1.0, max_y - min_y + 2.3))
        self._screen = {}
        for number, item in self._data.items():
            x = 18 + (item['coord'][0] - min_x + 1.15) * scale
            y = 18 + (item['coord'][1] - min_y + 1.15) * scale
            radius = max(12.0, min(24.0, scale * 0.72))
            self._screen[number] = (x, y, radius)
            selected = number in self._selected
            painter.setPen(QPen(QColor('#1177cc') if selected else
                                QColor('#708090'), 2.4 if selected else 1.2))
            painter.setBrush(QBrush(QColor('#d7ecff') if selected else
                                    QColor('#ffffff')))
            painter.drawEllipse(QRectF(x - radius, y - radius,
                                       2 * radius, 2 * radius))
            painter.setPen(QColor('#12324a'))
            painter.drawText(QRectF(x - radius, y - radius,
                                    2 * radius, 2 * radius),
                             Qt.AlignmentFlag.AlignCenter, str(number))
        if self._lasso is not None:
            painter.setPen(QPen(QColor('#087dcc'), 1.5, Qt.PenStyle.DashLine))
            painter.setBrush(QColor(40, 145, 230, 42))
            painter.drawRect(self._lasso)

    def mousePressEvent(self, event):
        hit = None
        for number, (x, y, radius) in self._screen.items():
            if (event.position().x() - x) ** 2 + \
                    (event.position().y() - y) ** 2 <= radius ** 2:
                hit = number
                break
        modifiers = event.modifiers()
        additive = bool(modifiers & (Qt.KeyboardModifier.ControlModifier |
                                     Qt.KeyboardModifier.MetaModifier |
                                     Qt.KeyboardModifier.ShiftModifier))
        if hit is None:
            self._lassoOrigin = event.position()
            self._lasso = QRectF(self._lassoOrigin, self._lassoOrigin)
            self._lassoBase = set(self._selected) if additive else set()
            self.update()
            return
        if not additive:
            self._selected = {hit}
        elif hit in self._selected:
            self._selected.remove(hit)
        else:
            self._selected.add(hit)
        self.selectionChanged.emit()
        self.update()

    def mouseMoveEvent(self, event):
        if self._lassoOrigin is not None:
            self._lasso = QRectF(self._lassoOrigin,
                                 event.position()).normalized()
            self.update()

    def mouseReleaseEvent(self, event):
        if self._lassoOrigin is None:
            return
        selected = set(self._lassoBase)
        for number, (x, y, unused_radius) in self._screen.items():
            if self._lasso.contains(x, y):
                selected.add(number)
        self._selected = selected
        self._lasso = None
        self._lassoOrigin = None
        self._lassoBase = set()
        self.selectionChanged.emit()
        self.update()


class RangePicker(QWidget):
    rangeChanged = pyqtSignal(int, int)

    def __init__(self, maximum, data=None, parent=None):
        super(RangePicker, self).__init__(parent)
        self._maximum = max(1, int(maximum))
        self._data = data or {}
        self._selected = []
        self._start, self._end = 0, self._maximum
        self._dragHandle = None
        self.setMinimumHeight(150)
        self.setToolTip('显示选中 helix 的实际 Path；拖动左右把手设置起点和终点')

    def setSelected(self, selected):
        self._selected = [int(number) for number in selected
                          if int(number) in self._data]
        self.setMinimumHeight(max(150, 64 + len(self._selected) * 24))
        self.updateGeometry()
        self.update()

    def setRange(self, start, end):
        self._start, self._end = sorted((max(0, int(start)),
                                         min(self._maximum, int(end))))
        self.update()

    def _index(self, x):
        return max(0, min(self._maximum,
                          int(round((x - 58) / max(1, self.width() - 82) *
                                    self._maximum))))

    def _x(self, index):
        return 58 + int(index) / float(self._maximum) * max(1, self.width()-82)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#fbfcfd'))
        left, right = 58.0, self.width() - 24.0
        x0, x1 = self._x(self._start), self._x(self._end)
        painter.fillRect(QRectF(x0, 28, max(1, x1-x0), self.height()-48),
                         QColor(28, 134, 229, 35))
        painter.setPen(QColor('#334455'))
        painter.drawText(int(left), 18, '0')
        painter.drawText(self.width() - 54, 18, str(self._maximum))
        rows = self._selected or []
        if not rows:
            painter.drawText(QRectF(left, 35, right-left, 45),
                             Qt.AlignmentFlag.AlignCenter,
                             '请先选择 helix，以显示对应 Path')
        for rank, number in enumerate(rows):
            data = self._data[number]
            y = 44 + rank * 24
            painter.setPen(QColor('#324b60'))
            painter.drawText(6, y+5, 'H%d' % number)
            painter.setPen(QPen(QColor('#d7dfe5'), 1))
            painter.drawLine(int(left), y, int(right), y)
            for low, high in data.get('scaffold_intervals', []):
                painter.setPen(QPen(QColor('#2878c7'), 5))
                painter.drawLine(int(self._x(low)), y-4,
                                 int(self._x(high)), y-4)
            for low, high in data.get('staple_intervals', []):
                painter.setPen(QPen(QColor('#d64545'), 4))
                painter.drawLine(int(self._x(low)), y+5,
                                 int(self._x(high)), y+5)
        painter.setPen(QPen(QColor('#087dcc'), 2))
        painter.drawLine(int(x0), 26, int(x0), self.height()-18)
        painter.drawLine(int(x1), 26, int(x1), self.height()-18)
        painter.setBrush(QColor('#087dcc'))
        painter.drawRoundedRect(QRectF(x0-9, 21, 18, 16), 4, 4)
        painter.drawRoundedRect(QRectF(x1-9, 21, 18, 16), 4, 4)
        painter.setPen(QColor('#ffffff'))
        painter.drawText(QRectF(x0-9, 21, 18, 16),
                         Qt.AlignmentFlag.AlignCenter, 'L')
        painter.drawText(QRectF(x1-9, 21, 18, 16),
                         Qt.AlignmentFlag.AlignCenter, 'R')
        painter.setPen(QColor('#334455'))
        painter.drawText(int((x0 + x1) / 2 - 62), self.height() - 4,
                         '%d–%d（%d positions）' %
                         (self._start, self._end,
                          self._end - self._start + 1))

    def mousePressEvent(self, event):
        x = event.position().x()
        self._dragHandle = ('left' if abs(x-self._x(self._start)) <=
                            abs(x-self._x(self._end)) else 'right')
        self._moveHandle(x)

    def _moveHandle(self, x):
        value = self._index(x)
        if self._dragHandle == 'left':
            self._start = min(value, self._end-1)
        elif self._dragHandle == 'right':
            self._end = max(value, self._start+1)
        self.rangeChanged.emit(self._start, self._end)
        self.update()

    def mouseMoveEvent(self, event):
        if self._dragHandle is not None:
            self._moveHandle(event.position().x())

    def mouseReleaseEvent(self, event):
        if self._dragHandle is not None:
            self._moveHandle(event.position().x())
        self._dragHandle = None


class CalculationDiagram(QWidget):
    """Compact in-page explanation of the integer-indel calculation."""

    def __init__(self, kind, parent=None):
        super(CalculationDiagram, self).__init__(parent)
        self._kind = kind
        self.setMinimumHeight(112)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#f5f8fb'))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        w, h = self.width(), self.height()
        if self._kind == 'twist':
            centers = ((55, h/2), (w-65, h/2))
            for rank, (cx, cy) in enumerate(centers):
                painter.setPen(QPen(QColor('#6b7d8d'), 2))
                painter.setBrush(QBrush(QColor('#ffffff')))
                painter.drawEllipse(QRectF(cx-29, cy-29, 58, 58))
                angle = 0 if rank == 0 else math.radians(55)
                painter.setPen(QPen(QColor('#087dcc'), 4))
                painter.drawLine(int(cx), int(cy),
                                 int(cx+25*math.cos(angle)),
                                 int(cy-25*math.sin(angle)))
                painter.drawText(int(cx-26), int(cy+45),
                                 'Start 0°' if rank == 0 else 'End Δθ')
            painter.setPen(QPen(QColor('#00a88f'), 3))
            painter.drawLine(92, int(h/2), w-105, int(h/2))
            painter.drawText(105, 24,
                '均匀分配整数 indel → 全部选区的平均旋转接近目标')
            painter.setPen(QPen(QColor('#f3a323'), 3))
            for x in range(120, max(121, w-120), max(24, int((w-240)/6) or 24)):
                painter.drawLine(x, int(h/2-7), x, int(h/2+7))
        else:
            painter.setPen(QPen(QColor('#738596'), 4))
            for offset in (-22, 0, 22):
                previous = None
                for step in range(31):
                    t = step / 30.0
                    angle = math.radians(55*t)
                    radius = 135 + offset
                    point = (55 + radius*(1-math.cos(angle)),
                             h-15-radius*math.sin(angle))
                    if previous:
                        painter.drawLine(int(previous[0]), int(previous[1]),
                                         int(point[0]), int(point[1]))
                    previous = point
            painter.setPen(QPen(QColor('#3498db'), 3))
            painter.drawText(w-175, 28, '外侧 + insertion')
            painter.setPen(QPen(QColor('#e5a000'), 3))
            painter.drawText(34, h-12, '内侧 − deletion')
            painter.setPen(QColor('#41566a'))
            painter.drawText(170, h-12,
                'ΔN ≈ θ·d / 0.34；2.8 nm 点阵坐标确定中性轴')


class AngleDial(QWidget):
    """Direct-manipulation angle selector used by Add Twist and Bend."""

    valueChanged = pyqtSignal(float)

    def __init__(self, label, parent=None):
        super(AngleDial, self).__init__(parent)
        self._label = str(label)
        self._value = 0.0
        self._dragging = False
        self.setMinimumSize(150, 150)
        self.setMaximumHeight(190)
        self.setToolTip('拖动圆周上的箭头选择角度；数值框仍可精确输入')

    def setValue(self, value):
        self._value = float(value) % 360.0
        self.update()

    def value(self):
        return self._value

    def _setFromPosition(self, position):
        cx, cy = self.width() / 2.0, self.height() / 2.0 + 8
        angle = math.degrees(math.atan2(cy - position.y(),
                                       position.x() - cx)) % 360.0
        if abs(angle - self._value) > 1e-6:
            self._value = angle
            self.valueChanged.emit(angle)
            self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#f5f8fb'))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        cx, cy = self.width() / 2.0, self.height() / 2.0 + 8
        radius = max(34.0, min(self.width(), self.height()) * .31)
        painter.setPen(QPen(QColor('#8395a7'), 2))
        painter.setBrush(QColor('#ffffff'))
        painter.drawEllipse(QRectF(cx-radius, cy-radius,
                                   radius*2, radius*2))
        painter.setPen(QPen(QColor('#d8e1e8'), 1))
        for angle in range(0, 360, 30):
            rad = math.radians(angle)
            x0, y0 = cx + (radius-6)*math.cos(rad), cy - (radius-6)*math.sin(rad)
            x1, y1 = cx + radius*math.cos(rad), cy - radius*math.sin(rad)
            painter.drawLine(int(x0), int(y0), int(x1), int(y1))
        rad = math.radians(self._value)
        tip = (cx + radius*math.cos(rad), cy - radius*math.sin(rad))
        painter.setPen(QPen(QColor('#087dcc'), 5))
        painter.drawLine(int(cx), int(cy), int(tip[0]), int(tip[1]))
        painter.setBrush(QColor('#087dcc'))
        painter.drawEllipse(QRectF(tip[0]-7, tip[1]-7, 14, 14))
        painter.setPen(QColor('#243b53'))
        painter.drawText(QRectF(0, 4, self.width(), 24),
                         Qt.AlignmentFlag.AlignCenter, self._label)
        painter.drawText(QRectF(cx-42, cy-13, 84, 26),
                         Qt.AlignmentFlag.AlignCenter, '%.1f°' % self._value)

    def mousePressEvent(self, event):
        self._dragging = True
        self._setFromPosition(event.position())

    def mouseMoveEvent(self, event):
        if self._dragging:
            self._setFromPosition(event.position())

    def mouseReleaseEvent(self, event):
        if self._dragging:
            self._setFromPosition(event.position())
        self._dragging = False


class BendDirectionPicker(QWidget):
    """Actual helix cross-section with a draggable bend-direction arrow."""

    valueChanged = pyqtSignal(float)

    def __init__(self, data, parent=None):
        super(BendDirectionPicker, self).__init__(parent)
        self._data = data
        self._selected = set()
        self._value = 0.0
        self._center = (self.width() / 2.0, self.height() / 2.0)
        self._dragging = False
        self.setMinimumHeight(260)
        self.setToolTip('截面与左侧 helix 位置一致；拖动箭头设置 Bend 方向')

    def setSelected(self, selected):
        self._selected = set(int(number) for number in selected)
        self.update()

    def setValue(self, value):
        self._value = float(value) % 360.0
        self.update()

    def _setFromPosition(self, position):
        cx, cy = self._center
        angle = math.degrees(math.atan2(cy - position.y(),
                                       position.x() - cx)) % 360.0
        self._value = angle
        self.valueChanged.emit(angle)
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#f5f8fb'))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setPen(QColor('#304b62'))
        painter.drawText(12, 22,
                         'Helix 截面（蓝色箭头指向弯曲内侧）')
        if not self._data:
            return
        xs = [item['coord'][0] for item in self._data.values()]
        ys = [item['coord'][1] for item in self._data.values()]
        min_x, max_x, min_y, max_y = min(xs), max(xs), min(ys), max(ys)
        available_w, available_h = self.width() - 50, self.height() - 70
        scale = min(available_w / max(2.8, max_x-min_x+2.8),
                    available_h / max(2.8, max_y-min_y+2.8))
        screen = {}
        for number, item in self._data.items():
            x = 25 + (item['coord'][0]-min_x+1.4)*scale
            y = 42 + (item['coord'][1]-min_y+1.4)*scale
            screen[number] = (x, y)
            selected = number in self._selected
            radius = max(7.0, min(17.0, scale*.38))
            painter.setPen(QPen(QColor('#087dcc') if selected else
                                QColor('#8294a5'), 2 if selected else 1))
            painter.setBrush(QColor('#bde8ff') if selected else
                             QColor('#ffffff'))
            painter.drawEllipse(QRectF(x-radius, y-radius,
                                       radius*2, radius*2))
            painter.setPen(QColor('#28445b'))
            painter.drawText(QRectF(x-radius, y-radius, radius*2, radius*2),
                             Qt.AlignmentFlag.AlignCenter, str(number))
        selected_points = [screen[number] for number in self._selected
                           if number in screen]
        if not selected_points:
            selected_points = list(screen.values())
        cx = sum(point[0] for point in selected_points) / len(selected_points)
        cy = sum(point[1] for point in selected_points) / len(selected_points)
        self._center = (cx, cy)
        length = max(42.0, min(78.0, min(self.width(), self.height())*.25))
        rad = math.radians(self._value)
        tip = (cx + length*math.cos(rad), cy - length*math.sin(rad))
        painter.setPen(QPen(QColor('#006dcc'), 5))
        painter.drawLine(int(cx), int(cy), int(tip[0]), int(tip[1]))
        painter.setBrush(QColor('#006dcc'))
        painter.drawEllipse(QRectF(tip[0]-7, tip[1]-7, 14, 14))
        painter.setPen(QColor('#173f5f'))
        painter.drawText(int(cx+8), int(cy-8), '%.1f°' % self._value)

    def mousePressEvent(self, event):
        self._dragging = True
        self._setFromPosition(event.position())

    def mouseMoveEvent(self, event):
        if self._dragging:
            self._setFromPosition(event.position())

    def mouseReleaseEvent(self, event):
        if self._dragging:
            self._setFromPosition(event.position())
        self._dragging = False


class DeformationPreview(QWidget):
    """Painter-based whole-design 3D preview with a live selected overlay."""

    def __init__(self, data, parent=None):
        super(DeformationPreview, self).__init__(parent)
        self._data = data
        self._selected, self._start, self._end = [], 0, 1
        self._mode, self._angle, self._direction = 'remove_twist', 0.0, 0.0
        self._pipeline = []
        self._selectionCenter = (0.0, 0.0)
        self._edits = []
        self._geometryCache = None
        self._simulationPaths = None
        self._centerCache, self._bendCache = {}, {}
        self._yaw, self._pitch, self._zoom = -0.62, -0.38, 1.0
        self._last = QPoint()
        self._drag = False
        self.setMinimumSize(480, 400)
        self.setToolTip('拖动旋转 · 滚轮缩放；灰色为原始整体，彩色为实时目标区域')

    def setState(self, selected, start, end, mode, angle, direction):
        self._selected = list(selected)
        self._start, self._end = int(start), int(end)
        self._mode, self._angle = str(mode), float(angle)
        self._direction = float(direction)
        coords = [self._data[number]['coord'] for number in self._selected
                  if number in self._data]
        if coords:
            self._selectionCenter = (
                sum(point[0] for point in coords) / len(coords),
                sum(point[1] for point in coords) / len(coords))
        self._geometryCache = None
        self.update()

    def setEdits(self, edits):
        self._edits = [dict(edit) for edit in edits]
        self.update()

    def setPipeline(self, transforms):
        """Display an ordered, not-yet-applied virtual deformation state."""
        updated = [dict(item) for item in (transforms or ())]
        if updated == self._pipeline:
            return
        self._pipeline = updated
        self._geometryCache = None
        self.update()

    def setSimulationGeometry(self, centerlines):
        """Use measured helix centre lines as the preview's base geometry."""
        if not centerlines:
            self._simulationPaths = None
        else:
            self._simulationPaths = dict(
                (int(number), tuple(((float(row[1]), float(row[2]),
                                      float(row[3])), int(row[0]))
                                    for row in rows))
                for number, rows in centerlines)
        self._geometryCache = None
        self.update()

    def _zForIndex(self, data, idx):
        """Actual input-axis coordinate including pre-existing indels.

        Bend-plan indels are the mechanism that realizes the different inner
        and outer arc lengths already represented by ``_deform``. Applying
        them here as well would count that length difference twice.
        """
        delta = sum(int(length) for position, length in
                    data.get('insertions', {}).items()
                    if int(position) <= int(idx))
        return (int(idx) + delta) * .34

    @staticmethod
    def _vadd(a, b):
        return tuple(a[i]+b[i] for i in range(3))

    @staticmethod
    def _vsub(a, b):
        return tuple(a[i]-b[i] for i in range(3))

    @staticmethod
    def _vscale(a, value):
        return tuple(component*value for component in a)

    @staticmethod
    def _dot(a, b):
        return sum(a[i]*b[i] for i in range(3))

    @staticmethod
    def _cross(a, b):
        return (a[1]*b[2]-a[2]*b[1],
                a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0])

    @classmethod
    def _unit(cls, value, fallback=(0.0, 0.0, 1.0)):
        length = math.sqrt(max(0.0, cls._dot(value, value)))
        return (cls._vscale(value, 1.0/length)
                if length > 1e-10 else tuple(fallback))

    @classmethod
    def _rotateAxis(cls, vector, axis, angle):
        """Rodrigues rotation used for Bishop-frame material rotation."""
        axis = cls._unit(axis)
        ca, sa = math.cos(angle), math.sin(angle)
        return cls._vadd(
            cls._vadd(cls._vscale(vector, ca),
                      cls._vscale(cls._cross(axis, vector), sa)),
            cls._vscale(axis, cls._dot(axis, vector)*(1.0-ca)))

    def _taskCenterRaw(self, task, idx):
        coords = [self._data[number]['coord'] for number in task['helices']
                  if number in self._data]
        if not coords:
            return (0.0, 0.0, float(idx)*.34)
        return (sum(point[0] for point in coords)/float(len(coords)),
                sum(point[1] for point in coords)/float(len(coords)),
                float(idx)*.34)

    def _centerAfter(self, task, idx, upto):
        key = (tuple(task.get('helices', ())), int(idx), int(upto))
        cached = self._centerCache.get(key)
        if cached is not None:
            return cached
        representative = next((number for number in task.get('helices', ())
                               if number in self._data), -1)
        point = self._applyPipeline(self._taskCenterRaw(task, idx),
                                    representative, idx, upto)
        self._centerCache[key] = point
        return point

    def _localFrameBefore(self, task, idx, upto):
        before = self._centerAfter(task, idx-1, upto)
        after = self._centerAfter(task, idx+1, upto)
        tangent = self._unit(self._vsub(after, before))
        direction = math.radians(float(task.get('direction', 0.0)))
        raw_center = self._taskCenterRaw(task, idx)
        raw_direction = self._vadd(
            raw_center, (math.cos(direction), -math.sin(direction), 0.0))
        representative = next((number for number in task.get('helices', ())
                               if number in self._data), -1)
        moved_direction = self._applyPipeline(raw_direction, representative,
                                              idx, upto)
        center = self._centerAfter(task, idx, upto)
        normal = self._vsub(moved_direction, center)
        normal = self._vsub(normal,
                            self._vscale(tangent, self._dot(normal, tangent)))
        normal = self._unit(normal, (1.0, 0.0, 0.0))
        side = self._unit(self._cross(tangent, normal), (0.0, 1.0, 0.0))
        return center, tangent, normal, side

    def _bendGeometry(self, task, idx, task_index):
        """Integrate curvature along the already-twisted material frame."""
        key = ('bend', task_index, int(idx))
        cached = self._bendCache.get(key)
        if cached is not None:
            return cached
        start, end = int(task['start']), int(task['end'])
        target = min(max(start, int(idx)), end)
        if target == start:
            center, tangent, unused_normal, unused_side = \
                self._localFrameBefore(task, start, task_index)
            result = (center, tangent, [])
        else:
            center, tangent, rotations = self._bendGeometry(
                task, target-1, task_index)
            rotations = list(rotations)
            unused_c, prior_tangent, material_normal, unused_s = \
                self._localFrameBefore(task, target-1, task_index)
            transport_axis = self._cross(prior_tangent, tangent)
            transport_sine = math.sqrt(max(
                0.0, self._dot(transport_axis, transport_axis)))
            if transport_sine > 1e-10:
                transport_angle = math.atan2(
                    transport_sine,
                    max(-1.0, min(1.0,
                        self._dot(prior_tangent, tangent))))
                material_normal = self._rotateAxis(
                    material_normal, transport_axis, transport_angle)
            material_normal = self._vsub(
                material_normal,
                self._vscale(tangent, self._dot(material_normal, tangent)))
            material_normal = self._unit(material_normal,
                                         (1.0, 0.0, 0.0))
            axis = self._unit(self._cross(tangent, material_normal),
                              (0.0, 1.0, 0.0))
            step_angle = (math.radians(float(task.get('angle', 0.0))) /
                          float(max(1, end-start)))
            midpoint_tangent = self._rotateAxis(
                tangent, axis, step_angle*.5)
            center = self._vadd(center,
                                self._vscale(midpoint_tangent, .34))
            tangent = self._rotateAxis(tangent, axis, step_angle)
            rotations.append((axis, step_angle))
            result = (center, tangent, rotations)
        if idx > end:
            center, tangent, rotations = result
            result = (self._vadd(center,
                                 self._vscale(tangent, (idx-end)*.34)),
                      tangent, rotations)
        self._bendCache[key] = result
        return result

    def _materialTwistAngleBefore(self, number, idx, upto):
        """Return accumulated material rotation before pipeline item *upto*.

        Ordinary Add/Baseline transforms are incremental.  A
        ``remove_twist_final`` transform instead replaces the material twist
        accumulated across its selected interval with an absolute residual
        angle.  Tracking that state explicitly avoids relying on cancellation
        between transforms whose selected ranges may not have identical
        boundaries.
        """
        total = 0.0
        for task_index, task in enumerate(self._pipeline[:int(upto)]):
            if number not in task.get('helices', ()):
                continue
            mode = task.get('mode')
            if mode not in ('add_twist', 'remove_twist', 'baseline_twist',
                            'remove_twist_final'):
                continue
            start, end = int(task['start']), int(task['end'])
            if idx < start or end <= start:
                continue
            fraction = min(1.0, max(0.0,
                           (idx-start)/float(end-start)))
            if mode != 'remove_twist_final':
                total += float(task.get('angle', 0.0))*fraction
                continue
            target = float(task.get('target_angle', task.get('angle', 0.0)))
            start_before = self._materialTwistAngleBefore(
                number, start, task_index)
            end_before = self._materialTwistAngleBefore(
                number, end, task_index)
            if idx <= end:
                total = start_before + target*fraction
            else:
                # Preserve all twist outside the selected interval, changing
                # only the downstream phase by the interval correction.
                total += start_before + target-end_before
        return total

    def _applyOne(self, point, number, idx, task, task_index):
        if number not in task.get('helices', ()) or idx < int(task['start']):
            return point
        start, end = int(task['start']), int(task['end'])
        if end <= start:
            return point
        fraction = min(1.0, max(0.0, (idx-start)/float(end-start)))
        mode = task.get('mode')
        if mode == 'remove_twist_final':
            prior = self._materialTwistAngleBefore(number, idx, task_index)
            at_start = self._materialTwistAngleBefore(
                number, start, task_index)
            target = float(task.get('target_angle', task.get('angle', 0.0)))
            if idx <= end:
                desired = at_start + target*fraction
            else:
                at_end = self._materialTwistAngleBefore(
                    number, end, task_index)
                desired = prior + at_start + target-at_end
            correction = desired-prior
            center, tangent, unused_normal, unused_side = \
                self._localFrameBefore(task, idx, task_index)
            return self._vadd(
                center, self._rotateAxis(self._vsub(point, center),
                                         tangent,
                                         math.radians(correction)))
        if mode in ('add_twist', 'remove_twist', 'baseline_twist'):
            center, tangent, unused_normal, unused_side = \
                self._localFrameBefore(task, idx, task_index)
            angle = math.radians(float(task.get('angle', 0.0))*fraction)
            return self._vadd(
                center, self._rotateAxis(self._vsub(point, center),
                                         tangent, angle))
        if mode != 'bend' or abs(float(task.get('angle', 0.0))) < 1e-9:
            return point
        current_center = self._centerAfter(task, idx, task_index)
        offset = self._vsub(point, current_center)
        bent_center, unused_tangent, rotations = self._bendGeometry(
            task, idx, task_index)
        for axis, rotation in rotations:
            offset = self._rotateAxis(offset, axis, rotation)
        return self._vadd(bent_center, offset)

    def _applyPipeline(self, point, number, idx, upto=None):
        limit = len(self._pipeline) if upto is None else int(upto)
        result = point
        for task_index, task in enumerate(self._pipeline[:limit]):
            if (self._simulationPaths is not None and
                    task.get('source') == 'simulation_measurement'):
                continue
            result = self._applyOne(result, number, idx, task, task_index)
        return result

    def _rotate(self, point):
        x, y, z = point
        cy, sy = math.cos(self._yaw), math.sin(self._yaw)
        x, z = x * cy + z * sy, -x * sy + z * cy
        cp, sp = math.cos(self._pitch), math.sin(self._pitch)
        return x, y * cp - z * sp, y * sp + z * cp

    def _buildGeometryCache(self):
        """Calculate deformed world coordinates once per target pipeline.

        Camera rotation and zoom must never reintegrate Bend or reapply Twist.
        They project these cached world-space paths instead.
        """
        self._centerCache = {}
        self._bendCache = {}
        endpoints, reference_points = [], []
        world_paths = []
        if self._simulationPaths:
            for number, measured_samples in self._simulationPaths.items():
                if not measured_samples:
                    continue
                for point, unused_idx in measured_samples:
                    reference_points.append(point)
                    endpoints.append(point)
                samples = [(self._applyPipeline(point, number, idx), idx)
                           for point, idx in measured_samples]
                world_paths.append((number, tuple(samples)))
        else:
            for number, data in self._data.items():
                low, high = data['low'], data['high']
                for idx in (low, high):
                    axis = (data['coord'][0], data['coord'][1],
                            self._zForIndex(data, idx))
                    reference_points.append(axis)
                    endpoints.append(self._applyPipeline(axis, number, idx))
                step = max(1, int((high-low+1)/90))
                indices = set(range(low, high+1, step))
                indices.update((low, high))
                if low <= self._start <= high:
                    indices.add(self._start)
                if low <= self._end <= high:
                    indices.add(self._end)
                samples = []
                for idx in sorted(indices):
                    raw = (data['coord'][0], data['coord'][1],
                           self._zForIndex(data, idx))
                    samples.append((self._applyPipeline(raw, number, idx), idx))
                if samples:
                    world_paths.append((number, tuple(samples)))
        center = tuple(
            (min(point[dimension] for point in endpoints) +
             max(point[dimension] for point in endpoints)) / 2.0
            for dimension in range(3))
        extent = max(
            max(point[dimension] for point in reference_points) -
            min(point[dimension] for point in reference_points)
            for dimension in range(3)) or 1.0
        self._geometryCache = {
            'center': center, 'extent': extent,
            'world_paths': tuple(world_paths)}
        return self._geometryCache

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#eef0f2'))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        if not self._data:
            return
        geometry = self._geometryCache or self._buildGeometryCache()
        center = geometry['center']
        # Keep magnification invariant while a deformation angle changes.
        # Auto-fit uses the undeformed design only; the user's wheel zoom is
        # the sole operation that changes scale after that.
        extent = geometry['extent']
        scale = min(self.width(), self.height()) * .68 / extent * self._zoom
        rod_width = max(2.0, 2.0 * scale)  # physical dsDNA diameter: 2 nm
        outline_width = rod_width + max(1.2, .18 * scale)

        def project(point):
            rotated = self._rotate(tuple(point[d] - center[d] for d in range(3)))
            return self.width()/2 + rotated[0]*scale, \
                self.height()/2 - rotated[1]*scale, rotated[2]

        def make_path(samples):
            path = QPainterPath()
            if not samples:
                return path
            path.moveTo(samples[0][0][0], samples[0][0][1])
            for point, unused_idx in samples[1:]:
                path.lineTo(point[0], point[1])
            return path

        def draw_path(path, color, width):
            pen = QPen(color)
            pen.setWidthF(width)
            # One QPainterPath gives a single capsule. Round caps now occur
            # only at its two ends, never at every sampled base/module.
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(path)

        # Only project the cached world-space paths. No physical deformation
        # calculation occurs while the user rotates or zooms the camera.
        overlay = []
        for number, world_samples in geometry['world_paths']:
            samples = [(project(point), idx) for point, idx in world_samples]
            if samples:
                overlay.append((sum(item[0][2] for item in samples) /
                                len(samples), number, samples))

        for unused_depth, number, samples in sorted(overlay,
                                                     key=lambda item: item[0]):
            full_path = make_path(samples)
            draw_path(full_path, QColor(8, 12, 16, 245), outline_width)
            draw_path(full_path, QColor('#aeb4ba'), rod_width)
            draw_path(full_path, QColor('#f5f6f7'),
                      max(1.0, rod_width * .52))
            if number in self._selected:
                selected_samples = [item for item in samples
                                    if self._start <= item[1] <= self._end]
                if len(selected_samples) >= 2:
                    selected_path = make_path(selected_samples)
                    draw_path(selected_path, QColor(8, 12, 16, 245),
                              outline_width)
                    draw_path(selected_path, QColor('#269db8'), rod_width)
                    draw_path(selected_path, QColor('#a9edf6'),
                              max(1.0, rod_width * .52))

        # Draw helix numbers last.  Use a solid black second pass so the label
        # cannot become visually white on the light rods.
        label_font = QFont()
        label_font.setBold(True)
        label_font.setPointSize(12)
        halo_pen = QPen(QColor(255, 255, 255, 225))
        halo_pen.setWidthF(1.6)
        halo_pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        for unused_depth, number, samples in sorted(
                overlay, key=lambda item: item[0]):
            endpoint = samples[0][0]
            label_path = QPainterPath()
            label_path.addText(endpoint[0] + rod_width*.55 + 4,
                               endpoint[1] - 4, label_font, str(number))
            painter.setPen(halo_pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(label_path)
            black_pen = QPen(QColor('#050505'))
            black_pen.setWidthF(.7)
            painter.setPen(black_pen)
            painter.setBrush(QBrush(QColor('#050505')))
            painter.drawPath(label_path)

        painter.setPen(QColor('#263746'))
        painter.drawText(
            12, 22,
            '累计目标结构：当前设计残余 twist ＋ 已加入任务 ＋ 当前预览')
        painter.drawText(12, self.height()-12, '拖动旋转  ·  滚轮缩放')

    def mousePressEvent(self, event):
        self._last, self._drag = event.pos(), True

    def mouseMoveEvent(self, event):
        if not self._drag:
            return
        delta, self._last = event.pos() - self._last, event.pos()
        self._yaw += delta.x() * .009
        self._pitch += delta.y() * .009
        self._pitch = ((self._pitch + math.pi) % (2.0*math.pi)) - math.pi
        self.update()

    def mouseReleaseEvent(self, event):
        self._drag = False

    def alignCrossSection(self):
        """Match the screen orientation of the left lattice/helix panel."""
        self._yaw = 0.0
        # Looking back along the helix axis flips projected Y into the same
        # screen-down convention used by the left Qt lattice panel.
        self._pitch = math.pi
        self.update()

    def wheelEvent(self, event):
        self._zoom = max(.2, min(8.0, self._zoom *
                         math.pow(1.0015, event.angleDelta().y())))
        self.update()


class TwistBendDialog(QDialog):
    """Configure multiple deformation regions and return an executable plan."""

    def __init__(self, document, part, parent=None):
        super(TwistBendDialog, self).__init__(parent)
        self._document, self._part = document, part
        self._data = _helix_data(part)
        self._savedMetadata = document.twistBendMetadata() or {}
        self._hasRecordedTasks = bool(self._savedMetadata.get('tasks'))
        self._previewBaseData = (
            _preview_data_before_recorded_tasks(
                self._data, self._savedMetadata)
            if self._hasRecordedTasks else deepcopy(self._data))
        self._tasks, self._plans, self._pendingPlans = [], [], []
        self._previewTransforms = []
        self._previewWorkingData = deepcopy(self._data)
        self._activeDraftTransform = None
        self._draftTaskIndex = None
        self._syncingParameters = False
        self._applyingParameters = False
        self._parametersDirty = True
        self._appliedParameterTask = None
        self._inputSource = {'add_twist': 'angle', 'bend': 'angle'}
        self._simulationMeasurement = deepcopy(
            self._savedMetadata.get('simulation_measurement'))
        self.setWindowTitle('Twist and Bend')
        self.resize(1220, 790)
        root = QVBoxLayout(self)
        title = QLabel('<b>Twist and Bend</b>　普通单点阵设计 · SNUPI 校正预测')
        title.setStyleSheet('font-size: 16px;')
        root.addWidget(title)
        self.stepTitle = QLabel()
        self.stepTitle.setWordWrap(True)
        self.stepTitle.setStyleSheet(
            'font-size:14px; color:#173f5f; background:#eaf3fb; padding:9px;')
        root.addWidget(self.stepTitle)
        self.steps = QStackedWidget()
        root.addWidget(self.steps, 1)

        def scroll_page(widget):
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            scroll.setHorizontalScrollBarPolicy(
                Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll.setWidget(widget)
            return scroll

        # Step 1 — helix selection.
        helix_page = QWidget()
        helix_page_layout = QVBoxLayout(helix_page)
        helix_intro = QLabel(
            '单击选择单个 helix；在空白区域按住鼠标拖出矩形，可一次圈选大范围 helix。'
            '按住 Command/Ctrl/Shift 可将新圈选区域追加到当前选择。')
        helix_intro.setWordWrap(True)
        helix_page_layout.addWidget(helix_intro)
        helix_box = QGroupBox('选择 helix（可多选、可拖框）')
        helix_layout = QVBoxLayout(helix_box)
        self.helixPicker = HelixPicker(part, self._data)
        self.helixPicker.setMinimumHeight(430)
        helix_layout.addWidget(self.helixPicker)
        helix_page_layout.addWidget(helix_box, 1)
        self.steps.addWidget(scroll_page(helix_page))

        # Step 2 — inclusive base range.
        range_page = QWidget()
        range_page_layout = QVBoxLayout(range_page)
        range_intro = QLabel(
            '在标尺上拖动选择区域，也可以直接输入 Start/End。长度按首尾均包含计算，'
            '并实时计入当前设计已有的 insertion/deletion。')
        range_intro.setWordWrap(True)
        range_page_layout.addWidget(range_intro)
        range_box = QGroupBox('定义碱基区域（首尾均包含）')
        range_layout = QVBoxLayout(range_box)
        self.rangePicker = RangePicker(part.maxBaseIdx(), self._data)
        self.rangePicker.setMinimumHeight(180)
        range_layout.addWidget(self.rangePicker)
        row = QHBoxLayout()
        self.startSpin, self.endSpin = QSpinBox(), QSpinBox()
        for spin in (self.startSpin, self.endSpin):
            spin.setRange(0, part.maxBaseIdx())
        self.endSpin.setValue(part.indexOfRightmostNonemptyBase())
        row.addWidget(QLabel('Start base'))
        row.addWidget(self.startSpin)
        row.addWidget(QLabel('End base'))
        row.addWidget(self.endSpin)
        self.actualLabel = QLabel()
        self.actualLabel.setWordWrap(True)
        row.addWidget(self.actualLabel, 1)
        range_layout.addLayout(row)
        range_page_layout.addWidget(range_box)
        range_page_layout.addStretch(1)
        self.steps.addWidget(scroll_page(range_page))

        # Step 3 — measurement plus three independent deformation channels.
        mode_page = QWidget()
        mode_page_layout = QVBoxLayout(mode_page)
        mode_box = QGroupBox('选择并设置变形方式')
        mode_layout = QVBoxLayout(mode_box)
        self.channels = QTabWidget()
        mode_layout.addWidget(self.channels)
        analysis_page = QWidget()
        clear_page = QWidget()
        remove_page, add_page, bend_page = QWidget(), QWidget(), QWidget()
        self.channels.addTab(analysis_page, 'Accurate Twist/Bending')
        self.channels.addTab(clear_page, 'Remove Insertion/Deletion')
        self.channels.addTab(remove_page, 'Remove Twist')
        self.channels.addTab(add_page, 'Add Twist')
        self.channels.addTab(bend_page, 'Add Bending')

        analysis_form = QFormLayout(analysis_page)
        analysis_title = QLabel(
            '<b>准确计算 Twist/Bending：上传当前设计的模拟坐标</b>')
        analysis_title.setWordWrap(True)
        analysis_form.addRow(analysis_title)
        analysis_note = QLabel(
            '使用上方选中的 helix 和碱基区域进行计算。软件先读取 JSON 的 indel '
            '横截面分布：无 indel 或各 helix 完全均匀时使用整体最佳轴；只有明确的'
            '设计性 bending 才采用平滑 Bishop frame 扣除弯曲坐标架旋转。模拟中'
            '观察到的中心线弯曲始终单独报告。计算成功后，'
            '实测 twist（°/base）会自动成为 '
            'Remove Twist、Add Twist 和 Add Bending 的共同起始基线。')
        analysis_note.setWordWrap(True)
        analysis_form.addRow(analysis_note)
        self.loadSimulationButton = QPushButton('选择模拟坐标并准确计算')
        self.loadSimulationButton.setToolTip(
            '支持带 helix 分块标签的 XYZ、本软件导出的 PDB/mmCIF，以及同名 '
            'DAT+TOP。映射会按当前 JSON 的双链碱基、insertion 和 deletion 严格'
            '核对；外部结构不能唯一映射时会拒绝而不会猜测。')
        analysis_form.addRow(self.loadSimulationButton)
        self.clearSimulationButton = QPushButton('清除实测基线')
        self.clearSimulationButton.setEnabled(
            self._simulationMeasurement is not None)
        analysis_form.addRow(self.clearSimulationButton)
        self.simulationStatus = QLabel()
        self.simulationStatus.setWordWrap(True)
        self.simulationStatus.setStyleSheet(
            'color:#475569; background:#f2f5f7; padding:8px; '
            'border:1px solid #d7dfe5;')
        analysis_form.addRow(self.simulationStatus)

        clear_form = QFormLayout(clear_page)
        clear_title = QLabel(
            '<b>Remove Insertion/Deletion：移除选区内已有的 indel</b>')
        clear_title.setWordWrap(True)
        clear_form.addRow(clear_title)
        clear_note = QLabel(
            '按 helix 和 base 范围精确列出当前 insertion/deletion。可以只移除 '
            'insertion、只移除 deletion，或同时移除两者；不会改动选区外的 indel。')
        clear_note.setWordWrap(True)
        clear_form.addRow(clear_note)
        self.clearInsertionsCheck = QCheckBox('移除 insertion（正值）')
        self.clearInsertionsCheck.setChecked(True)
        self.clearDeletionsCheck = QCheckBox('移除 deletion/skip（负值）')
        self.clearDeletionsCheck.setChecked(True)
        clear_form.addRow(self.clearInsertionsCheck)
        clear_form.addRow(self.clearDeletionsCheck)
        self.existingIndelStatus = QLabel()
        self.existingIndelStatus.setWordWrap(True)
        self.existingIndelStatus.setStyleSheet(
            'color:#334e68; background:#f4f8fb; padding:7px; '
            'border:1px solid #c8d8e6;')
        clear_form.addRow('当前选区已有 indel', self.existingIndelStatus)

        remove_form = QFormLayout(remove_page)
        remove_title = QLabel('<b>Remove Twist：目标为预测宏观 twist 最接近 0</b>')
        remove_title.setWordWrap(True)
        remove_form.addRow(remove_title)
        remove_note = QLabel(
            '对全部选中的 helix 和区域进行联合整数优化。软件同时考虑当前螺距、'
            '截面刚度、crossover 连通性与选区长度；deletion 数量和位置按预测宏观 '
            'twist 最接近零且尽量均匀的原则分配。点击 Apply Parameters 后，'
            '3D 预览直接显示预测的 Remove Twist 最终宏观构型。')
        remove_note.setWordWrap(True)
        remove_form.addRow(remove_note)
        self.removeBaselineLabel = QLabel()
        self.removeBaselineLabel.setWordWrap(True)
        self.removeBaselineLabel.setStyleSheet(
            'font-weight:700; color:#173f5f; background:#f4f8fb; '
            'border:1px solid #b8cfdf; padding:6px;')
        baseline_title = QLabel('<b>当前 Twist 基线</b>')
        remove_form.addRow(baseline_title)
        # Use the full form width.  In the previous two-column row the long
        # measured/calibrated baseline description was clipped on narrower
        # screens even though word wrapping was enabled.
        remove_form.addRow(self.removeBaselineLabel)
        add_form = QFormLayout(add_page)
        self.twistAngle = QDoubleSpinBox()
        self.twistAngle.setRange(0, 2160)
        self.twistAngle.setValue(90)
        self.twistAngle.setSuffix('°')
        self.handedness = QComboBox()
        self.handedness.addItems(['Right-handed', 'Left-handed'])
        self.twistPitch = QDoubleSpinBox()
        self.twistPitch.setRange(5.0, 20.0)
        self.twistPitch.setDecimals(4)
        self.twistPitch.setSingleStep(0.01)
        self.twistPitch.setValue(float(part.helicalPitch()))
        self.twistPitch.setSuffix(' base/turn')
        self.twistPitch.setToolTip(
            '输入目标平均螺距；软件会换算为最接近的整数 addition/deletion 分配。')
        self.currentTwistPitch = QLabel()
        self.currentTwistPitch.setWordWrap(True)
        self.currentTwistPitch.setMinimumWidth(170)
        self.currentTwistPitch.setText('请先选择 helix 和碱基区域')
        self.currentTwistPitch.setStyleSheet(
            'font-weight:700; color:#173f5f; background:#f4f8fb; '
            'border:1px solid #b8cfdf; padding:5px;')
        self.twistPhysicalLabel = QLabel('请先选择 helix 和碱基区域')
        self.twistPhysicalLabel.setWordWrap(True)
        self.twistPhysicalLabel.setStyleSheet(
            'font-weight:700; color:#173f5f; background:#f4f8fb; '
            'border:1px solid #b8cfdf; padding:5px;')
        self.twistIndels = QDoubleSpinBox()
        self.twistIndels.setRange(-100.0, 100.0)
        self.twistIndels.setDecimals(3)
        self.twistIndels.setSingleStep(0.1)
        self.twistIndels.setValue(90.0 / 360.0 * TARGET_PITCH)
        self.twistIndels.setSuffix(' base/helix')
        self.twistIndels.setToolTip(
            '正数为 addition（右手），负数为 deletion（左手）。小数会在全部选中 helix 间均匀取整分配。')
        self.twistDial = AngleDial('拖动设置末端相对角度')
        self.twistDial.setValue(90)
        add_form.addRow('预测/目标宏观 twist', self.twistAngle)
        add_form.addRow('旋转方向', self.handedness)
        add_form.addRow('当前选区平均螺距', self.currentTwistPitch)
        add_form.addRow('目标平均螺距', self.twistPitch)
        add_form.addRow('平均 addition/deletion', self.twistIndels)
        add_form.addRow('物理预测', self.twistPhysicalLabel)
        add_form.addRow(self.twistDial)
        add_note = QLabel(
            '角度为基于文献数据、截面刚度与 crossover 连通性估算，并经 SNUPI Static '
            '校正的整体宏观 twist。'
            '拖动圆盘或输入目标角度后，软件反求最接近的整数 indel；也可直接输入螺距'
            '或平均 indel。结果属于快速经验预测，不等同于有限元或分子动力学模拟。')
        add_note.setWordWrap(True)
        add_form.addRow(add_note)

        bend_form = QFormLayout(bend_page)
        self.bendAngle, self.bendDirection = QDoubleSpinBox(), QDoubleSpinBox()
        self.bendAngle.setRange(0, 720)
        self.bendAngle.setValue(45)
        self.bendAngle.setSuffix('°')
        self.bendDirection.setRange(0, 359.9)
        self.bendDirection.setValue(0)
        self.bendDirection.setSuffix('°')
        self.bendIndels = QDoubleSpinBox()
        self.bendIndels.setRange(0.0, 100.0)
        self.bendIndels.setDecimals(3)
        self.bendIndels.setSingleStep(0.1)
        self.bendIndels.setSuffix(' |base|/helix')
        self.bendIndels.setToolTip(
            '平均绝对 indel 强度；内侧自动 deletion、外侧自动 addition，并按截面距离成比例分配。')
        self.bendAngleSlider = QSlider(Qt.Orientation.Horizontal)
        self.bendAngleSlider.setRange(0, 360)
        self.bendAngleSlider.setValue(45)
        self.bendAngleSlider.setTickInterval(15)
        self.bendAngleSlider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.bendDirectionDial = BendDirectionPicker(self._data)
        self.bendDirectionDial.setValue(0)
        bend_form.addRow('总弯曲角', self.bendAngle)
        bend_form.addRow('平均 addition/deletion', self.bendIndels)
        bend_form.addRow('拖动调整弯曲角', self.bendAngleSlider)
        bend_form.addRow('截面弯曲方向', self.bendDirection)
        bend_form.addRow(self.bendDirectionDial)
        self.bendSafetyLabel = QLabel()
        self.bendSafetyLabel.setWordWrap(True)
        self.bendSafetyLabel.setStyleSheet('color:#8a5a00;')
        bend_form.addRow(self.bendSafetyLabel)
        bend_note = QLabel('0° 向右，90° 向上；蓝色箭头指向弯曲内侧。使用实际 '
                           'Square/Honeycomb 坐标和 2.8 nm helix 中心间距计算中性轴，'
                           '箭头侧 deletion、相反的外侧 insertion。')
        bend_note.setWordWrap(True)
        bend_form.addRow(bend_note)
        # Keep the live result immediately below the inputs of each channel,
        # rather than below the entire tab widget.
        self.livePlanSummaries = [self.simulationStatus]
        for form in (clear_form, remove_form, add_form, bend_form):
            label = QLabel('当前参数尚未形成可计算方案。')
            label.setWordWrap(True)
            label.setStyleSheet(
                'color:#284b63; background:#edf5fa; padding:8px; border:1px solid #cadde9;')
            form.addRow(label)
            self.livePlanSummaries.append(label)
        apply_row = QHBoxLayout()
        apply_row.addStretch(1)
        self.applyParameters = QPushButton(
            'Apply Parameters / 更新计算与预览')
        self.applyParameters.setToolTip(
            '根据当前输入计算联动参数、物理预测、indel 方案和完整 3D 目标结构。')
        self.applyParameters.setMinimumWidth(250)
        apply_row.addWidget(self.applyParameters)
        mode_layout.addLayout(apply_row)
        self.applyParameterStatus = QLabel(
            '参数尚未计算。设置任意参数后，请点击 Apply Parameters。')
        self.applyParameterStatus.setWordWrap(True)
        self.applyParameterStatus.setStyleSheet(
            'color:#8a5a00; background:#fff6dc; padding:7px; '
            'border:1px solid #ead39a;')
        mode_layout.addWidget(self.applyParameterStatus)
        mode_page_layout.addWidget(mode_box)
        live_box = QGroupBox('当前参数的 3D 预览')
        live_layout = QVBoxLayout(live_box)
        self.parameterPreview = DeformationPreview(self._previewBaseData)
        self.parameterPreview.setMinimumHeight(390)
        live_layout.addWidget(self.parameterPreview)
        live_hint = QLabel(
            '拖动 3D 视图可旋转，滚轮缩放；参数改变后点击 Apply Parameters 更新。')
        live_hint.setWordWrap(True)
        live_layout.addWidget(live_hint)
        mode_page_layout.addWidget(live_box)
        mode_page_layout.addStretch(1)
        self.steps.addWidget(scroll_page(mode_page))

        # Step 4 — full-width 3D review, task list, and Apply.
        preview_page = QWidget()
        preview_page_layout = QVBoxLayout(preview_page)
        preview_box = QGroupBox('3D 整体预测')
        preview_layout = QVBoxLayout(preview_box)
        self.preview = DeformationPreview(self._previewBaseData)
        self.preview.setMinimumHeight(470)
        view_buttons = QHBoxLayout()
        self.alignCrossSectionButton = QPushButton('对齐左侧截面视角')
        self.alignCrossSectionButton.setToolTip(
            '沿 helix 轴观察，并保持 helix 的上下、左右顺序与左侧截面完全一致。')
        view_buttons.addWidget(self.alignCrossSectionButton)
        view_buttons.addStretch(1)
        preview_layout.addLayout(view_buttons)
        preview_layout.addWidget(self.preview, 1)
        legend = QLabel(
            '灰白圆柱：每根 dsDNA（直径 2.0 nm）　青蓝色：选中的碱基区域　'
            'helix 中心间距：2.8 nm　轴向长度按实际碱基数计算')
        legend.setWordWrap(True)
        legend.setStyleSheet('color:#61758a;')
        preview_layout.addWidget(legend)
        preview_page_layout.addWidget(preview_box, 1)
        task_row = QHBoxLayout()
        self.addTask = QPushButton('更新当前区域')
        self.addAnotherTask = QPushButton('＋ 再添加一个区域')
        self.removeTask = QPushButton('删除所选任务')
        task_row.addWidget(self.addTask)
        task_row.addWidget(self.addAnotherTask)
        task_row.addWidget(self.removeTask)
        preview_page_layout.addLayout(task_row)
        self.taskList = QListWidget()
        self.taskList.setMaximumHeight(115)
        preview_page_layout.addWidget(self.taskList)
        self.summary = QLabel('尚未添加任务')
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet('color:#334e68; background:#eef5fb; padding:7px;')
        preview_page_layout.addWidget(self.summary)
        warning = QLabel('这是几何/弹性模型经 SNUPI Static 校正后的目标预览，并非分子动力学弛豫结果。Apply 不会删除或移动 crossover；'
                         '所有新增 indel 会避开现有 crossover、nick 和 indel。')
        warning.setWordWrap(True)
        warning.setStyleSheet('background:#fff5d7; color:#604b00; padding:8px;')
        preview_page_layout.addWidget(warning)
        self.steps.addWidget(scroll_page(preview_page))

        navigation = QHBoxLayout()
        self.backButton = QPushButton('← 上一步')
        self.nextButton = QPushButton('下一步 →')
        self.stepCounter = QLabel()
        self.stepCounter.setAlignment(Qt.AlignmentFlag.AlignCenter)
        navigation.addWidget(self.backButton)
        navigation.addWidget(self.stepCounter, 1)
        navigation.addWidget(self.nextButton)
        root.addLayout(navigation)

        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel |
                                        QDialogButtonBox.StandardButton.Ok)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText(
            'Apply（单步 Undo）')
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(False)
        root.addWidget(self.buttons)

        # A persistent single-page task editor makes it possible to add a
        # Twist task, switch channels, add a Bend task, and inspect both in
        # one place. The left side scrolls independently on smaller screens.
        self.stepTitle.hide()
        self.steps.hide()
        self.backButton.hide()
        self.nextButton.hide()
        self.stepCounter.hide()
        live_box.hide()  # the larger, always-visible preview is on the right
        self.addAnotherTask.hide()
        self.addTask.setText('＋ 添加当前任务')
        self.helixPicker.setMinimumHeight(260)
        self.rangePicker.setMinimumHeight(115)

        editor = QWidget()
        editor_layout = QHBoxLayout(editor)
        editor_layout.setContentsMargins(0, 0, 0, 0)
        controls_widget = QWidget()
        controls_layout = QVBoxLayout(controls_widget)
        controls_layout.setContentsMargins(4, 4, 8, 4)
        controls_layout.addWidget(helix_box)
        controls_layout.addWidget(range_box)
        controls_layout.addWidget(mode_box)
        task_box = QGroupBox('任务列表（可同时包含 Twist 和 Bend）')
        task_box_layout = QVBoxLayout(task_box)
        task_buttons = QHBoxLayout()
        task_buttons.addWidget(self.addTask)
        task_buttons.addWidget(self.removeTask)
        task_box_layout.addLayout(task_buttons)
        task_box_layout.addWidget(self.taskList)
        task_box_layout.addWidget(self.summary)
        controls_layout.addWidget(task_box)
        controls_layout.addStretch(1)
        controls_scroll = QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setFrameShape(QFrame.Shape.NoFrame)
        controls_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        controls_scroll.setMinimumWidth(500)
        controls_scroll.setWidget(controls_widget)
        editor_layout.addWidget(controls_scroll, 5)
        preview_column = QVBoxLayout()
        preview_column.addWidget(preview_box, 1)
        preview_column.addWidget(warning)
        editor_layout.addLayout(preview_column, 7)
        root.insertWidget(2, editor, 1)

        self.rangePicker.rangeChanged.connect(self._rangeDragged)
        self.startSpin.valueChanged.connect(self._controlsChanged)
        self.endSpin.valueChanged.connect(self._controlsChanged)
        self.helixPicker.selectionChanged.connect(self._controlsChanged)
        self.channels.currentChanged.connect(self._channelChanged)
        self.loadSimulationButton.clicked.connect(self._loadSimulationResult)
        self.clearSimulationButton.clicked.connect(
            self._clearSimulationMeasurement)
        self.clearInsertionsCheck.toggled.connect(self._controlsChanged)
        self.clearDeletionsCheck.toggled.connect(self._controlsChanged)
        self.twistAngle.valueChanged.connect(self._twistAngleChanged)
        self.twistDial.valueChanged.connect(self._twistDialChanged)
        self.handedness.currentIndexChanged.connect(self._twistAngleChanged)
        self.twistPitch.valueChanged.connect(self._twistPitchChanged)
        self.twistIndels.valueChanged.connect(self._twistIndelsChanged)
        self.bendAngle.valueChanged.connect(self._bendAngleChanged)
        self.bendAngle.valueChanged.connect(
            lambda value: self.bendAngleSlider.setValue(int(round(value))))
        self.bendAngleSlider.valueChanged.connect(self.bendAngle.setValue)
        self.bendIndels.valueChanged.connect(self._bendIndelsChanged)
        self.bendDirection.valueChanged.connect(self._bendDirectionChanged)
        self.bendDirectionDial.valueChanged.connect(
            self._bendDirectionDialChanged)
        self.applyParameters.clicked.connect(self._applyParameters)
        self.addTask.clicked.connect(self._addTask)
        self.addAnotherTask.clicked.connect(self._addAnotherTask)
        self.removeTask.clicked.connect(self._removeTask)
        self.alignCrossSectionButton.clicked.connect(
            self.preview.alignCrossSection)
        self.alignCrossSectionButton.clicked.connect(
            self.parameterPreview.alignCrossSection)
        self.backButton.clicked.connect(self._previousStep)
        self.nextButton.clicked.connect(self._nextStep)
        self.buttons.rejected.connect(self.reject)
        self.buttons.accepted.connect(self._accept)
        self._restoreMetadata()
        self._applySimulationGeometry()
        self._refreshSimulationStatus()
        self._controlsChanged()
        self.addTask.setEnabled(False)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            bool(self._tasks))

    def _showStep(self, index):
        index = max(0, min(3, int(index)))
        self.steps.setCurrentIndex(index)
        titles = (
            '第 1 步：选择要处理的 helix',
            '第 2 步：定义 helix 上的碱基起止区域',
            '第 3 步：选择 Remove Twist、Add Twist 或 Bend 并设置参数',
            '第 4 步：检查 3D 目标形状、indel 方案和任务汇总')
        self.stepTitle.setText(titles[index])
        self.stepCounter.setText('%d / 4' % (index + 1))
        self.backButton.setEnabled(index > 0)
        self.nextButton.setVisible(index < 3)
        self.addTask.setEnabled(index == 3 and
                                self._draftTaskIndex is not None)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            index == 3 and bool(self._tasks))

    def _previousStep(self):
        self._showStep(self.steps.currentIndex() - 1)

    def _nextStep(self):
        index = self.steps.currentIndex()
        if index == 0 and not self.helixPicker.selected():
            QMessageBox.information(self, 'Twist and Bend',
                                    '请先单击或拖框选择至少一个 helix。')
            return
        if index == 1 and self.startSpin.value() == self.endSpin.value():
            QMessageBox.information(self, 'Twist and Bend',
                                    '碱基区域至少需要包含两个位置。')
            return
        if index == 2:
            if not self._storeDraftTask():
                return
            self._showStep(3)
            return
        self._showStep(index + 1)

    def _restoreMetadata(self):
        metadata = self._savedMetadata
        tasks = metadata.get('tasks', [])
        saved_plans = metadata.get('last_plans', [])
        lattice = ('honeycomb' if self._part._step == 21 else 'square')
        upgrade_predictions = (
            metadata.get('prediction_calibration') !=
            SNUPI_CALIBRATION_VERSION)
        if tasks:
            first = tasks[0]
            self.helixPicker.setSelected(first.get('helices', []))
            self.startSpin.setValue(int(first.get('start', 0)))
            self.endSpin.setValue(int(first.get('end', self.endSpin.value())))
            mode = first.get('mode', 'remove_twist')
            self.channels.setCurrentIndex(
                {'remove_indels': 1, 'remove_twist': 2,
                 'add_twist': 3, 'bend': 4}.get(mode, 2))
            if mode == 'add_twist':
                self.twistAngle.setValue(float(first.get('angle', 90)))
                self.handedness.setCurrentIndex(
                    0 if first.get('handedness', 'right') == 'right' else 1)
                if 'indels_per_helix' in first:
                    self.twistIndels.setValue(
                        float(first['indels_per_helix']))
                if 'target_pitch' in first:
                    self.twistPitch.setValue(float(first['target_pitch']))
            elif mode == 'bend':
                self.bendAngle.setValue(float(first.get('angle', 45)))
                self.bendDirection.setValue(float(first.get('direction', 0)))
                if 'indels_per_helix' in first:
                    self.bendIndels.setValue(
                        abs(float(first['indels_per_helix'])))
        for index, task in enumerate(tasks):
            restored = dict(task)
            restored['_applied'] = True
            if index < len(saved_plans):
                saved_plan = dict(saved_plans[index])
                if upgrade_predictions:
                    saved_plan = calibrate_saved_plan(saved_plan, lattice)
                restored['_saved_plan'] = saved_plan
            self._appendTask(restored, rebuild=False)
        if tasks:
            self._rebuildPlans()

    def _rangeDragged(self, start, end):
        self.startSpin.blockSignals(True)
        self.endSpin.blockSignals(True)
        self.startSpin.setValue(start)
        self.endSpin.setValue(end)
        self.startSpin.blockSignals(False)
        self.endSpin.blockSignals(False)
        self._controlsChanged()

    def _currentMode(self):
        return ('simulation_analysis', 'remove_indels', 'remove_twist',
                'add_twist', 'bend')[
            self.channels.currentIndex()]

    def _channelChanged(self, unused=None):
        analysis = self._currentMode() == 'simulation_analysis'
        self.applyParameters.setEnabled(not analysis and
                                        self._parametersDirty)
        if analysis:
            self.addTask.setEnabled(False)
            self._refreshSimulationStatus()
        self._controlsChanged()

    def _selectedExistingIndels(self, selected=None, start=None, end=None):
        if selected is None:
            selected = self.helixPicker.selected()
        if start is None or end is None:
            start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        rows = []
        for number in sorted(selected):
            for idx, length in sorted(
                    self._data[number].get('insertions', {}).items()):
                if start <= int(idx) <= end:
                    rows.append((int(number), int(idx), int(length)))
        return rows

    def _refreshExistingIndels(self):
        rows = self._selectedExistingIndels()
        if not rows:
            self.existingIndelStatus.setText('无')
            return
        grouped = {}
        for number, idx, length in rows:
            grouped.setdefault(number, []).append(
                '%d:%+d' % (idx, length))
        self.existingIndelStatus.setText(
            '；'.join('H%d [%s]' % (number, ', '.join(values))
                     for number, values in sorted(grouped.items())))

    def _measurementMatches(self, selected=None, start=None, end=None):
        measurement = self._simulationMeasurement
        if not measurement:
            return False
        if selected is None:
            selected = self.helixPicker.selected()
        if start is None or end is None:
            start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        return (sorted(int(value) for value in selected) ==
                sorted(int(value) for value in measurement.get('helices', ()))
                and int(start) == int(measurement.get('start', -1))
                and int(end) == int(measurement.get('end', -1)))

    def _refreshSimulationStatus(self):
        measurement = self._simulationMeasurement
        if not measurement:
            self.simulationStatus.setText(
                '尚未加载模拟结果。Remove Twist 将使用当前弹性/校准预测基线。')
            self.simulationStatus.setStyleSheet(
                'color:#475569; background:#f2f5f7; padding:8px; '
                'border:1px solid #d7dfe5;')
            self.removeBaselineLabel.setText(
                '未加载匹配的模拟结果：使用当前弹性/校准预测。')
            return
        radius = measurement.get('bend_radius_nm')
        if radius is None and measurement.get(
                'bend_radius_coordinate_units') is not None:
            radius = (float(measurement['bend_radius_coordinate_units']) *
                      float(measurement.get('coordinate_scale_nm', 1.0)))
        radius_text = ('∞' if radius is None else '%.3f nm' % radius)
        match = self._measurementMatches()
        bend_design = measurement.get('design_bend_classification', {})
        bend_class = bend_design.get('classification', 'unknown')
        bend_class_text = {
            'none': '无设计性 bending',
            'bending': '存在设计性 bending',
            'complex': '复杂局部变形',
        }.get(bend_class, '旧结果/尚未分类')
        bishop_text = ('已启用平滑 Bishop frame 扣除' if measurement.get(
            'bishop_correction_applied') else '未启用 Bishop bending 扣除')
        self.simulationStatus.setText(
            '已加载：%s\nhelix %s；base %d–%d；实测 twist %+.6f°/base'
            '（Accumulated material twist %+.3f°；拟合 RMS %.3f°）。\n'
            '设计图判断：%s；%s。%s\n'
            'Observed bending analysis：角度 %.3f°，方向 %.2f°，曲率半径 %s。\n'
            '整体轴/Bishop Twist：%+.6f / %+.6f°/base。\n%s' %
            (measurement.get('source_path', ''),
             ','.join(str(value) for value in measurement.get('helices', ())),
             measurement.get('start', 0), measurement.get('end', 0),
             measurement.get('twist_per_base_deg', 0.0),
             measurement.get('total_twist_deg', 0.0),
             measurement.get('twist_fit_rms_deg', 0.0),
             bend_class_text, bishop_text, bend_design.get('reason', ''),
             measurement.get('bend_angle_deg', 0.0),
             measurement.get('bend_direction_deg', 0.0), radius_text,
             measurement.get('global_axis_twist_per_base_deg',
                             measurement.get('twist_per_base_deg', 0.0)),
             measurement.get('bishop_twist_per_base_deg',
                             measurement.get('twist_per_base_deg', 0.0)),
             ('当前选区匹配：Remove/Add Twist 与 Add Bending 将自动使用此实测基线。'
              if match else
              '当前选区不匹配：请恢复同一 helix/base 区域后再使用实测基线。')))
        self.simulationStatus.setStyleSheet(
            ('color:#234e34; background:#e9f7ef; padding:8px; '
             'border:1px solid #b9ddc5;' if match else
             'color:#8a5a00; background:#fff6dc; padding:8px; '
             'border:1px solid #ead39a;'))
        if match:
            self.removeBaselineLabel.setText(
                '实测 %+.6f°/base；Remove Twist 将以 0°/base 为目标反求 indel。' %
                measurement.get('twist_per_base_deg', 0.0))
        else:
            self.removeBaselineLabel.setText(
                '模拟结果与当前选区不匹配：暂时使用弹性/校准预测。')

    def _loadSimulationResult(self):
        selected = self.helixPicker.selected()
        start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        if len(selected) < 2:
            QMessageBox.information(
                self, 'Accurate Twist/Bending',
                '请先在上方选择至少两个 helix，并定义碱基区域。')
            return
        path, unused_filter = QFileDialog.getOpenFileName(
            self, '选择模拟坐标', '',
            'Labelled XYZ (*.xyz);;PDB/mmCIF (*.pdb *.cif *.mmcif);;oxDNA DAT (*.dat);;All files (*)')
        if not path:
            return
        try:
            measurement = analyze_simulation_file(
                path, self._document, self._data, selected, start, end)
        except (OSError, SimulationAnalysisError) as error:
            QMessageBox.warning(self, 'Accurate Twist/Bending', str(error))
            return
        self._simulationMeasurement = measurement
        self._applySimulationGeometry()
        self.clearSimulationButton.setEnabled(True)
        self._refreshSimulationStatus()
        # Put Remove Twist immediately after the measurement utilities and
        # make the newly measured baseline visible without another step.
        self.channels.setCurrentIndex(2)
        self._markParametersDirty()

    def _clearSimulationMeasurement(self):
        self._simulationMeasurement = None
        self._applySimulationGeometry()
        self.clearSimulationButton.setEnabled(False)
        self._refreshSimulationStatus()
        self._markParametersDirty()

    def _applySimulationGeometry(self):
        centerlines = ((self._simulationMeasurement or {}).get(
            'preview_centerlines'))
        self.preview.setSimulationGeometry(centerlines)
        self.parameterPreview.setSimulationGeometry(centerlines)

    def _twistDialChanged(self, angle):
        current = self.twistAngle.value()
        revolutions = int(current // 360.0)
        value = min(self.twistAngle.maximum(), revolutions * 360.0 + angle)
        self.twistAngle.setValue(value)

    def _mechanicalBaseIndices(self, number, start, end):
        """Return paired dsDNA indices, respecting staggered helix ends."""
        data = self._data[number]
        scaffold_intervals = data.get('scaffold_intervals', ())
        staple_intervals = data.get('staple_intervals', ())
        if scaffold_intervals or staple_intervals:
            scaffold = set()
            staple = set()
            for low, high in scaffold_intervals:
                scaffold.update(range(max(int(start), int(low)),
                                      min(int(end), int(high))+1))
            for low, high in staple_intervals:
                staple.update(range(max(int(start), int(low)),
                                   min(int(end), int(high))+1))
            return sorted(scaffold & staple)
        return list(range(int(start), int(end)+1))

    def _mechanicalLength(self, number, start, end):
        existing = self._data[number]['insertions']
        return sum(max(0, 1+int(existing.get(idx, 0)))
                   for idx in self._mechanicalBaseIndices(
                       number, start, end))

    def _twistPitchMetrics(self):
        selected = self.helixPicker.selected()
        start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        if not selected or end < start:
            return 0.0, 0, float(self._part.helicalPitch())
        native_total = sum(len(self._mechanicalBaseIndices(
            number, start, end)) for number in selected)
        turns = native_total / float(self._part.helicalPitch())
        actual = sum(self._mechanicalLength(number, start, end)
                     for number in selected)
        current_pitch = actual / turns if turns > 1e-12 else 0.0
        return turns, actual, current_pitch

    def _setTwistPitchForDelta(self, total_delta):
        turns, actual, current_pitch = self._twistPitchMetrics()
        self.currentTwistPitch.setText(
            ('%.4f base/turn' % current_pitch)
            if turns > 0 else '请先选择 helix 和碱基区域')
        if turns <= 0:
            return
        self.twistPitch.blockSignals(True)
        self.twistPitch.setValue((actual + int(total_delta)) / turns)
        self.twistPitch.blockSignals(False)

    def _twistRegion(self):
        start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        return {'mode': 'add_twist',
                'helices': self.helixPicker.selected(),
                'start': start, 'end': end}

    def _twistEstimateForDelta(self, total_delta):
        region = self._twistRegion()
        if not region['helices'] or region['end'] <= region['start']:
            return None
        prediction = estimate_global_twist(
            region, self._data, self._part.helicalPitch(),
            extra_base_delta=int(total_delta))
        if self._measurementMatches(region['helices'], region['start'],
                                    region['end']):
            baseline = estimate_global_twist(
                region, self._data, self._part.helicalPitch(),
                extra_base_delta=0)
            prediction = self._anchorTwistPrediction(
                prediction, baseline, self._simulationMeasurement)
        return prediction

    @staticmethod
    def _anchorTwistPrediction(prediction, model_baseline, measurement):
        """Apply a model-predicted indel increment to a measured baseline."""
        anchored = dict(prediction)
        measured = float(measurement['twist_per_base_deg'])
        delta = (prediction['twist_per_base_deg']-
                 model_baseline['twist_per_base_deg'])
        value = measured+delta
        anchored['model_twist_per_base_deg'] = prediction[
            'twist_per_base_deg']
        anchored['twist_per_base_deg'] = value
        anchored['total_twist_deg'] = value*prediction['length_bp']
        anchored['handedness'] = ('右手' if value > 1e-9 else
                                  '左手' if value < -1e-9 else
                                  '近似无扭转')
        anchored['measurement_anchored'] = True
        anchored['measured_baseline_twist_per_base_deg'] = measured
        anchored['confidence'] = '实测基线'
        anchored['note'] = (
            '基线来自上传模拟坐标；新增 indel 的响应使用最新校正模型。')
        return anchored

    def _setTwistPrediction(self, prediction):
        if prediction is None:
            self.twistPhysicalLabel.setText('请先选择 helix 和碱基区域')
            return
        self.twistPhysicalLabel.setText(
            ('%+.4f°/base；总计 %+.2f°（%s）；J=%.1f nm⁴；'
             'crossover 连通 %.0f%%；可信度：%s。%s') %
            (prediction['twist_per_base_deg'],
             prediction['total_twist_deg'], prediction['handedness'],
             prediction['polar_moment_nm4'],
             prediction['connectivity_fraction']*100.0,
             prediction['confidence'], prediction['note']))

    def _syncTwistFromDelta(self, total_delta):
        selected_count = max(1, len(self.helixPicker.selected()))
        prediction = self._twistEstimateForDelta(total_delta)
        self._syncingParameters = True
        self.twistIndels.blockSignals(True)
        self.twistAngle.blockSignals(True)
        self.handedness.blockSignals(True)
        self.twistIndels.setValue(int(total_delta) / float(selected_count))
        if prediction is not None:
            self.twistAngle.setValue(abs(prediction['total_twist_deg']))
            if abs(prediction['total_twist_deg']) > 1e-9:
                self.handedness.setCurrentIndex(
                    0 if prediction['total_twist_deg'] > 0 else 1)
        self.twistIndels.blockSignals(False)
        self.twistAngle.blockSignals(False)
        self.handedness.blockSignals(False)
        self._setTwistPitchForDelta(total_delta)
        self._setTwistPrediction(prediction)
        self._syncingParameters = False

    def _deltaForTargetTwist(self, target_angle):
        selected = self.helixPicker.selected()
        if not selected:
            return 0
        count = len(selected)
        base = self._twistEstimateForDelta(0)
        layer = self._twistEstimateForDelta(count)
        if base is None or layer is None:
            return 0
        slope = ((layer['total_twist_deg']-base['total_twist_deg']) /
                 float(count))
        if abs(slope) < 1e-9:
            return 0
        candidate = int(round((float(target_angle)-
                               base['total_twist_deg']) / slope))
        region = self._twistRegion()
        capacity = sum(len([
            idx for idx in range(region['start']+1, region['end'])
            if idx not in self._data[number].get('forbidden', set()) and
            idx not in self._data[number].get('insertions', {})])
            for number in selected)
        candidate = max(-capacity, min(capacity, candidate))
        radius = max(6, min(12, count//4 + 2))
        candidates = range(max(-capacity, candidate-radius),
                           min(capacity, candidate+radius)+1)
        return min(candidates, key=lambda delta: (
            abs(self._twistEstimateForDelta(delta)['total_twist_deg']-
                float(target_angle)), abs(delta)))

    def _twistAngleChanged(self, unused=None):
        if self._syncingParameters:
            return
        self._inputSource['add_twist'] = 'angle'
        self._controlsChanged()

    def _twistPitchChanged(self, value):
        if self._syncingParameters:
            return
        self._inputSource['add_twist'] = 'pitch'
        self._controlsChanged()

    def _twistIndelsChanged(self, value):
        if self._syncingParameters:
            return
        self._inputSource['add_twist'] = 'indels'
        self._controlsChanged()

    def _bendMeanOffset(self):
        selected = self.helixPicker.selected()
        if not selected:
            return 0.0
        direction = math.radians(self.bendDirection.value())
        ux, uy = math.cos(direction), -math.sin(direction)
        values = [self._data[number]['coord'][0] * ux +
                  self._data[number]['coord'][1] * uy for number in selected]
        neutral = (max(values) + min(values)) / 2.0
        return sum(abs(value-neutral) for value in values) / float(len(values))

    def _bendAngleChanged(self, unused=None):
        if self._syncingParameters:
            return
        self._inputSource['bend'] = 'angle'
        self._controlsChanged()

    def _bendIndelsChanged(self, value):
        if self._syncingParameters:
            return
        self._inputSource['bend'] = 'indels'
        self._controlsChanged()

    def _bendDirectionChanged(self, unused=None):
        if self._syncingParameters:
            return
        self._controlsChanged()

    def _bendDirectionDialChanged(self, angle):
        self.bendDirection.setValue(angle)

    def _currentTask(self):
        start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        mode = self._currentMode()
        task = {'mode': mode, 'helices': self.helixPicker.selected(),
                'start': start, 'end': end,
                'angle': (self.bendAngle.value() if mode == 'bend'
                          else self.twistAngle.value()),
                'direction': self.bendDirection.value(),
                'handedness': ('right' if self.handedness.currentIndex() == 0
                               else 'left')}
        if mode == 'remove_indels':
            task['remove_insertions'] = self.clearInsertionsCheck.isChecked()
            task['remove_deletions'] = self.clearDeletionsCheck.isChecked()
        elif mode == 'add_twist':
            task['indels_per_helix'] = self.twistIndels.value()
            task['target_pitch'] = self.twistPitch.value()
        elif mode == 'bend':
            task['indels_per_helix'] = self.bendIndels.value()
        elif mode == 'remove_twist' and self._measurementMatches(
                task['helices'], start, end):
            task['measured_twist_per_base_deg'] = float(
                self._simulationMeasurement['twist_per_base_deg'])
            task['measurement_source'] = self._simulationMeasurement.get(
                'source_path', '')
        return task

    @staticmethod
    def _previewAxialBounds(data, numbers):
        """Return mechanical helices plus union/common paired-dsDNA bounds."""
        paired_by_helix = {}
        for number in numbers:
            item = data.get(number)
            if item is None:
                continue
            scaffold_intervals = item.get('scaffold_intervals', ())
            staple_intervals = item.get('staple_intervals', ())
            if scaffold_intervals or staple_intervals:
                scaffold = set()
                staple = set()
                for low, high in scaffold_intervals:
                    scaffold.update(range(int(low), int(high)+1))
                for low, high in staple_intervals:
                    staple.update(range(int(low), int(high)+1))
                paired = scaffold & staple
            else:
                paired = set(range(int(item['low']), int(item['high'])+1))
            if paired:
                paired_by_helix[int(number)] = paired
        if not paired_by_helix:
            return [], 0, 0, 0, 0
        union = set().union(*paired_by_helix.values())
        common = set.intersection(*paired_by_helix.values())
        if not common:
            # No common full cross-section: use the union only as a safe
            # display interval and let the predictor report low confidence.
            common = set(union)
        # Use the longest continuous common run; disconnected paired islands
        # must not be treated as one calibrated axial span.
        runs = []
        for idx in sorted(common):
            if not runs or idx != runs[-1][-1]+1:
                runs.append([idx])
            else:
                runs[-1].append(idx)
        longest = max(runs, key=lambda values: (len(values), -values[0]))
        return (sorted(paired_by_helix), min(union), max(union),
                longest[0], longest[-1])

    def _baselinePreviewTransform(self):
        if self._simulationMeasurement:
            measurement = self._simulationMeasurement
            measured_numbers = list(measurement.get('helices', ()))
            mechanical, display_start, display_end, unused_start, unused_end = \
                self._previewAxialBounds(self._previewBaseData,
                                         measured_numbers)
            if not mechanical:
                mechanical = measured_numbers
                display_start = int(measurement.get('start', 0))
                display_end = int(measurement.get('end', display_start))
            measured_rate = measurement.get('twist_per_base_deg')
            if measured_rate is None:
                span = max(1, int(measurement.get('end', 0))-
                           int(measurement.get('start', 0)))
                measured_rate = float(
                    measurement.get('total_twist_deg', 0.0))/span
            return {
                'mode': 'baseline_twist',
                'helices': mechanical,
                'start': display_start,
                'end': display_end,
                'angle': float(measured_rate)*max(
                    0, display_end-display_start),
                'direction': 0.0,
                'source': 'simulation_measurement'}
        numbers = sorted(self._previewBaseData)
        if not numbers:
            return None
        mechanical, display_start, display_end, start, end = \
            self._previewAxialBounds(self._previewBaseData, numbers)
        if not mechanical or end <= start:
            return None
        region = {'mode': 'baseline_twist', 'helices': numbers,
                  'start': start, 'end': end}
        try:
            prediction = estimate_global_twist(
                region, self._previewBaseData,
                self._part.helicalPitch())
        except TwistBendError:
            return None
        # Calibrate the angular rate only where the full mechanical
        # cross-section exists, then extend that rate through the entire
        # displayed dsDNA length so the rods show one continuous global twist.
        return {
            'mode': 'baseline_twist', 'helices': mechanical,
            'start': display_start, 'end': display_end,
            'angle': prediction['twist_per_base_deg']*max(
                0, display_end-display_start),
            'direction': 0.0,
            'calibration_start': start, 'calibration_end': end}

    def _baselineBendTransform(self):
        measurement = self._simulationMeasurement
        if not measurement or abs(float(measurement.get(
                'bend_angle_deg', 0.0))) <= 1e-9:
            return None
        return {
            'mode': 'bend',
            'helices': list(measurement.get('helices', ())),
            'start': int(measurement.get('start', 0)),
            'end': int(measurement.get('end', 0)),
            'angle': float(measurement.get('bend_angle_deg', 0.0)),
            'direction': float(measurement.get('bend_direction_deg', 0.0)),
            'source': 'simulation_measurement'}

    @staticmethod
    def _applyPlanToData(working_data, plan):
        for edit in plan.get('edits', ()):
            number, idx = int(edit['helix']), int(edit['idx'])
            if number not in working_data:
                continue
            existing = working_data[number].setdefault('insertions', {})
            value = int(existing.get(idx, 0)) + int(edit['length'])
            if value:
                existing[idx] = value
            else:
                existing.pop(idx, None)
            working_data[number].setdefault('forbidden', set()).add(idx)

    def _predictionForTask(self, task, working_data):
        prediction = estimate_global_twist(
            task, working_data, self._part.helicalPitch())
        if self._measurementMatches(task.get('helices'), task.get('start'),
                                    task.get('end')):
            prediction = self._anchorTwistPrediction(
                prediction, prediction, self._simulationMeasurement)
        return prediction

    def _planForTask(self, task, working_data):
        if task['mode'] == 'remove_indels':
            edits = []
            for number in task['helices']:
                for idx, length in sorted(working_data[number].get(
                        'insertions', {}).items()):
                    if not task['start'] <= int(idx) <= task['end']:
                        continue
                    if ((int(length) > 0 and task.get('remove_insertions')) or
                            (int(length) < 0 and task.get('remove_deletions'))):
                        edits.append({'helix': int(number), 'idx': int(idx),
                                      'length': -int(length),
                                      'operation': 'remove_existing',
                                      'original_length': int(length)})
            plan = {'kind': 'remove_indels', 'edits': edits,
                    'removed_insertions': sum(
                        edit['original_length'] > 0 for edit in edits),
                    'removed_deletions': sum(
                        edit['original_length'] < 0 for edit in edits),
                    'preview_transform': None}
            self._applyPlanToData(working_data, plan)
            return plan
        model_before = estimate_global_twist(
            task, working_data, self._part.helicalPitch())
        before = self._predictionForTask(task, working_data)
        if task['mode'] == 'remove_twist':
            measured = task.get('measured_twist_per_base_deg')
            plan = plan_remove_twist(
                task, working_data, self._part.helicalPitch(), measured)
            if measured is not None:
                before = dict(before)
                before['model_twist_per_base_deg'] = before[
                    'twist_per_base_deg']
                before['twist_per_base_deg'] = float(measured)
                before['total_twist_deg'] = (
                    float(measured)*before['length_bp'])
                before['measurement_anchored'] = True
        elif task['mode'] == 'add_twist':
            plan = plan_add_twist(
                task, working_data, self._part.helicalPitch(), task['angle'],
                task['handedness'], task.get('indels_per_helix'))
        else:
            requested_angle = task['angle']
            requested_direction = task['direction']
            if self._measurementMatches(task.get('helices'),
                                        task.get('start'), task.get('end')):
                measured_angle = float(
                    self._simulationMeasurement.get('bend_angle_deg', 0.0))
                measured_direction = math.radians(float(
                    self._simulationMeasurement.get(
                        'bend_direction_deg', requested_direction)))
                target_direction = math.radians(float(requested_direction))
                dx = (requested_angle*math.cos(target_direction)-
                      measured_angle*math.cos(measured_direction))
                dy = (requested_angle*math.sin(target_direction)-
                      measured_angle*math.sin(measured_direction))
                requested_angle = math.hypot(dx, dy)
                requested_direction = (math.degrees(math.atan2(dy, dx)) %
                                       360.0 if requested_angle > 1e-9 else
                                       task['direction'])
            plan = plan_bend(
                task, working_data, requested_angle, requested_direction,
                material_twist_degrees=before['total_twist_deg'],
                incremental=True, elastic_compensation=True,
                lattice=('honeycomb' if self._part._step == 21
                         else 'square'))
            if self._measurementMatches(task.get('helices'),
                                        task.get('start'), task.get('end')):
                plan['measured_baseline_bend_angle_deg'] = measured_angle
                plan['measured_baseline_bend_direction_deg'] = (
                    math.degrees(measured_direction) % 360.0)
                plan['final_target_bend_angle_deg'] = task['angle']
                plan['final_target_bend_direction_deg'] = task['direction']
        self._applyPlanToData(working_data, plan)
        if task['mode'] in ('remove_twist', 'add_twist'):
            after = (plan.get('twist_prediction')
                     if task['mode'] == 'remove_twist' and
                     task.get('measured_twist_per_base_deg') is not None
                     else estimate_global_twist(
                         task, working_data, self._part.helicalPitch()))
            if (task['mode'] == 'add_twist' and
                    self._measurementMatches(task.get('helices'),
                                             task.get('start'),
                                             task.get('end'))):
                after = self._anchorTwistPrediction(
                    after, model_before, self._simulationMeasurement)
            # The preview begins with the current design's residual twist;
            # each ordered task contributes only its physical change.
            if task['mode'] == 'remove_twist':
                # Preview the absolute final material twist, rather than a
                # negative increment.  The latter can fail to cancel when the
                # baseline and selected regions have different boundaries.
                transform = {
                    'mode': 'remove_twist_final',
                    'helices': list(task['helices']),
                    'start': task['start'], 'end': task['end'],
                    'angle': after['total_twist_deg'],
                    'target_angle': after['total_twist_deg'],
                    'direction': 0.0}
            else:
                transform = {
                    'mode': task['mode'], 'helices': list(task['helices']),
                    'start': task['start'], 'end': task['end'],
                    'angle': (after['total_twist_deg']-
                              before['total_twist_deg']),
                    'direction': 0.0}
            plan['twist_prediction'] = after
            plan['twist_before_prediction'] = before
        else:
            predicted_angle = plan['elastic_prediction']['angle_degrees']
            transform = {
                'mode': 'bend', 'helices': list(task['helices']),
                'start': task['start'], 'end': task['end'],
                'angle': predicted_angle,
                'direction': plan.get('direction', task['direction'])}
            plan['twist_before_prediction'] = before
        plan['preview_transform'] = transform
        return plan

    @staticmethod
    def _previewTransformForPlan(plan):
        """Normalize legacy Remove Twist previews to absolute final state."""
        transform = deepcopy(plan.get('preview_transform'))
        if (transform and plan.get('kind') == 'remove_twist' and
                plan.get('twist_prediction')):
            final_angle = float(
                plan['twist_prediction'].get('total_twist_deg', 0.0))
            transform['mode'] = 'remove_twist_final'
            transform['angle'] = final_angle
            transform['target_angle'] = final_angle
        return transform

    def _displayPipeline(self, extra_transform=None):
        transforms = []
        baseline = self._baselinePreviewTransform()
        if baseline is not None:
            transforms.append(baseline)
        baseline_bend = self._baselineBendTransform()
        if baseline_bend is not None:
            transforms.append(baseline_bend)
        transforms.extend(deepcopy(self._previewTransforms))
        if extra_transform is not None:
            transforms.append(dict(extra_transform))
        self.preview.setPipeline(transforms)
        self.parameterPreview.setPipeline(transforms)

    def _updateBendSafety(self, selected, start, end):
        if not selected or end <= start:
            maximum = 360.0
        else:
            direction = math.radians(self.bendDirection.value())
            ux, uy = math.cos(direction), -math.sin(direction)
            projections = [self._data[number]['coord'][0] * ux +
                           self._data[number]['coord'][1] * uy
                           for number in selected]
            neutral = (max(projections) + min(projections)) / 2.0
            inner_extent = max(projections) - neutral
            lengths = [self._mechanicalLength(number, start, end)
                       for number in selected]
            length_nm = sum(lengths) / float(len(lengths)) * .34
            maximum = math.degrees(
                length_nm / max(.7, inner_extent + .7)) * .92
            maximum = max(1.0, min(720.0, maximum))
        current = min(self.bendAngle.value(), maximum)
        self.bendAngle.blockSignals(True)
        self.bendAngle.setMaximum(maximum)
        self.bendAngle.setValue(current)
        self.bendAngle.blockSignals(False)
        self.bendAngleSlider.blockSignals(True)
        self.bendAngleSlider.setMaximum(max(1, int(min(360, maximum))))
        self.bendAngleSlider.setValue(int(round(current)))
        self.bendAngleSlider.blockSignals(False)
        self.bendSafetyLabel.setText(
            '当前选区最大安全弯曲角约 %.1f°；超过此值会导致内侧半径反转。' %
            maximum)

    def _markParametersDirty(self):
        if self._applyingParameters:
            return
        # An Apply result is only a temporary parameter preview until the
        # user explicitly adds it as a task.  Every mode/selection/parameter
        # change must therefore remove that draft from both 3D canvases.
        # Previously the simulation-analysis early return and the ``wasDirty``
        # guard could leave the last mode's deformation in the pipeline.
        self._activeDraftTransform = None
        self.parameterPreview.setEdits([])
        # Always synchronize both canvases with the committed pipeline.  The
        # preview widget itself skips rebuilding when this pipeline is already
        # current, so slider motion stays inexpensive while a stale draft can
        # never survive because of an out-of-sync bookkeeping flag.
        self._displayPipeline()
        if self._currentMode() == 'simulation_analysis':
            self._parametersDirty = True
            self._appliedParameterTask = None
            self.applyParameters.setEnabled(False)
            self.addTask.setEnabled(False)
            self._refreshSimulationStatus()
            return
        self._parametersDirty = True
        self._appliedParameterTask = None
        self.applyParameters.setEnabled(True)
        self.addTask.setEnabled(False)
        self.applyParameterStatus.setText(
            '参数已更改。请点击 Apply Parameters 更新联动参数、物理预测和 3D 预览。')
        self.applyParameterStatus.setStyleSheet(
            'color:#8a5a00; background:#fff6dc; padding:7px; '
            'border:1px solid #ead39a;')
        live_summary = self.livePlanSummaries[self.channels.currentIndex()]
        live_summary.setText(
            '当前显示结果已过期；点击 Apply Parameters 后重新计算。')
        live_summary.setStyleSheet(
            'color:#6b7280; background:#f1f3f5; padding:8px; '
            'border:1px solid #d4d8dc;')
        # Repeated slider ticks remain inexpensive because setPipeline returns
        # immediately when the committed pipeline has not changed.

    def _controlsChanged(self, unused=None):
        """Update inexpensive selectors only; defer physics and 3D work."""
        start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        self.rangePicker.setRange(start, end)
        selected = self.helixPicker.selected()
        self.rangePicker.setSelected(selected)
        lengths = [self._mechanicalLength(number, start, end)
                   for number in selected]
        self.actualLabel.setText('实际长度：%s' %
            (', '.join('H%d=%d nt' % pair for pair in zip(selected, lengths))
             if lengths else '请选择 helix'))
        self.bendDirectionDial.setSelected(selected)
        self.bendDirectionDial.setValue(self.bendDirection.value())
        self._refreshSimulationStatus()
        self._refreshExistingIndels()
        self._markParametersDirty()

    def _resolveLinkedParameters(self):
        """Resolve the active channel's dependent values on explicit Apply."""
        mode = self._currentMode()
        if mode == 'add_twist':
            source = self._inputSource.get('add_twist', 'angle')
            if source == 'pitch':
                turns, actual, unused_current = self._twistPitchMetrics()
                total = (int(round(self.twistPitch.value()*turns-actual))
                         if turns > 0 else 0)
            elif source == 'indels':
                count = max(1, len(self.helixPicker.selected()))
                total = int(round(abs(self.twistIndels.value())*count))
                if self.twistIndels.value() < 0:
                    total *= -1
            else:
                sign = (1.0 if self.handedness.currentIndex() == 0
                        else -1.0)
                total = self._deltaForTargetTwist(
                    sign*self.twistAngle.value())
            self._syncTwistFromDelta(total)
        elif mode == 'bend':
            mean_offset = self._bendMeanOffset()
            self._syncingParameters = True
            if self._inputSource.get('bend', 'angle') == 'indels':
                angle = (math.degrees(self.bendIndels.value()*.34 /
                                      mean_offset)
                         if mean_offset > 1e-12 else 0.0)
                self.bendAngle.blockSignals(True)
                self.bendAngle.setValue(angle)
                self.bendAngle.blockSignals(False)
            else:
                strength = (math.radians(self.bendAngle.value()) *
                            mean_offset / .34
                            if mean_offset > 1e-12 else 0.0)
                self.bendIndels.blockSignals(True)
                self.bendIndels.setValue(strength)
                self.bendIndels.blockSignals(False)
            self._syncingParameters = False

    def _applyParameters(self):
        self._applyingParameters = True
        try:
            self._resolveLinkedParameters()
        finally:
            self._applyingParameters = False
        return self._calculateAndPreview()

    def _calculateAndPreview(self):
        """Run linked physics, indel planning and full 3D only on Apply."""
        start, end = sorted((self.startSpin.value(), self.endSpin.value()))
        self.rangePicker.setRange(start, end)
        selected = self.helixPicker.selected()
        self.rangePicker.setSelected(selected)
        mode = self._currentMode()
        if mode != 'remove_indels':
            selected_count_for_twist = max(1, len(selected))
            twist_total = int(round(abs(self.twistIndels.value()) *
                                    selected_count_for_twist))
            signed_total = (twist_total if self.twistIndels.value() >= 0
                            else -twist_total)
            self._setTwistPitchForDelta(signed_total)
            self._setTwistPrediction(
                self._twistEstimateForDelta(signed_total))
        self._updateBendSafety(selected, start, end)
        lengths = [self._mechanicalLength(number, start, end)
                   for number in selected]
        self.actualLabel.setText('实际长度：%s' %
            (', '.join('H%d=%d nt' % pair for pair in zip(selected, lengths))
             if lengths else '请选择 helix'))
        live_summary = self.livePlanSummaries[self.channels.currentIndex()]
        angle = self.bendAngle.value() if mode == 'bend' else self.twistAngle.value()
        if mode == 'remove_twist':
            angle = 0.0
        elif mode == 'add_twist' and self.handedness.currentIndex() == 1:
            angle *= -1
        self.twistDial.setValue(self.twistAngle.value())
        self.bendDirectionDial.setSelected(selected)
        self.bendDirectionDial.setValue(self.bendDirection.value())
        self.preview.setState(selected, start, end, mode, angle,
                              self.bendDirection.value())
        self.parameterPreview.setState(selected, start, end, mode, angle,
                                       self.bendDirection.value())
        self._displayPipeline()
        try:
            task = self._currentTask()
            live_working = deepcopy(self._previewWorkingData)
            live_plan = self._planForTask(task, live_working)
            self.parameterPreview.setEdits(live_plan.get('edits', []))
            self._activeDraftTransform = (
                deepcopy(live_plan.get('preview_transform')))
            self._displayPipeline(self._activeDraftTransform)
            insertions = sum(
                1 for edit in live_plan.get('edits', ())
                if edit.get('operation') != 'remove_existing' and
                edit['length'] > 0)
            deletions = sum(
                1 for edit in live_plan.get('edits', ())
                if edit.get('operation') != 'remove_existing' and
                edit['length'] < 0)
            helix_count = max(1, len(selected))
            if mode == 'remove_indels':
                detail = ('将移除选区内 %d 个 insertion 位点和 %d 个 '
                          'deletion/skip 位点；其他 indel 保持不变。') % (
                    live_plan['removed_insertions'],
                    live_plan['removed_deletions'])
            elif mode == 'remove_twist':
                prediction = live_plan['twist_prediction']
                final_preview_angle = float(
                    (live_plan.get('preview_transform') or {}).get(
                        'angle', 0.0))
                detail = (('优先撤销已有 indel：%d insertion 位点、%d '
                           'deletion/skip 位点；随后新增 %d insertion、%d '
                           'deletion。最终预览残余 Twist θ=%+.2f°；结果螺距 '
                           '%.4f base/turn；预测 %+.4f°/base，'
                           '总计 %+.2f°（%s）；J=%.1f nm⁴；连通 %.0f%%；'
                           '可信度：%s。%s') %
                          (live_plan.get('removed_insertions', 0),
                           live_plan.get('removed_deletions', 0),
                           live_plan.get('added_insertions', 0),
                           live_plan.get('added_deletions', 0),
                           final_preview_angle, live_plan['achieved_pitch'],
                           prediction['twist_per_base_deg'],
                           prediction['total_twist_deg'],
                           prediction['handedness'],
                           prediction['polar_moment_nm4'],
                           prediction['connectivity_fraction']*100.0,
                           prediction['confidence'], prediction['note']))
            elif mode == 'add_twist':
                signed_average = (insertions-deletions) / float(helix_count)
                turns, actual, current_pitch = self._twistPitchMetrics()
                achieved_pitch = ((actual + insertions - deletions) / turns
                                  if turns > 1e-12 else 0.0)
                prediction = live_plan['twist_prediction']
                angle = prediction['total_twist_deg']
                self.preview.setState(selected, start, end, mode, angle,
                                      self.bendDirection.value())
                self.parameterPreview.setState(
                    selected, start, end, mode, angle,
                    self.bendDirection.value())
                detail = ('当前/结果螺距 %.4f → %.4f base/turn；'
                          '整数分配后平均 %+.3f base/helix；预测 %+.4f°/base，'
                          '总计 %+.2f°（%s）；J=%.1f nm⁴；连通 %.0f%%；可信度：%s。%s' %
                          (current_pitch, achieved_pitch,
                           signed_average, prediction['twist_per_base_deg'],
                           prediction['total_twist_deg'],
                           prediction['handedness'],
                           prediction['polar_moment_nm4'],
                           prediction['connectivity_fraction']*100.0,
                           prediction['confidence'], prediction['note']))
            else:
                counts = list(live_plan['per_helix_counts'].values())
                mean_absolute = (sum(abs(value) for value in counts) /
                                 float(helix_count))
                radius = ('∞' if live_plan['radius_nm'] is None else
                          '%.2f nm' % live_plan['radius_nm'])
                prior_twist = live_plan['twist_before_prediction'][
                    'total_twist_deg']
                elastic = live_plan['elastic_prediction']
                detail = ('整数分配后平均 |indel| %.3f base/helix；目标 %.2f°；'
                          'SNUPI 校正 Bend 预测 %.2f°，半径 %s；内部几何预补偿 %.2f°；'
                          '已按此前材料 twist %+.2f° 重新计算内外侧分配；可信度：%s。%s' %
                          (mean_absolute, live_plan['requested_angle'],
                           elastic['angle_degrees'], radius,
                           live_plan['geometric_design_angle'], prior_twist,
                           elastic['confidence'], elastic['note']))
            if mode == 'remove_twist':
                live_summary.setText('当前参数预计：%s' % detail)
            else:
                live_summary.setText(
                    '当前参数预计：%d addition、%d deletion。%s' %
                    (insertions, deletions, detail))
            live_summary.setStyleSheet(
                'color:#284b63; background:#edf5fa; padding:8px; border:1px solid #cadde9;')
        except TwistBendError as error:
            self.parameterPreview.setEdits([])
            live_summary.setText('当前参数无法应用：%s' % error)
            live_summary.setStyleSheet(
                'color:#8b1e1e; background:#ffecec; padding:8px; border:1px solid #efb8b8;')
            self._parametersDirty = True
            self._appliedParameterTask = None
            self.applyParameters.setEnabled(True)
            self.addTask.setEnabled(False)
            self.applyParameterStatus.setText('计算失败：%s' % error)
            self.applyParameterStatus.setStyleSheet(
                'color:#8b1e1e; background:#ffecec; padding:7px; '
                'border:1px solid #efb8b8;')
            return False
        self._parametersDirty = False
        self._appliedParameterTask = deepcopy(self._currentTask())
        self.applyParameters.setEnabled(False)
        self.addTask.setEnabled(True)
        self.applyParameterStatus.setText(
            '计算与 3D 预览已更新。可以添加当前任务；再次修改参数后需要重新 Apply。')
        self.applyParameterStatus.setStyleSheet(
            'color:#234e34; background:#e9f7ef; padding:7px; '
            'border:1px solid #b9ddc5;')
        return True

    def _taskText(self, task):
        label = {'remove_indels': 'Remove Insertion/Deletion',
                 'remove_twist': 'Remove Twist', 'add_twist': 'Add Twist',
                 'bend': 'Bend'}[task['mode']]
        suffix = ''
        if task['mode'] == 'add_twist':
            suffix = (' · %s %.1f° · %+.3f base/helix · pitch %.4f' %
                      (task['handedness'], task['angle'],
                       task.get('indels_per_helix', 0.0),
                       task.get('target_pitch', self._part.helicalPitch())))
        elif task['mode'] == 'bend':
            suffix = (' · %.1f° @ %.1f° · |indel| %.3f/helix' %
                      (task['angle'], task['direction'],
                       task.get('indels_per_helix', 0.0)))
        return '%s · helix %s · base %d–%d%s' % (
            label, ','.join(str(v) for v in task['helices']),
            task['start'], task['end'], suffix)

    def _appendTask(self, task, rebuild=True):
        self._tasks.append(task)
        self.taskList.addItem(self._taskText(task))
        if rebuild:
            self._rebuildPlans()

    def _storeDraftTask(self):
        if self._parametersDirty or self._appliedParameterTask is None:
            QMessageBox.information(
                self, 'Twist and Bend',
                '当前参数尚未计算，或计算后又发生了修改。\n'
                '请先点击 Apply Parameters / 更新计算与预览。')
            return False
        task = deepcopy(self._appliedParameterTask)
        candidates = list(self._tasks)
        if self._draftTaskIndex is None:
            candidate_index = len(candidates)
            candidates.append(task)
        else:
            candidate_index = self._draftTaskIndex
            candidates[candidate_index] = task
        try:
            validate_regions(candidates)
        except TwistBendError as error:
            QMessageBox.warning(self, 'Twist and Bend', str(error))
            return False
        if self._draftTaskIndex is None:
            self._appendTask(task, rebuild=False)
            self._draftTaskIndex = candidate_index
        else:
            self._tasks[candidate_index] = task
            self.taskList.item(candidate_index).setText(self._taskText(task))
        self._rebuildPlans()
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(True)
        return True

    def _addTask(self):
        self._draftTaskIndex = None
        if self._storeDraftTask():
            # The next click creates another independent task, allowing a
            # Twist task and a Bend task to coexist in the same design.
            self._draftTaskIndex = None
            self._parametersDirty = True
            self._appliedParameterTask = None
            self.addTask.setEnabled(False)
            self.applyParameters.setEnabled(True)
            self.applyParameterStatus.setText(
                '当前任务已加入。若要继续添加或修改任务，请设置参数后重新 Apply。')
            self.applyParameterStatus.setStyleSheet(
                'color:#334e68; background:#eef5fb; padding:7px; '
                'border:1px solid #c8d8e6;')

    def _addAnotherTask(self):
        self._draftTaskIndex = None
        self._showStep(0)

    def _removeTask(self):
        row = self.taskList.currentRow()
        if row >= 0:
            self.taskList.takeItem(row)
            del self._tasks[row]
            if self._draftTaskIndex == row:
                self._draftTaskIndex = None
            elif self._draftTaskIndex is not None and row < self._draftTaskIndex:
                self._draftTaskIndex -= 1
            self._rebuildPlans()
            self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
                bool(self._tasks))

    def _rebuildPlans(self):
        plans, pending_plans, preview_transforms = [], [], []
        try:
            validate_regions(self._tasks)
            working_data = deepcopy(self._data)
            for task in self._tasks:
                if task.get('_applied'):
                    if task.get('_saved_plan'):
                        plan = dict(task['_saved_plan'])
                    else:
                        plan = self._planForTask(task, deepcopy(working_data))
                else:
                    plan = self._planForTask(task, working_data)
                plans.append(plan)
                preview_transform = self._previewTransformForPlan(plan)
                if preview_transform:
                    preview_transforms.append(preview_transform)
                if not task.get('_applied'):
                    pending_plans.append(plan)
            edits = merge_plans(plans)
        except TwistBendError as error:
            self._plans = []
            self._pendingPlans = []
            self.preview.setEdits([])
            self.summary.setText('无法应用：%s' % error)
            self.summary.setStyleSheet('color:#8b1e1e; background:#ffecec; padding:7px;')
            return False
        self._plans = plans
        self._pendingPlans = pending_plans
        self._previewTransforms = preview_transforms
        self._previewWorkingData = working_data
        self._activeDraftTransform = None
        self.preview.setEdits(edits)
        self._displayPipeline()
        insertions = sum(
            1 for edit in edits
            if edit.get('operation') != 'remove_existing' and
            edit['length'] > 0)
        deletions = sum(
            1 for edit in edits
            if edit.get('operation') != 'remove_existing' and
            edit['length'] < 0)
        removed_insertions = sum(
            1 for edit in edits
            if edit.get('operation') == 'remove_existing' and
            edit.get('original_length', 0) > 0)
        removed_deletions = sum(
            1 for edit in edits
            if edit.get('operation') == 'remove_existing' and
            edit.get('original_length', 0) < 0)
        details = []
        for plan in plans:
            if plan['kind'] == 'remove_indels':
                details.append(
                    'Remove Insertion/Deletion：移除 %d insertion、%d deletion/skip' %
                    (plan['removed_insertions'], plan['removed_deletions']))
            elif plan['kind'] == 'remove_twist':
                prediction = plan['twist_prediction']
                details.append(
                    'Remove Twist%s：基线 %+.4f°/base；%.4f base/turn，'
                    '撤销已有 %d insertion/%d deletion，新增 %d insertion/'
                    '%d deletion；残余 %+.4f°/base，%+.1f°（%s可信度）' %
                    (('（实测）' if plan.get('baseline_source') ==
                      'simulation_measurement' else ''),
                     plan.get('baseline_twist_per_base_deg', 0.0),
                     plan['achieved_pitch'],
                     plan.get('removed_insertions', 0),
                     plan.get('removed_deletions', 0),
                     plan.get('added_insertions', 0),
                     plan.get('added_deletions', 0),
                     prediction['twist_per_base_deg'],
                     prediction['total_twist_deg'],
                     prediction['confidence']))
            elif plan['kind'] == 'add_twist':
                prediction = plan['twist_prediction']
                details.append(
                    'Add Twist：%+.4f°/base，%+.1f°（%s，%s可信度）' %
                    (prediction['twist_per_base_deg'],
                     prediction['total_twist_deg'], prediction['handedness'],
                     prediction['confidence']))
            else:
                prior = plan.get('twist_before_prediction', {})
                elastic = plan.get('elastic_prediction', {})
                details.append(
                    'Bend SNUPI 校正预测 %.1f°，半径 %s（目标 %.1f°；'
                    '基于此前材料 twist %+.1f°；%s可信度）' %
                    (elastic.get('angle_degrees', plan['requested_angle']),
                     ('∞' if plan['radius_nm'] is None else
                      '%.2f nm' % plan['radius_nm']),
                     plan['requested_angle'],
                     prior.get('total_twist_deg', 0.0),
                     elastic.get('confidence', '低')))
        self.summary.setText(
            '%d 个任务；新增 %d insertion、%d deletion；撤销已有 %d '
            'insertion、%d deletion。%s' %
            (len(plans), insertions, deletions, removed_insertions,
             removed_deletions, '；'.join(details)))
        self.summary.setStyleSheet('color:#234e34; background:#e9f7ef; padding:7px;')
        return True

    def _accept(self):
        if not self._tasks:
            self._addTask()
            if not self._tasks:
                return
        if self._rebuildPlans():
            self.accept()

    def resultData(self):
        edits = merge_plans(self._pendingPlans)
        saved_tasks = []
        for task in self._tasks:
            saved = dict((key, value) for key, value in task.items()
                         if not key.startswith('_'))
            saved['_applied'] = True
            saved_tasks.append(saved)
        metadata = {'version': 4,
                    'prediction_calibration': SNUPI_CALIBRATION_VERSION,
                    'lattice': ('honeycomb' if self._part._step == 21
                                              else 'square'),
                    'tasks': saved_tasks, 'last_plans': self._plans,
                    'simulation_measurement': deepcopy(
                        self._simulationMeasurement)}
        return {'edits': edits, 'metadata': metadata}

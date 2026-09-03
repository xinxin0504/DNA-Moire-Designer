"""A dependency-free interactive 3D view rendered with QPainter.

The view deliberately avoids OpenGL.  cadnano 2.4 is commonly installed with
Qt combinations for which an OpenGL widget is unavailable; the painter based
renderer keeps the feature portable and lets us release all scene data as soon
as its dock is closed.
"""

import math

import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), ['QPoint', 'QTimer', 'Qt'])
util.qtWrapImport('QtGui', globals(), ['QBrush', 'QColor', 'QFont',
                                       'QPainter', 'QPen'])
util.qtWrapImport('QtWidgets', globals(), ['QComboBox', 'QHBoxLayout',
                                           'QLabel', 'QVBoxLayout', 'QWidget'])


_BASE_COLORS = {
    'A': QColor('#e74c3c'),
    'T': QColor('#2ecc71'),
    'C': QColor('#3498db'),
    'G': QColor('#f39c12'),
}
_UNKNOWN_COLOR = QColor('#aeb6bf')
_SCAF_OFFSET = 0.36
_STAP_OFFSET = -0.36
_RISE_PER_BASE = 0.34


class ThreeDCanvas(QWidget):
    """Rotatable helix/base renderer with model-object hit testing."""

    def __init__(self, document, partProvider, parent=None):
        super(ThreeDCanvas, self).__init__(parent)
        self._document = document
        self._partProvider = partProvider
        self._parts = []
        self._sceneSignature = None
        self._selectionSignature = None
        self._rods = []
        self._segments = []
        self._bases = []
        self._projectedBases = []
        self._activeBase = None
        self._yaw = -0.58
        self._pitch = -0.42
        self._zoom = 1.0
        self._fitScale = 8.0
        self._panX = 0.0
        self._panY = 0.0
        self._lastPos = QPoint()
        self._dragMode = None
        self._colorMode = 'sequence'
        self.setMinimumSize(320, 260)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setToolTip('Drag: rotate  |  Shift-drag/right-drag: pan  |  Wheel: zoom')

        self._timer = QTimer(self)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self.refreshFromModel)
        self._timer.start()
        self.refreshFromModel(force=True)

    def releaseResources(self):
        """Stop all background work and drop references/cached geometry."""
        self._timer.stop()
        self._rods = []
        self._segments = []
        self._bases = []
        self._projectedBases = []
        self._activeBase = None
        self._parts = []
        self._document = None
        self._partProvider = None

    def setColorMode(self, index):
        self._colorMode = 'sequence' if index == 0 else 'oligo'
        self.update()

    def resetCamera(self):
        self._yaw = -0.58
        self._pitch = -0.42
        self._zoom = 1.0
        self._panX = self._panY = 0.0
        self.update()

    def _modelSignature(self, parts):
        if not parts:
            return None
        values = []
        for part in parts:
            values.extend((id(part), part._step))
            for vh in sorted(part.getVirtualHelices(),
                             key=lambda item: item.number()):
                values.extend((vh.number(), vh.coord()))
                for strandSet in vh.getStrandSets():
                    for strand in strandSet:
                        values.extend((id(strand), strand.idxs(),
                                       strand.sequence(),
                                       strand.oligo().color()))
        return tuple(values)

    def _selectedState(self):
        if self._document is None:
            return (), set(), {}
        selected = set()
        endpointValues = {}
        signature = []
        for strandSet, strandDict in self._document.selectionDict().items():
            for strand, value in strandDict.items():
                selected.add(strand)
                endpointValues[strand] = value
                signature.append((id(strand), bool(value[0]), bool(value[1])))
        signature.sort()
        return tuple(signature), selected, endpointValues

    def refreshFromModel(self, force=False):
        if self._partProvider is None:
            return
        provided = self._partProvider()
        parts = (list(provided) if isinstance(provided, (list, tuple))
                 else ([provided] if provided is not None else []))
        signature = self._modelSignature(parts)
        if force or tuple(parts) != tuple(self._parts) or \
                signature != self._sceneSignature:
            self._parts = parts
            self._sceneSignature = signature
            self._buildScene()
        selectionSignature, unused, unused2 = self._selectedState()
        if selectionSignature != self._selectionSignature:
            self._selectionSignature = selectionSignature
            self.update()

    @staticmethod
    def _sequenceByBase(strand):
        """Map lattice indices to bases while accounting for skips/inserts."""
        sequence = strand.sequence() or ''
        if not strand.isDrawn5to3():
            sequence = sequence[::-1]
        insertions = dict((item.idx(), item.length())
                          for item in strand.insertionsOnStrand())
        result = {}
        cursor = 0
        low, high = strand.idxs()
        for idx in range(low, high + 1):
            insertionLength = insertions.get(idx, 0)
            if insertionLength < 0:
                result[idx] = ' '
                continue
            result[idx] = sequence[cursor:cursor + 1].upper() or ' '
            cursor += 1 + max(0, insertionLength)
        return result

    def _buildScene(self):
        self._rods = []
        self._segments = []
        self._bases = []
        self._projectedBases = []
        parts = self._parts
        if not parts:
            self.update()
            return

        rawRods = []
        rawSegments = []
        rawBases = []
        allPoints = []
        latticeOffsetX = 0.0
        for part in parts:
            twist = math.radians(part._twistPerBase)
            twistOffset = math.radians(part._twistOffset)
            prefix = 'H' if part._step == 21 else 'S'
            virtualHelices = sorted(part.getVirtualHelices(),
                                    key=lambda item: item.number())
            localAxes = [part.latticeCoordToPositionXY(*vh.coord())
                         for vh in virtualHelices]
            minLocalX = min((point[0] for point in localAxes), default=0.0)
            minLocalY = min((point[1] for point in localAxes), default=0.0)
            maxLocalX = max((point[0] for point in localAxes), default=0.0)
            partWidth = max(maxLocalX - minLocalX, part.radius() * 2.0)
            for vh in virtualHelices:
                row, col = vh.coord()
                localX, localY = part.latticeCoordToPositionXY(row, col)
                axisX = localX - minLocalX + latticeOffsetX
                axisY = localY - minLocalY
                strandSets = vh.getStrandSets()
                strands = [strand for strandSet in strandSets
                           for strand in strandSet]
                if strands:
                    low = min(strand.lowIdx() for strand in strands)
                    high = max(strand.highIdx() for strand in strands)
                else:
                    low, high = 0, max(part.stepSize() - 1, 1)
                p0 = (axisX, axisY, low * _RISE_PER_BASE)
                p1 = (axisX, axisY, high * _RISE_PER_BASE)
                rawRods.append((p0, p1, '%s%d' % (prefix, vh.number())))
                allPoints.extend((p0, p1))

                for strandSet in strandSets:
                    isScaffold = strandSet.isScaffold()
                    radialSign = (_SCAF_OFFSET if isScaffold
                                  else _STAP_OFFSET)
                    strandType = 'scaffold' if isScaffold else 'staple'
                    for strand in strandSet:
                        seqByBase = self._sequenceByBase(strand)
                        oligoColor = QColor(strand.oligo().color())
                        strandPoints = []
                        for idx in range(strand.lowIdx(),
                                         strand.highIdx() + 1):
                            angle = twistOffset + idx * twist
                            radius = radialSign
                            point = (axisX + radius * math.cos(angle),
                                     axisY + radius * math.sin(angle),
                                     idx * _RISE_PER_BASE)
                            strandPoints.append(point)
                            rawBases.append((point, strand, idx,
                                             seqByBase[idx], oligoColor,
                                             '%s%d' % (prefix, vh.number()),
                                             strandType))
                            allPoints.append(point)
                        if strandPoints:
                            rawSegments.append((strandPoints[0],
                                                strandPoints[-1],
                                                oligoColor, strand))
            latticeOffsetX += partWidth + max(part.radius() * 6.0, 6.0)

        basePoints = dict(((strand, idx), point)
                          for point, strand, idx, unused_base,
                          unused_color, unused_number, unused_kind
                          in rawBases)
        for part in parts:
            for vh in part.getVirtualHelices():
                for strandSet in vh.getStrandSets():
                    for strand in strandSet:
                        connected = strand.connection3p()
                        if connected is None or connected.part() is part:
                            continue
                        point5 = basePoints.get((strand,
                                                strand.idx3Prime()))
                        point3 = basePoints.get((connected,
                                                connected.idx5Prime()))
                        if point5 is not None and point3 is not None:
                            rawSegments.append((point5, point3,
                                                QColor(strand.oligo().color()),
                                                strand))

        if allPoints:
            minX = min(point[0] for point in allPoints)
            maxX = max(point[0] for point in allPoints)
            minY = min(point[1] for point in allPoints)
            maxY = max(point[1] for point in allPoints)
            minZ = min(point[2] for point in allPoints)
            maxZ = max(point[2] for point in allPoints)
            center = ((minX + maxX) / 2.0, (minY + maxY) / 2.0,
                      (minZ + maxZ) / 2.0)
            extent = max(maxX - minX, maxY - minY, maxZ - minZ, 1.0)
            self._fitScale = 0.68 * min(max(self.width(), 320),
                                        max(self.height(), 260)) / extent
        else:
            center = (0.0, 0.0, 0.0)
            self._fitScale = 8.0

        def centered(point):
            return (point[0] - center[0], point[1] - center[1],
                    point[2] - center[2])
        self._rods = [(centered(a), centered(b), number)
                      for a, b, number in rawRods]
        self._segments = [(centered(a), centered(b), color, strand)
                          for a, b, color, strand in rawSegments]
        self._bases = [(centered(point), strand, idx, base, color, number, kind)
                       for point, strand, idx, base, color, number, kind in rawBases]
        self.update()

    def _rotate(self, point):
        x, y, z = point
        cy, sy = math.cos(self._yaw), math.sin(self._yaw)
        x, z = x * cy + z * sy, -x * sy + z * cy
        cp, sp = math.cos(self._pitch), math.sin(self._pitch)
        y, z = y * cp - z * sp, y * sp + z * cp
        return x, y, z

    def _project(self, point):
        x, y, depth = self._rotate(point)
        scale = self._fitScale * self._zoom
        return (self.width() / 2.0 + self._panX + x * scale,
                self.height() / 2.0 + self._panY - y * scale, depth)

    @staticmethod
    def _lighter(color, factor=125):
        return color.lighter(factor)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor('#182028'))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing,
                              self._dragMode is None)
        if not self._parts:
            painter.setPen(QColor('#d5d8dc'))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             'No active design')
            return

        selectedSignature, selectedStrands, endpointValues = self._selectedState()
        rodPen = QPen(QColor('#566573'))
        rodPen.setWidthF(8.0)
        rodPen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(rodPen)
        projectedRods = []
        for p0, p1, number in self._rods:
            a, b = self._project(p0), self._project(p1)
            projectedRods.append((a, b, number))
            painter.drawLine(int(a[0]), int(a[1]), int(b[0]), int(b[1]))

        # Exact 2D oligo colors are retained on the strand rods in both modes.
        for p0, p1, color, strand in self._segments:
            a, b = self._project(p0), self._project(p1)
            pen = QPen(color)
            pen.setWidthF(3.2 if strand not in selectedStrands else 5.0)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.drawLine(int(a[0]), int(a[1]), int(b[0]), int(b[1]))

        projected = []
        for item in self._bases:
            point, strand, idx, base, oligoColor, number, kind = item
            projected.append((self._project(point), item))
        projected.sort(key=lambda entry: entry[0][2])
        self._projectedBases = []
        radius = max(2.4, min(6.5, 3.2 * math.sqrt(self._zoom)))
        for projection, item in projected:
            x, y, depth = projection
            point, strand, idx, base, oligoColor, number, kind = item
            color = (_BASE_COLORS.get(base, _UNKNOWN_COLOR)
                     if self._colorMode == 'sequence' else oligoColor)
            isActive = self._activeBase == (strand, idx)
            endpoint = endpointValues.get(strand, (False, False))
            isEndpointSelected = ((idx == strand.lowIdx() and endpoint[0]) or
                                  (idx == strand.highIdx() and endpoint[1]))
            isSelected = strand in selectedStrands and (endpoint == (True, True)
                                                         or isEndpointSelected)
            if isActive or isSelected:
                painter.setPen(QPen(QColor('#fff176'), 2.4))
                painter.setBrush(QBrush(color))
                painter.drawEllipse(int(x - radius - 2), int(y - radius - 2),
                                    int(2 * radius + 4), int(2 * radius + 4))
            else:
                painter.setPen(QPen(self._lighter(color), 0.8))
                painter.setBrush(QBrush(color))
                painter.drawEllipse(int(x - radius), int(y - radius),
                                    int(2 * radius), int(2 * radius))
            self._projectedBases.append((x, y, radius + 4, strand, idx,
                                         number, base, kind, depth))

        labelFont = QFont()
        labelFont.setPointSize(8)
        painter.setFont(labelFont)
        painter.setPen(QColor('#ecf0f1'))
        for a, b, number in projectedRods:
            label = a if a[1] < b[1] else b
            painter.drawText(int(label[0] + 6), int(label[1] - 4), str(number))

        painter.setPen(QColor('#aab7b8'))
        legend = ('A red  T green  C blue  G orange' if self._colorMode == 'sequence'
                  else 'Colors match 2D oligos')
        painter.drawText(10, self.height() - 10, legend)

    def _baseAt(self, pos):
        best = None
        bestDistance = None
        # reverse order favors the visually nearest sphere
        for x, y, radius, strand, idx, number, base, kind, depth in reversed(self._projectedBases):
            distance = (pos.x() - x) ** 2 + (pos.y() - y) ** 2
            if distance <= radius ** 2 and (bestDistance is None or distance < bestDistance):
                bestDistance = distance
                best = (strand, idx, number, base, kind)
        return best

    def _selectInDocument(self, strand, idx, wholeOligo=False):
        doc = self._document
        if doc is None:
            return
        # Deselect through model signals so the existing 2D items update.
        for strandSet, strandDict in list(doc.selectionDict().items()):
            for oldStrand in list(strandDict.keys()):
                doc.removeStrandFromSelection(oldStrand)
        doc.updateSelection()
        if wholeOligo:
            strand5p = strand.oligo().strand5p()
            targets = (list(strand5p.generator3pStrand())
                       if strand5p is not None else [strand])
        else:
            targets = [strand]
        for target in targets:
            doc.addStrandToSelection(target, (True, True))
        doc.updateSelection()
        self._activeBase = (strand, idx)
        self._selectionSignature = None
        self.update()

    def mousePressEvent(self, event):
        self._lastPos = event.pos()
        if event.button() == Qt.MouseButton.LeftButton:
            hit = self._baseAt(event.pos())
            if hit is not None:
                strand, idx, number, base, kind = hit
                self._selectInDocument(strand, idx, wholeOligo=False)
                self.setToolTip('Helix %s, base %s, %s, %s' %
                                (number, idx, kind, base.strip() or 'unassigned'))
                return
            self._dragMode = ('pan' if event.modifiers() & Qt.KeyboardModifier.ShiftModifier
                              else 'rotate')
        elif event.button() in (Qt.MouseButton.RightButton,
                                Qt.MouseButton.MiddleButton):
            self._dragMode = 'pan'

    def mouseDoubleClickEvent(self, event):
        hit = self._baseAt(event.pos())
        if hit is not None:
            strand, idx, number, base, kind = hit
            self._selectInDocument(strand, idx, wholeOligo=True)

    def mouseMoveEvent(self, event):
        if self._dragMode is None:
            hit = self._baseAt(event.pos())
            if hit is not None:
                strand, idx, number, base, kind = hit
                self.setToolTip('Helix %s, base %s, %s, %s' %
                                (number, idx, kind, base.strip() or 'unassigned'))
            return
        delta = event.pos() - self._lastPos
        self._lastPos = event.pos()
        if self._dragMode == 'rotate':
            self._yaw += delta.x() * 0.009
            self._pitch = max(-1.52, min(1.52, self._pitch + delta.y() * 0.009))
        else:
            self._panX += delta.x()
            self._panY += delta.y()
        self.update()

    def mouseReleaseEvent(self, event):
        self._dragMode = None
        self.update()

    def wheelEvent(self, event):
        factor = math.pow(1.0015, event.angleDelta().y())
        self._zoom = max(0.12, min(18.0, self._zoom * factor))
        self.update()


class ThreeDView(QWidget):
    """Small control bar plus the 3D canvas used inside a dock widget."""

    def __init__(self, document, partProvider, parent=None):
        super(ThreeDView, self).__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        controls = QHBoxLayout()
        controls.addWidget(QLabel('Color:'))
        self.colorMode = QComboBox()
        self.colorMode.addItem('Sequence colors')
        self.colorMode.addItem('2D oligo colors')
        controls.addWidget(self.colorMode)
        controls.addStretch(1)
        hint = QLabel('Drag rotate · Shift/right drag pan · Wheel zoom')
        hint.setStyleSheet('color: #6c757d;')
        controls.addWidget(hint)
        layout.addLayout(controls)
        self.canvas = ThreeDCanvas(document, partProvider, self)
        layout.addWidget(self.canvas, 1)
        self.colorMode.currentIndexChanged.connect(self.canvas.setColorMode)

    def releaseResources(self):
        self.canvas.releaseResources()

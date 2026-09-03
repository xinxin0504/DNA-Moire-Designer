#!/usr/bin/env python
# encoding: utf-8

"""Shared, scalable artwork for connections between Hybrid path panels."""

from cadnano2.views import styles
import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), ['QPoint', 'QPointF', 'QTimer', 'Qt'])
util.qtWrapImport('QtGui', globals(), ['QColor', 'QFont', 'QFontMetricsF',
                                       'QPainter', 'QPainterPath', 'QPen'])
util.qtWrapImport('QtWidgets', globals(), ['QWidget'])


def hybridConnections(document):
    """Yield every physical 3'-to-5' cross-part connection exactly once."""
    for part in document.parts():
        for virtualHelix in part.getVirtualHelices():
            for strandSet in (virtualHelix.scaffoldStrandSet(),
                              virtualHelix.stapleStrandSet()):
                for strand3Prime in strandSet:
                    strand5Prime = strand3Prime.connection3p()
                    if strand5Prime is not None and \
                            strand5Prime.part() is not part:
                        yield strand3Prime, strand5Prime


def hybridEndpointLabel(strand, end):
    """Return the stable lattice-qualified label for one physical end."""
    prefix = 'H' if strand.part()._step == 21 else 'S'
    index = strand.idx3Prime() if end == '3p' else strand.idx5Prime()
    return '%s: %d[%d]' % (prefix, strand.virtualHelix().number(), index)


def hybridSceneEndpoint(root, strand, end):
    """Map a model endpoint to the scene containing ``root``."""
    try:
        partItem = root.partItemForPart(strand.part())
        virtualHelixItem = partItem.itemForVirtualHelix(
                                            strand.virtualHelix())
    except (KeyError, AttributeError):
        return None
    index = strand.idx3Prime() if end == '3p' else strand.idx5Prime()
    point = virtualHelixItem.upperLeftCornerOfBaseType(
                                    index, strand.strandType())
    return virtualHelixItem.mapToScene(
        point[0] + 0.5 * styles.PATH_BASE_WIDTH,
        point[1] + 0.5 * styles.PATH_BASE_WIDTH)


def hybridCurvePath(threePrimePoint, fivePrimePoint, scale=1.0):
    """Construct a smooth cubic whose endpoints remain base-anchored."""
    start = QPointF(threePrimePoint)
    end = QPointF(fivePrimePoint)
    verticalDistance = abs(end.y() - start.y())
    bend = max(24.0 * scale, verticalDistance * 0.38)
    direction = 1.0 if end.y() >= start.y() else -1.0
    firstControl = start + QPointF(0.0, direction * bend)
    secondControl = end - QPointF(0.0, direction * bend)
    path = QPainterPath(start)
    path.cubicTo(firstControl, secondControl, end)
    return path


def drawHybridConnection(painter, threePrimePoint, fivePrimePoint,
                         strand3Prime, selected=False,
                         drawLabels=True, scale=1.0):
    """Draw one curve plus its short, rounded 3' direction accent."""
    color = QColor(styles.selected_color if selected
                   else strand3Prime.oligo().color())
    baseWidth = (styles.PATH_SCAFFOLD_STROKE_WIDTH
                 if strand3Prime.isScaffold()
                 else styles.PATH_STRAND_STROKE_WIDTH)
    mainWidth = max(1.0, float(baseWidth) * scale)
    path = hybridCurvePath(threePrimePoint, fivePrimePoint, scale)

    pen = QPen(color, mainWidth)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    painter.drawPath(path)

    # The final visual cue is deliberately short: the whole curve remains
    # uniformly readable, while only the physical 3' end receives emphasis.
    length = max(1.0, path.length())
    accentFraction = min(0.16, (0.60 * styles.PATH_BASE_WIDTH * scale) /
                                      length)
    accentEnd = path.pointAtPercent(accentFraction)
    accent = QPainterPath(QPointF(threePrimePoint))
    accent.lineTo(accentEnd)
    accentPen = QPen(color, max(mainWidth + 1.0 * scale,
                               mainWidth * 1.8))
    accentPen.setCapStyle(Qt.PenCapStyle.RoundCap)
    painter.setPen(accentPen)
    painter.drawPath(accent)

    if not drawLabels:
        return
    font = QFont('Monaco')
    font.setPointSizeF(max(4.0, 8.0 * scale))
    painter.setFont(font)
    painter.setPen(QPen(color))
    metrics = QFontMetricsF(font)
    labels = ((hybridEndpointLabel(strand3Prime.connection3p(), '5p'),
               QPointF(threePrimePoint), True),
              (hybridEndpointLabel(strand3Prime, '3p'),
               QPointF(fivePrimePoint), False))
    gap = max(5.0, 0.35 * styles.PATH_BASE_WIDTH * scale)
    for text, point, placeLeft in labels:
        bounds = metrics.boundingRect(text)
        x = (point.x() - bounds.width() - gap if placeLeft
             else point.x() + gap)
        y = point.y() - gap
        painter.drawText(QPointF(x, y), text)


class HybridConnectionOverlay(QWidget):
    """Transparent artwork spanning the two independent Hybrid Path Views."""

    def __init__(self, parent, geometryWidget, document,
                 honeyView, honeyRoot, squareView, squareRoot):
        super(HybridConnectionOverlay, self).__init__(parent)
        self._geometryWidget = geometryWidget
        self._document = document
        self._views = {21: (honeyView, honeyRoot),
                       32: (squareView, squareRoot)}
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents,
                          True)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.setAutoFillBackground(False)
        self.hide()
        self._timer = QTimer(self)
        self._timer.setInterval(80)
        self._timer.timeout.connect(self.update)
        document.documentHybridConnectionChangedSignal.connect(
                                            lambda unused: self.update())

    def setActive(self, active):
        if active:
            self._syncGeometry()
            self.raise_()
            self.show()
            self._timer.start()
            self.update()
        else:
            self._timer.stop()
            self.hide()

    def _syncGeometry(self):
        """Cover the Path splitter without becoming one of its panes."""
        target = self._geometryWidget
        # QWidget.mapFrom(widget, point) is only safe for the documented
        # parent relationship.  The overlay and splitter are siblings on
        # macOS, where passing the splitter directly can segfault inside Qt.
        globalTopLeft = target.mapToGlobal(QPoint(0, 0))
        topLeft = self.parentWidget().mapFromGlobal(globalTopLeft)
        self.setGeometry(topLeft.x(), topLeft.y(),
                         target.width(), target.height())

    def _widgetEndpoint(self, strand, end):
        view, root = self._views.get(strand.part()._step, (None, None))
        if view is None:
            return None
        scenePoint = hybridSceneEndpoint(root, strand, end)
        if scenePoint is None:
            return None
        viewportPoint = view.mapFromScene(scenePoint)
        globalPoint = view.viewport().mapToGlobal(viewportPoint)
        widgetPoint = self.mapFromGlobal(globalPoint)
        return QPointF(widgetPoint)

    def paintEvent(self, event):
        self._syncGeometry()
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        viewScales = [abs(view.transform().m11())
                      for view, unusedRoot in self._views.values()]
        displayScale = (sum(viewScales) / len(viewScales)
                        if viewScales else 1.0)
        selectedObjects = getattr(self._document, '_selectionDict', {})
        for strand3Prime, strand5Prime in hybridConnections(self._document):
            point3 = self._widgetEndpoint(strand3Prime, '3p')
            point5 = self._widgetEndpoint(strand5Prime, '5p')
            if point3 is None or point5 is None:
                continue
            selected = (strand3Prime in selectedObjects or
                        strand5Prime in selectedObjects or
                        strand3Prime.oligo() in selectedObjects)
            drawHybridConnection(painter, point3, point5, strand3Prime,
                                 selected=selected, drawLabels=True,
                                 scale=displayScale)
        painter.end()

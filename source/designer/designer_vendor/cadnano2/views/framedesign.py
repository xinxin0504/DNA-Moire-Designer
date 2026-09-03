"""Frame Design dialog: straight edges with vertex-local curvature."""

import math
import os

import cadnano2.util as util
from ..model.io.frame import (FRAME_PRESETS, plan_frame, polygon_from_spec,
                              polygon_metrics, safe_name)
from ..model.io.curved import CURVED_SCAFFOLD_MAX_BASES

util.qtWrapImport('QtCore', globals(), ['QDir', 'QPointF', 'QRectF', 'Qt'])
util.qtWrapImport('QtGui', globals(), [
    'QBrush', 'QColor', 'QImage', 'QPainter', 'QPainterPath', 'QPen'])
util.qtWrapImport('QtWidgets', globals(), [
    'QComboBox', 'QDialog', 'QDialogButtonBox', 'QDoubleSpinBox',
    'QFileDialog', 'QFormLayout', 'QHBoxLayout', 'QLabel', 'QLineEdit',
    'QMessageBox', 'QPushButton', 'QScrollArea', 'QSpinBox', 'QVBoxLayout',
    'QWidget'])


class FramePreview(QWidget):
    def __init__(self, parent=None):
        super(FramePreview, self).__init__(parent)
        self._plan = None
        self.setMinimumSize(590, 410)

    def setPlan(self, plan):
        self._plan = plan
        self.update()

    def paintEvent(self, unused_event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor(248, 249, 251))
        if not self._plan:
            return
        vertices = self._plan['vertices_nm']
        minimum_x = min(point[0] for point in vertices)
        maximum_x = max(point[0] for point in vertices)
        minimum_y = min(point[1] for point in vertices)
        maximum_y = max(point[1] for point in vertices)
        width = max(1e-6, maximum_x - minimum_x)
        height = max(1e-6, maximum_y - minimum_y)
        # Reserve explicit annotation gutters around the geometry.  Edge
        # dimensions sit close to the polygon, while overall W/H dimensions
        # use the outer gutter; neither should overlap vertex-angle labels.
        bounds = self.rect().adjusted(92, 58, -112, -142)
        scale = min(bounds.width() / width, bounds.height() / height)
        centre_x, centre_y = bounds.center().x(), bounds.center().y()

        def mapped(point):
            return QPointF(
                centre_x + (point[0] - 0.5 * (minimum_x + maximum_x)) * scale,
                centre_y - (point[1] - 0.5 * (minimum_y + maximum_y)) * scale)

        painter.setPen(QPen(QColor(112, 124, 138), 10.0,
                            Qt.PenStyle.SolidLine,
                            Qt.PenCapStyle.RoundCap,
                            Qt.PenJoinStyle.RoundJoin))
        path = QPainterPath(mapped(vertices[0]))
        for point in vertices[1:]:
            path.lineTo(mapped(point))
        path.closeSubpath()
        painter.drawPath(path)

        # Orange vertex zones make the localised-curvature design explicit;
        # the unmarked grey spans are the zero-curvature frame edges.
        painter.setPen(QPen(QColor(226, 119, 54), 12.0,
                            Qt.PenStyle.SolidLine,
                            Qt.PenCapStyle.RoundCap))
        for index, vertex in enumerate(vertices):
            previous = vertices[index - 1]
            following = vertices[(index + 1) % len(vertices)]
            trim = min(self._plan['tangent_trim_nm'][index],
                       0.35 * min(math.hypot(vertex[0]-previous[0],
                                            vertex[1]-previous[1]),
                                  math.hypot(following[0]-vertex[0],
                                             following[1]-vertex[1])))
            first_length = math.hypot(vertex[0]-previous[0],
                                      vertex[1]-previous[1])
            second_length = math.hypot(following[0]-vertex[0],
                                       following[1]-vertex[1])
            first = (vertex[0] + (previous[0]-vertex[0]) * trim/first_length,
                     vertex[1] + (previous[1]-vertex[1]) * trim/first_length)
            second = (vertex[0] + (following[0]-vertex[0]) * trim/second_length,
                      vertex[1] + (following[1]-vertex[1]) * trim/second_length)
            local = QPainterPath(mapped(first))
            local.quadTo(mapped(vertex), mapped(second))
            painter.drawPath(local)

        painter.setPen(QPen(QColor(35, 53, 70), 1.0))
        font = painter.font()
        font.setPointSize(9)
        painter.setFont(font)

        def boxed_text(centre, text, color, width_px=116.0,
                       height_px=34.0):
            centre_x = max(width_px / 2.0 + 5.0,
                           min(self.width() - width_px / 2.0 - 5.0,
                               centre.x()))
            centre_y = max(height_px / 2.0 + 5.0,
                           min(self.height() - height_px / 2.0 - 5.0,
                               centre.y()))
            rect = QRectF(centre_x - width_px / 2.0,
                          centre_y - height_px / 2.0,
                          width_px, height_px)
            painter.fillRect(rect, QColor(248, 249, 251, 232))
            painter.setPen(QPen(color, 1.0))
            painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, text)

        def draw_tick(point, normal, color, half_length=5.0):
            painter.setPen(QPen(color, 1.15))
            painter.drawLine(
                QPointF(point.x() - normal[0] * half_length,
                        point.y() - normal[1] * half_length),
                QPointF(point.x() + normal[0] * half_length,
                        point.y() + normal[1] * half_length))

        def edge_dimension(first, second, outward, label):
            """Draw a conventional dimension line outside one frame side."""
            offset = 20.0
            first_dim = QPointF(first.x() + outward[0] * offset,
                                first.y() + outward[1] * offset)
            second_dim = QPointF(second.x() + outward[0] * offset,
                                 second.y() + outward[1] * offset)
            color = QColor(43, 103, 145)
            painter.setPen(QPen(color, 1.15))
            painter.drawLine(first, QPointF(
                first.x() + outward[0] * (offset + 5.0),
                first.y() + outward[1] * (offset + 5.0)))
            painter.drawLine(second, QPointF(
                second.x() + outward[0] * (offset + 5.0),
                second.y() + outward[1] * (offset + 5.0)))
            painter.drawLine(first_dim, second_dim)
            draw_tick(first_dim, outward, color)
            draw_tick(second_dim, outward, color)
            middle = QPointF(
                0.5 * (first_dim.x() + second_dim.x()) +
                outward[0] * 11.0,
                0.5 * (first_dim.y() + second_dim.y()) +
                outward[1] * 11.0)
            boxed_text(middle, label, color, 158.0, 20.0)

        centroid = (
            sum(point[0] for point in vertices) / float(len(vertices)),
            sum(point[1] for point in vertices) / float(len(vertices)))
        for index, vertex in enumerate(vertices):
            point = mapped(vertex)
            inward_x = centroid[0] - vertex[0]
            inward_y = centroid[1] - vertex[1]
            inward_length = max(1e-9, math.hypot(inward_x, inward_y))
            # Convert the physical inward vector to screen coordinates.
            inward_screen = (inward_x / inward_length,
                              -inward_y / inward_length)
            label_centre = QPointF(
                point.x() + inward_screen[0] * 43.0,
                point.y() + inward_screen[1] * 43.0)
            painter.setPen(QPen(QColor(183, 84, 24), 1.0))
            painter.drawLine(
                QPointF(point.x() + inward_screen[0] * 9.0,
                        point.y() + inward_screen[1] * 9.0),
                QPointF(label_centre.x() - inward_screen[0] * 18.0,
                        label_centre.y() - inward_screen[1] * 18.0))
            turn = float(self._plan['turn_angles_degrees'][index])
            interior = float(self._plan.get(
                'interior_angles_degrees', [180.0 - turn] * len(vertices))[index])
            boxed_text(label_centre,
                       'V%d  弯曲转角 %.1f°\n(内角 %.1f°)' %
                       (index + 1, turn, interior),
                       QColor(160, 69, 20), 142.0, 36.0)

        # Label every neutral-axis edge with its realised physical length.
        # The polygon may be scaled slightly when the closed centreline is
        # snapped to a complete lattice period, so these are intentionally
        # the realised values rather than merely repeating the input fields.
        for index, (first, second) in enumerate(zip(
                vertices, vertices[1:] + vertices[:1])):
            first_screen, second_screen = mapped(first), mapped(second)
            dx, dy = second[0] - first[0], second[1] - first[1]
            length = max(1e-9, math.hypot(dx, dy))
            outward = (dy / length, -dx / length)
            # Physical y is inverted by mapped(), hence the screen-space
            # outward normal has the opposite y component.
            outward_screen = (outward[0], -outward[1])
            edge_dimension(
                first_screen, second_screen, outward_screen,
                '边%d（顶点间） %.1f nm' %
                (index + 1, self._plan['side_lengths_nm'][index]))

        # Overall sharp-vertex envelope.  Separate W and H dimension lines
        # make it unambiguous that these are not individual side lengths.
        screen_vertices = [mapped(point) for point in vertices]
        left_x = min(point.x() for point in screen_vertices)
        right_x = max(point.x() for point in screen_vertices)
        top_y = min(point.y() for point in screen_vertices)
        bottom_y = max(point.y() for point in screen_vertices)
        overall_color = QColor(0, 123, 132)
        overall_y = bottom_y + 53.0
        painter.setPen(QPen(overall_color, 1.35))
        painter.drawLine(QPointF(left_x, bottom_y + 5.0),
                         QPointF(left_x, overall_y + 5.0))
        painter.drawLine(QPointF(right_x, bottom_y + 5.0),
                         QPointF(right_x, overall_y + 5.0))
        painter.drawLine(QPointF(left_x, overall_y),
                         QPointF(right_x, overall_y))
        draw_tick(QPointF(left_x, overall_y), (0.0, 1.0), overall_color)
        draw_tick(QPointF(right_x, overall_y), (0.0, 1.0), overall_color)
        boxed_text(QPointF(0.5 * (left_x + right_x), overall_y + 12.0),
                   '总宽 W（顶点外包络）= %.1f nm' % width,
                   overall_color, 178.0, 20.0)

        overall_x = right_x + 55.0
        painter.setPen(QPen(overall_color, 1.35))
        painter.drawLine(QPointF(right_x + 5.0, top_y),
                         QPointF(overall_x + 5.0, top_y))
        painter.drawLine(QPointF(right_x + 5.0, bottom_y),
                         QPointF(overall_x + 5.0, bottom_y))
        painter.drawLine(QPointF(overall_x, top_y),
                         QPointF(overall_x, bottom_y))
        draw_tick(QPointF(overall_x, top_y), (1.0, 0.0), overall_color)
        draw_tick(QPointF(overall_x, bottom_y), (1.0, 0.0), overall_color)
        painter.save()
        painter.translate(overall_x + 13.0, 0.5 * (top_y + bottom_y))
        painter.rotate(-90.0)
        boxed_text(QPointF(0.0, 0.0),
                   '总高 H（顶点外包络）= %.1f nm' % height,
                   overall_color, 178.0, 20.0)
        painter.restore()

        legend_font = painter.font()
        legend_font.setPointSize(8)
        painter.setFont(legend_font)
        painter.setPen(QPen(QColor(22, 91, 135), 1.2))
        painter.drawText(
            QRectF(14, self.height() - 45, self.width() - 28, 18),
            Qt.AlignmentFlag.AlignCenter,
            '蓝色尺寸线：单边顶点间长度   青色尺寸线：总宽/总高   '
            'Vn：顶点弯曲转角（括号内为几何内角）')
        painter.drawText(
            self.rect().adjusted(20, 0, -20, -12),
            Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
            '灰色：直边（0曲率）   橙色：局部弯曲窗口   '
            '最大indel/domain：%d / %d' % (
                max(self._plan['realized_max_indel_per_domain'] or [0]),
                self._plan['maximum_indel_per_domain_allowed']))


class FrameCrossSectionPreview(QWidget):
    """Physical helix cross-section and fixed inward bend direction."""

    def __init__(self, parent=None):
        super(FrameCrossSectionPreview, self).__init__(parent)
        self._plan = None
        self.setMinimumSize(360, 410)

    def setPlan(self, plan):
        self._plan = plan
        self.update()

    def paintEvent(self, unused_event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor(248, 249, 251))
        if not self._plan or not self._plan.get('rings'):
            return
        rings = self._plan['rings']
        # x is the in-plane frame normal (positive = outside); y is the
        # out-of-plane wall direction.  Each circle is a 2-nm dsDNA helix.
        points = [(float(row.get('frame_normal_offset_nm', 0.0)),
                   float(row.get('frame_binormal_offset_nm', 0.0)))
                  for row in rings]
        helix_radius_nm = 1.0
        minimum_x = min(point[0] for point in points) - helix_radius_nm
        maximum_x = max(point[0] for point in points) + helix_radius_nm
        minimum_y = min(point[1] for point in points) - helix_radius_nm
        maximum_y = max(point[1] for point in points) + helix_radius_nm
        width_nm = max(1e-6, maximum_x - minimum_x)
        height_nm = max(1e-6, maximum_y - minimum_y)
        bounds = self.rect().adjusted(58, 68, -64, -78)
        scale = min(bounds.width() / width_nm, bounds.height() / height_nm)
        centre_x, centre_y = bounds.center().x(), bounds.center().y()
        physical_centre_x = 0.5 * (minimum_x + maximum_x)
        physical_centre_y = 0.5 * (minimum_y + maximum_y)

        def mapped(point):
            return QPointF(
                centre_x + (point[0] - physical_centre_x) * scale,
                centre_y - (point[1] - physical_centre_y) * scale)

        painter.setPen(QColor(35, 53, 70))
        font = painter.font()
        font.setPointSize(9)
        painter.setFont(font)
        painter.drawText(
            QRectF(8, 7, self.width() - 16, 40),
            Qt.AlignmentFlag.AlignCenter,
            'Helix截面（与Curved Design相同弯曲面；'
            '蓝色箭头指向顶点内侧）')

        layer_colors = (QColor(210, 229, 247), QColor(177, 205, 233),
                        QColor(142, 180, 218))
        radius_px = max(5.5, helix_radius_nm * scale)
        for row, point in zip(rings, points):
            screen = mapped(point)
            layer = int(row.get('layer', 0))
            painter.setPen(QPen(QColor(63, 86, 105), 1.2))
            painter.setBrush(QBrush(layer_colors[
                min(layer, len(layer_colors) - 1)]))
            painter.drawEllipse(QRectF(
                screen.x() - radius_px, screen.y() - radius_px,
                radius_px * 2.0, radius_px * 2.0))
            if radius_px >= 9.0:
                label_font = painter.font()
                label_font.setPointSize(7)
                painter.setFont(label_font)
                painter.setPen(QColor(31, 53, 69))
                painter.drawText(
                    QRectF(screen.x() - radius_px,
                           screen.y() - radius_px,
                           radius_px * 2.0, radius_px * 2.0),
                    Qt.AlignmentFlag.AlignCenter, str(int(row['index'])))
                painter.setFont(font)

        # Frame polygons are counter-clockwise and convex.  Positive normal
        # is transported to the outside, hence the bend centre is to the left
        # (negative normal) in this cross-section view.
        arrow_start = mapped((physical_centre_x, physical_centre_y))
        arrow_tip = QPointF(max(18.0, bounds.left() - 30.0),
                            arrow_start.y())
        painter.setPen(QPen(QColor(0, 109, 204), 4.0))
        painter.drawLine(arrow_start, arrow_tip)
        painter.drawLine(arrow_tip, QPointF(arrow_tip.x() + 12,
                                            arrow_tip.y() - 7))
        painter.drawLine(arrow_tip, QPointF(arrow_tip.x() + 12,
                                            arrow_tip.y() + 7))
        painter.setPen(QColor(0, 85, 160))
        painter.drawText(QRectF(3, arrow_tip.y() - 31, 76, 22),
                         Qt.AlignmentFlag.AlignCenter, '弯曲内侧')

        # Physical dimension lines include the 1-nm radius on both sides.
        left = mapped((minimum_x, physical_centre_y)).x()
        right = mapped((maximum_x, physical_centre_y)).x()
        dimension_y = self.height() - 43
        painter.setPen(QPen(QColor(42, 66, 84), 1.0))
        painter.drawLine(QPointF(left, dimension_y),
                         QPointF(right, dimension_y))
        painter.drawLine(QPointF(left, dimension_y - 5),
                         QPointF(left, dimension_y + 5))
        painter.drawLine(QPointF(right, dimension_y - 5),
                         QPointF(right, dimension_y + 5))
        painter.drawText(QRectF(left, dimension_y + 5,
                               max(80.0, right - left), 20),
                         Qt.AlignmentFlag.AlignCenter,
                         '截面宽 %.1f nm' % width_nm)
        top = mapped((physical_centre_x, maximum_y)).y()
        bottom = mapped((physical_centre_x, minimum_y)).y()
        dimension_x = self.width() - 31
        painter.drawLine(QPointF(dimension_x, top),
                         QPointF(dimension_x, bottom))
        painter.drawLine(QPointF(dimension_x - 5, top),
                         QPointF(dimension_x + 5, top))
        painter.drawLine(QPointF(dimension_x - 5, bottom),
                         QPointF(dimension_x + 5, bottom))
        painter.save()
        painter.translate(dimension_x + 19, 0.5 * (top + bottom))
        painter.rotate(-90.0)
        painter.drawText(QRectF(-60, -10, 120, 20),
                         Qt.AlignmentFlag.AlignCenter,
                         '壁厚 %.1f nm' % height_nm)
        painter.restore()


def _convex_hull(points):
    points = sorted(set(points))
    if len(points) <= 3:
        return points

    def cross(origin, first, second):
        return ((first[0] - origin[0]) * (second[1] - origin[1]) -
                (first[1] - origin[1]) * (second[0] - origin[0]))
    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def _simplify_hull(points, maximum_vertices=12):
    points = list(points)
    if len(points) <= 3:
        return points
    span = max(max(point[0] for point in points)-min(point[0] for point in points),
               max(point[1] for point in points)-min(point[1] for point in points))
    tolerance = max(1.0, 0.008 * span)
    while len(points) > 3:
        scores = []
        for index, point in enumerate(points):
            first, second = points[index-1], points[(index+1) % len(points)]
            denominator = math.hypot(second[0]-first[0], second[1]-first[1])
            distance = (
                abs((second[1]-first[1])*point[0] -
                    (second[0]-first[0])*point[1] +
                    second[0]*first[1] - second[1]*first[0]) /
                max(denominator, 1e-9))
            scores.append((distance, index))
        distance, index = min(scores)
        if len(points) <= maximum_vertices and distance >= tolerance:
            break
        points.pop(index)
    return points


def detect_image_polygon(path):
    """Detect the dominant convex polygon without optional CV packages."""
    image = QImage(path)
    if image.isNull():
        raise ValueError('无法读取所选图片。')
    image = image.convertToFormat(QImage.Format.Format_RGB32)
    maximum = 320
    if max(image.width(), image.height()) > maximum:
        image = image.scaled(maximum, maximum,
                             Qt.AspectRatioMode.KeepAspectRatio,
                             Qt.TransformationMode.SmoothTransformation)
    values = []
    border = []
    for y in range(image.height()):
        for x in range(image.width()):
            color = QColor(image.pixel(x, y))
            value = (color.red()*299 + color.green()*587 +
                     color.blue()*114) // 1000
            values.append(value)
            if x in (0, image.width()-1) or y in (0, image.height()-1):
                border.append(value)
    # Otsu threshold.
    histogram = [0] * 256
    for value in values:
        histogram[value] += 1
    total = len(values)
    all_sum = sum(index * count for index, count in enumerate(histogram))
    background_weight = background_sum = best_variance = 0.0
    threshold = 127
    for level, count in enumerate(histogram):
        background_weight += count
        if not background_weight:
            continue
        foreground_weight = total - background_weight
        if not foreground_weight:
            break
        background_sum += level * count
        mean_background = background_sum / background_weight
        mean_foreground = (all_sum - background_sum) / foreground_weight
        variance = background_weight * foreground_weight * \
            (mean_background - mean_foreground) ** 2
        if variance > best_variance:
            best_variance, threshold = variance, level
    border_mean = sum(border) / float(len(border) or 1)
    foreground_is_dark = border_mean > threshold
    points = []
    position = 0
    for y in range(image.height()):
        for x in range(image.width()):
            value = values[position]
            position += 1
            if (value <= threshold) == foreground_is_dark:
                points.append((x, image.height()-1-y))
    if len(points) < 16:
        raise ValueError('图片中没有识别到足够清晰的封闭轮廓。')
    hull = _simplify_hull(_convex_hull(points))
    if len(hull) < 3:
        raise ValueError('无法从图片构建多边形。')
    return [(float(x), float(y)) for x, y in hull]


class FrameDesignDialog(QDialog):
    def __init__(self, parent=None):
        super(FrameDesignDialog, self).__init__(parent)
        self.setWindowTitle('Frame Design')
        self.resize(1080, 900)
        self._imageVertices = []
        self._lastPlan = None
        self._projectEdited = False

        outer = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        page = QWidget(scroll)
        layout = QVBoxLayout(page)
        scroll.setWidget(page)
        outer.addWidget(scroll, 1)

        intro = QLabel(
            '生成闭合平面多边形Frame：直边保持0曲率，只在顶点窗口内通过'
            'insertion/deletion形成目标转角。Scaffold、AutoCS和crossover'
            '完全沿用Curved Design稳定规则，Frame不会新增、删除或移动crossover。',
            self)
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.shapeBox = QComboBox(self)
        for key, label in FRAME_PRESETS:
            self.shapeBox.addItem(label, key)
        form.addRow('Frame形状：', self.shapeBox)

        self.sideSpin = QDoubleSpinBox(self)
        self.sideSpin.setRange(15.0, 1000.0)
        self.sideSpin.setValue(200.0)
        self.sideSpin.setSuffix(' nm')
        form.addRow('正多边形边长：', self.sideSpin)
        self.firstSideSpin = QDoubleSpinBox(self)
        self.firstSideSpin.setRange(15.0, 1000.0)
        self.firstSideSpin.setValue(200.0)
        self.firstSideSpin.setSuffix(' nm')
        form.addRow('第一边长：', self.firstSideSpin)
        self.secondSideSpin = QDoubleSpinBox(self)
        self.secondSideSpin.setRange(15.0, 1000.0)
        self.secondSideSpin.setValue(160.0)
        self.secondSideSpin.setSuffix(' nm')
        form.addRow('第二边长：', self.secondSideSpin)
        self.angleSpin = QDoubleSpinBox(self)
        self.angleSpin.setRange(15.0, 165.0)
        self.angleSpin.setValue(70.0)
        self.angleSpin.setSuffix('°')
        form.addRow('平行四边形内角：', self.angleSpin)

        imageRow = QHBoxLayout()
        self.imageButton = QPushButton('上传轮廓图片…', self)
        self.imageLabel = QLabel('尚未上传', self)
        self.imageReferenceSpin = QDoubleSpinBox(self)
        self.imageReferenceSpin.setRange(15.0, 1000.0)
        self.imageReferenceSpin.setValue(200.0)
        self.imageReferenceSpin.setSuffix(' nm')
        imageRow.addWidget(self.imageButton)
        imageRow.addWidget(self.imageLabel, 1)
        imageRow.addWidget(QLabel('检测后第一条边：', self))
        imageRow.addWidget(self.imageReferenceSpin)
        form.addRow('任意多边形：', imageRow)
        self.imageLengthsEdit = QLineEdit(self)
        self.imageLengthsEdit.setPlaceholderText(
            '识别图片后输入每条边长，例如 80, 95, 70, 110')
        form.addRow('图片多边形各边长（nm）：', self.imageLengthsEdit)

        self.latticeBox = QComboBox(self)
        self.latticeBox.addItem('Square（8-bp domain）', 'square')
        self.latticeBox.addItem('Honeycomb（7-bp domain）', 'honeycomb')
        form.addRow('点阵通道：', self.latticeBox)
        self.densityBox = QComboBox(self)
        self.densityBox.addItem('最大合法密度（原生周期）', ('periodic', 1))
        self.densityBox.addItem('2×原生周期', ('periodic', 2))
        self.densityBox.addItem('最低合法密度', ('minimum', 0))
        form.addRow('Scaffold crossover密度：', self.densityBox)
        self.heightSpin = QDoubleSpinBox(self)
        self.heightSpin.setRange(5.2, 200.0)
        self.heightSpin.setDecimals(1)
        self.heightSpin.setValue(31.2)
        self.heightSpin.setSuffix(' nm')
        form.addRow('边束截面总高度：', self.heightSpin)
        self.layersSpin = QSpinBox(self)
        self.layersSpin.setRange(1, 3)
        self.layersSpin.setValue(1)
        form.addRow('边束壁层数：', self.layersSpin)

        self.maximumIndelBox = QComboBox(self)
        self.maximumIndelBox.addItem('最大 ±1 / domain', 1)
        self.maximumIndelBox.addItem('最大 ±2 / domain', 2)
        self.maximumIndelBox.addItem('最大 ±3 / domain', 3)
        self.maximumIndelBox.setCurrentIndex(2)
        form.addRow('最大增删密度：', self.maximumIndelBox)
        self.bendModeBox = QComboBox(self)
        self.bendModeBox.addItem('Auto：按最大indel/domain求最短窗口', 'auto')
        self.bendModeBox.addItem('自定义弯曲基础长度', 'custom')
        form.addRow('顶点弯曲长度：', self.bendModeBox)
        self.bendLengthSpin = QSpinBox(self)
        self.bendLengthSpin.setRange(7, 2000)
        self.bendLengthSpin.setValue(32)
        self.bendLengthSpin.setSuffix(' bp')
        form.addRow('自定义长度（自动吸附完整domain）：',
                    self.bendLengthSpin)
        self.nameEdit = QLineEdit('frame', self)
        form.addRow('设计名称：', self.nameEdit)
        projectRow = QHBoxLayout()
        self.projectEdit = QLineEdit(self)
        self.projectButton = QPushButton('选择…', self)
        projectRow.addWidget(self.projectEdit, 1)
        projectRow.addWidget(self.projectButton)
        form.addRow('项目文件夹：', projectRow)
        layout.addLayout(form)

        self.reportLabel = QLabel(self)
        self.reportLabel.setWordWrap(True)
        self.reportLabel.setMargin(9)
        layout.addWidget(self.reportLabel)

        previewRow = QHBoxLayout()
        frameColumn = QVBoxLayout()
        frameTitle = QLabel('Frame中性轴与顶点弯曲窗口：', self)
        frameColumn.addWidget(frameTitle)
        self.preview = FramePreview(self)
        frameColumn.addWidget(self.preview)
        previewRow.addLayout(frameColumn, 3)
        crossColumn = QVBoxLayout()
        crossTitle = QLabel('边束helix截面与弯曲方向：', self)
        crossColumn.addWidget(crossTitle)
        self.crossSectionPreview = FrameCrossSectionPreview(self)
        crossColumn.addWidget(self.crossSectionPreview)
        previewRow.addLayout(crossColumn, 2)
        layout.addLayout(previewRow)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel, self)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText(
            '运行Frame Design')
        outer.addWidget(self.buttons)

        for widget in (self.shapeBox, self.latticeBox, self.densityBox,
                       self.maximumIndelBox, self.bendModeBox):
            widget.currentIndexChanged.connect(self._refresh)
        for widget in (self.sideSpin, self.firstSideSpin,
                       self.secondSideSpin, self.angleSpin,
                       self.imageReferenceSpin, self.heightSpin,
                       self.layersSpin, self.bendLengthSpin):
            widget.valueChanged.connect(self._refresh)
        self.imageButton.clicked.connect(self._chooseImage)
        self.imageLengthsEdit.textChanged.connect(self._refresh)
        self.nameEdit.textChanged.connect(self._nameChanged)
        self.projectEdit.textEdited.connect(self._projectChanged)
        self.projectButton.clicked.connect(self._chooseProject)
        self.buttons.accepted.connect(self._acceptIfValid)
        self.buttons.rejected.connect(self.reject)
        self._refresh()

    def _defaultProject(self):
        return os.path.join(QDir.homePath(), 'Desktop',
                            safe_name(self.nameEdit.text() or 'frame'))

    def _nameChanged(self, unused):
        if not self._projectEdited:
            self.projectEdit.setText(self._defaultProject())

    def _projectChanged(self, unused):
        self._projectEdited = True

    def _chooseProject(self):
        parent = QFileDialog.getExistingDirectory(
            self, '选择Frame Design项目的父文件夹',
            os.path.dirname(self.projectEdit.text()) or QDir.homePath())
        if parent:
            self._projectEdited = True
            self.projectEdit.setText(os.path.join(
                str(parent), safe_name(self.nameEdit.text() or 'frame')))

    def _chooseImage(self):
        path, unused_filter = QFileDialog.getOpenFileName(
            self, '上传多边形轮廓图片', QDir.homePath(),
            'Images (*.png *.jpg *.jpeg *.bmp *.tif *.tiff)')
        if not path:
            return
        try:
            self._imageVertices = detect_image_polygon(path)
        except Exception as error:
            QMessageBox.critical(self, 'Frame Design', str(error))
            return
        self.imageLabel.setText('%s；识别%d个顶点' %
                                (os.path.basename(path),
                                 len(self._imageVertices)))
        source_lengths = [math.hypot(
            self._imageVertices[(index+1) % len(self._imageVertices)][0] -
            point[0],
            self._imageVertices[(index+1) % len(self._imageVertices)][1] -
            point[1]) for index, point in enumerate(self._imageVertices)]
        reference = source_lengths[0] or 1.0
        self.imageLengthsEdit.setText(', '.join(
            '%.1f' % (200.0 * value/reference) for value in source_lengths))
        self.shapeBox.setCurrentIndex(self.shapeBox.findData('image'))
        self._refresh()

    def _rawSpec(self):
        density_mode, density_multiple = self.densityBox.currentData()
        image_lengths = []
        for value in self.imageLengthsEdit.text().replace('，', ',').split(','):
            if value.strip():
                image_lengths.append(float(value.strip()))
        return {
            'frame_shape': str(self.shapeBox.currentData()),
            'side_nm': self.sideSpin.value(),
            'first_side_nm': self.firstSideSpin.value(),
            'second_side_nm': self.secondSideSpin.value(),
            'corner_angle_degrees': self.angleSpin.value(),
            'image_vertices': list(self._imageVertices),
            'image_reference_side_nm': self.imageReferenceSpin.value(),
            'image_side_lengths_nm': image_lengths,
            'lattice': str(self.latticeBox.currentData()),
            'layers': self.layersSpin.value(),
            'cross_section_height_nm': self.heightSpin.value(),
            'bend_length_mode': str(self.bendModeBox.currentData()),
            'bend_length_bp': self.bendLengthSpin.value(),
            'maximum_indel_per_domain': int(
                self.maximumIndelBox.currentData()),
            'scaffold_crossover_density_mode': density_mode,
            'scaffold_crossover_density_multiple': density_multiple,
            'name': self.nameEdit.text().strip() or 'frame',
            'project_root': self.projectEdit.text().strip() or
                            self._defaultProject()}

    def _refresh(self, unused=None):
        shape = str(self.shapeBox.currentData())
        regular = shape in ('triangle', 'square', 'pentagon', 'hexagon')
        two_side = shape in ('rectangle', 'parallelogram')
        self.sideSpin.setEnabled(regular)
        self.firstSideSpin.setEnabled(two_side)
        self.secondSideSpin.setEnabled(two_side)
        self.angleSpin.setEnabled(shape == 'parallelogram')
        self.imageButton.setEnabled(shape == 'image')
        self.imageReferenceSpin.setEnabled(shape == 'image')
        self.imageLengthsEdit.setEnabled(shape == 'image')
        self.bendLengthSpin.setEnabled(
            self.bendModeBox.currentData() == 'custom')
        if not self._projectEdited:
            self.projectEdit.setText(self._defaultProject())
        try:
            plan = plan_frame(self._rawSpec())
            self._lastPlan = plan
            self.preview.setPlan(plan)
            self.crossSectionPreview.setPlan(plan)
            targets = [int(row['frame_target_bp'])
                       for row in plan['rings']]
            required = sum(targets)
            capacity_ok = required <= CURVED_SCAFFOLD_MAX_BASES
            fully_feasible = bool(plan['feasible'] and capacity_ok)
            if fully_feasible:
                status = '合法，可生成'
            elif not plan['feasible']:
                status = '不可生成：' + '；'.join(plan['failure_reasons'])
            else:
                status = ('不可生成：scaffold需要%d nt，超过%d nt上限' %
                          (required, CURVED_SCAFFOLD_MAX_BASES))
            vertices = plan['vertices_nm']
            span_x = max(point[0] for point in vertices) - \
                min(point[0] for point in vertices)
            span_y = max(point[1] for point in vertices) - \
                min(point[1] for point in vertices)
            normal_offsets = [float(row['frame_normal_offset_nm'])
                              for row in plan['rings']]
            binormal_offsets = [float(row['frame_binormal_offset_nm'])
                                for row in plan['rings']]
            cross_width = max(normal_offsets) - min(normal_offsets) + 2.0
            cross_thickness = (max(binormal_offsets) -
                               min(binormal_offsets) + 2.0)
            density_mode, density_multiple = self.densityBox.currentData()
            native = int(plan['native_period_bp'])
            parent_period = int(plan.get('parent_period_bp', native))
            density_text = ('最低合法密度' if density_mode == 'minimum'
                            else '目标1/%d bp' %
                            (native * int(density_multiple)))
            indel_totals = [int(row.get('frame_total_indel', 0))
                            for row in plan['rings']]
            rows = []
            for index in range(len(plan['vertices_nm'])):
                rows.append(
                    'V%d：转角%.1f°，%d domain / %d bp，R=%.2f nm，'
                    '最大indel=%d，平均密度=%.2f/domain，实际峰值=%d/domain' %
                    (index+1, plan['turn_angles_degrees'][index],
                     plan['bend_domain_count'][index],
                     plan['bend_length_bp'][index],
                     plan['bend_radius_nm'][index],
                     plan['maximum_abs_indel_by_vertex'][index],
                     plan['average_max_indel_per_domain'][index],
                     plan['realized_max_indel_per_domain'][index]))
            self.reportLabel.setText(
                '<b>%s</b><br>'
                '实际Frame中性轴外包络：%.2f × %.2f nm；边数：%d；'
                '中性轴原生闭环：%d bp（%.2f nm，吸附%d-bp周期）。<br>'
                '边束截面设置高度：%.2f nm；实际截面：宽%.2f nm × '
                '壁厚%.2f nm；helix：%d；层数：%d。蓝色箭头指向所有'
                '凸顶点的弯曲内侧。<br>'
                '各helix目标闭环：%d–%d bp；相对中性轴indel：%+d–%+d；'
                '最大实际indel/domain：%d（硬上限±%d）；%s。<br>'
                'scaffold容量：%d/%d nt。所有helix先使用同一中性轴原生'
                '闭环，顶点窗口只重分配indel，不新增、删除或移动AutoCS '
                'crossover。<br>'
                '<b>为什么截面高度会影响原生闭环：</b>截面越高，内外侧'
                'helix在同一转角下的长度差越大；为保持所选最大indel/domain，自动'
                '弯曲窗口必须加长，圆角替代的路径长度随之变化，最终中性轴'
                '闭环再吸附到完整%d-bp周期，因此长度会离散改变。<br>%s' %
                (status, span_x, span_y, len(vertices),
                 plan['nominal_perimeter_bp'], plan['nominal_perimeter_nm'],
                 parent_period, float(self.heightSpin.value()), cross_width,
                 cross_thickness, len(plan['rings']),
                 int(self.layersSpin.value()), min(targets), max(targets),
                 min(indel_totals), max(indel_totals),
                 max(plan['realized_max_indel_per_domain'] or [0]),
                 plan['maximum_indel_per_domain_allowed'], density_text,
                 required, CURVED_SCAFFOLD_MAX_BASES, parent_period,
                 '<br>'.join(rows)))
            if fully_feasible:
                style = ('QLabel { color:#1f6130; background:#eff8f1; '
                         'border:1px solid #9ac8a5; border-radius:5px; }')
            else:
                style = ('QLabel { color:#8f1d1d; background:#fff1f0; '
                         'border:1px solid #e5a3a0; border-radius:5px; }')
            self.reportLabel.setStyleSheet(style)
            self.reportLabel.setProperty('frameRequiredBases', required)
            self.buttons.button(
                QDialogButtonBox.StandardButton.Ok).setEnabled(
                    fully_feasible)
        except Exception as error:
            self._lastPlan = None
            self.preview.setPlan(None)
            self.crossSectionPreview.setPlan(None)
            self.reportLabel.setText('<b>尚不可生成：</b>%s' % error)
            self.reportLabel.setProperty('frameRequiredBases', 0)
            self.reportLabel.setStyleSheet(
                'QLabel {background:#fff0f0; border:1px solid #c87979;}')
            self.buttons.button(
                QDialogButtonBox.StandardButton.Ok).setEnabled(False)

    def _acceptIfValid(self):
        self._refresh()
        required = int(self.reportLabel.property(
            'frameRequiredBases') or 0)
        if (self._lastPlan and self._lastPlan['feasible'] and
                required <= CURVED_SCAFFOLD_MAX_BASES):
            self.accept()

    def spec(self):
        return self._rawSpec()

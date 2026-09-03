"""Lightweight interactive previews with no OpenGL dependency."""

from __future__ import annotations

import json
import math
from pathlib import Path

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import (
    QColor, QFont, QPainter, QPainterPath, QPen, QPixmap, QPolygonF)
from PyQt6.QtWidgets import QSizePolicy, QWidget

from .lattice_preview_geometry import (
    lattice_graph, rotated_graph_outside_square)


PREVIEW_BACKGROUND = QColor("#070b10")
PREVIEW_CONTENT_TOP = 12.0
PREVIEW_CONTENT_BOTTOM = 42.0
# Use one apparent-height normalization for the default 2D and 3D models.
# The extra room prevents the rotated 3D structure from clipping.
PREVIEW_MODEL_HEIGHT_FACTOR = 1.48


def _right_handed_preview_rotation(angle_deg):
    """Map the design sign convention to the preview coordinate system."""
    # The on-screen/projection axes reverse the apparent handedness.  Keep
    # all design and calibration values unchanged and invert only drawing.
    return -float(angle_deg)


def _default_seed_cells():
    return {(row, col) for row in range(8) for col in range(8)
            if not (2 <= row <= 5 and 2 <= col <= 5)}


def _rotate(point, yaw, pitch):
    x, y, z = point
    cy, sy = math.cos(yaw), math.sin(yaw)
    x, y = cy*x-sy*y, sy*x+cy*y
    cp, sp = math.cos(pitch), math.sin(pitch)
    y, z = cp*y-sp*z, sp*y+cp*z
    return x, y, z


def _draw_scale_bar(painter, width, height, pixels_per_nm):
    """Draw a legible 10-nm overlay at the current physical screen scale."""
    length_px = 10.0*float(pixels_per_nm)
    card_width = max(88.0, length_px+20.0)
    card_right = float(width)-12.0
    card_left = card_right-card_width
    center = (card_left+card_right)/2.0
    left = center-length_px/2.0
    right = center+length_px/2.0
    baseline = float(height)-32.0
    background = QRectF(card_left, baseline-24.0, card_width, 32.0)
    painter.save()
    painter.setPen(QPen(QColor(40, 49, 58, 70), .8))
    painter.setBrush(QColor(255, 255, 255, 220))
    painter.drawRoundedRect(background, 4.0, 4.0)
    pen = QPen(QColor("#111820"), 2.4)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.drawLine(QPointF(left, baseline), QPointF(right, baseline))
    painter.drawLine(QPointF(left, baseline-4.0),
                     QPointF(left, baseline+4.0))
    painter.drawLine(QPointF(right, baseline-4.0),
                     QPointF(right, baseline+4.0))
    scale_font = QFont("Arial")
    scale_font.setPixelSize(11)
    scale_font.setWeight(QFont.Weight.DemiBold)
    painter.setFont(scale_font)
    painter.drawText(QRectF(card_left+4.0, baseline-22.0,
                            card_width-8.0, 16.0),
                     Qt.AlignmentFlag.AlignCenter, "10.0 nm")
    painter.restore()
    return length_px


def _draw_model_coverage_frame(painter, width, height):
    """Outline the complete drawable extent of an individual preview."""
    # This is the actual widget boundary: panning or zooming may place model
    # geometry anywhere inside it.  Do not reuse the default-fit safe margins
    # (which reserve room for annotations), because those are not the maximum
    # drawable/model extent requested by the UI.
    rect = QRectF(
        .5, .5,
        max(1.0, float(width)-1.0),
        max(1.0, float(height)-1.0))
    painter.save()
    pen = QPen(QColor("#52616f"), 1.0, Qt.PenStyle.DashLine)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRect(rect)
    painter.restore()
    return rect


def _preview_text_font(weight=QFont.Weight.Normal):
    """Return the shared 11-pixel body font used by preview annotations."""
    font = QFont("Arial")
    font.setPixelSize(11)
    font.setWeight(weight)
    return font


def _format_preview_bp(value):
    """Format an axial length compactly without hiding fractional spacing."""
    numeric = float(value)
    if math.isclose(numeric, round(numeric), abs_tol=1e-6):
        return "%d bp" % int(round(numeric))
    return "%.1f bp" % numeric


def _draw_double_arrow(painter, first, second, color):
    """Draw a cosmetic screen-space dimension line with two arrowheads."""
    dx = second.x()-first.x()
    dy = second.y()-first.y()
    length = math.hypot(dx, dy)
    if length < 2.0:
        return
    ux, uy = dx/length, dy/length
    px, py = -uy, ux
    size = max(3.5, min(6.0, length*.16))
    pen = QPen(QColor(color), 1.15)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(QColor(color))
    painter.drawLine(first, second)
    for tip, sign in ((first, 1.0), (second, -1.0)):
        base = QPointF(
            tip.x()+sign*ux*size,
            tip.y()+sign*uy*size)
        painter.drawPolygon(QPolygonF([
            tip,
            QPointF(base.x()+px*size*.55, base.y()+py*size*.55),
            QPointF(base.x()-px*size*.55, base.y()-py*size*.55),
        ]))


class SeedCrossSectionPicker(QWidget):
    """cadnano-like, number-free Square-grid Seed cross-section picker."""

    selectionChanged = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        # A smaller vertical minimum lets the surrounding selection page
        # reserve independent rows for status text and buttons.  The grid
        # remains square because _grid_geometry uses the available short side.
        self.setMinimumSize(180, 145)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self._size = 8
        self._selected = _default_seed_cells()
        self._interactive = True
        self.setToolTip("点击圆圈添加或移除Seed helix；仅使用Square网格。")

    def set_interactive(self, enabled):
        """Enable expert point editing without dimming a preset preview."""
        self._interactive = bool(enabled)
        self.setToolTip(
            "点击圆圈添加或移除Seed helix；仅使用Square网格。"
            if self._interactive else
            "当前截面由预设选择；请使用上方截面选项切换。")

    def cells(self):
        return [list(cell) for cell in sorted(self._selected)]

    def set_cells(self, cells):
        selected = {(int(row), int(col)) for row, col in cells
                    if 0 <= int(row) < self._size and
                    0 <= int(col) < self._size}
        if selected != self._selected:
            self._selected = selected
            self.selectionChanged.emit()
            self.update()

    def reset_default(self):
        self.set_cells(_default_seed_cells())

    def _grid_geometry(self):
        side = max(1.0, min(self.width(), self.height()) - 34.0)
        step = side / self._size
        radius = min(16.0, step * .39)
        left = (self.width() - side) / 2.0 + step / 2.0
        top = (self.height() - side) / 2.0 + step / 2.0
        return left, top, step, radius

    def mousePressEvent(self, event):
        if not self._interactive:
            return
        if event.button() != Qt.MouseButton.LeftButton:
            return
        left, top, step, radius = self._grid_geometry()
        col = int(round((event.position().x() - left) / step))
        row = int(round((event.position().y() - top) / step))
        if not (0 <= row < self._size and 0 <= col < self._size):
            return
        center = QPointF(left + col * step, top + row * step)
        delta = event.position() - center
        if delta.x() * delta.x() + delta.y() * delta.y() > \
                (radius * 1.25) ** 2:
            return
        cell = (row, col)
        if cell in self._selected:
            self._selected.remove(cell)
        else:
            self._selected.add(cell)
        self.selectionChanged.emit()
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#f5f7fa"))
        left, top, step, radius = self._grid_geometry()
        for row in range(self._size):
            for col in range(self._size):
                center = QPointF(left + col * step, top + row * step)
                selected = (row, col) in self._selected
                painter.setPen(QPen(QColor(
                    "#155da8" if selected else "#8794a3"),
                    1.8 if selected else 1.0))
                painter.setBrush(QColor(
                    "#4c9be8" if selected else "#e2e7ec"))
                painter.drawEllipse(center, radius, radius)


class BilayerPreview(QWidget):
    """Drag-to-rotate painter preview of a Seed-S templated bilayer."""

    # Square and Kagome share one physical square preview field.  Kagome's
    # larger lattice constant changes helix count, not the displayed extent.
    # Keep the 3-D SST field physically larger than the Seed: doubling the
    # half extent doubles the complete square side length without scaling the
    # Seed rods themselves.
    LATTICE_HALF_EXTENT_NM = 11.0*2.8
    HELIX_DIAMETER_NM = 2.0

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(360, 260)
        self.setMouseTracking(True)
        self._preview_background = QColor(PREVIEW_BACKGROUND)
        # The right preview opens in a true orthographic side view.  Users may
        # still rotate freely; double-click restores this same side view.
        self._yaw = 0.0
        self._pitch = 0.0
        self._zoom = 1.0
        self._last = None
        self._drag_button = None
        self._pan = QPointF(0.0, 0.0)
        self._angle = 3.2967555036483183
        self._period_nm = None
        self._sst_z1 = 128
        self._nominal_spacing = 32
        self._spacing = 32.0
        self._sst_z3 = 128
        self._seed_z1 = 128
        self._seed_z3 = 128
        self._seed_overlap_z1 = 128
        self._seed_overlap_z3 = 128
        self._seed_total = 288
        self._sst_layer_ranges = ((48, 175), (208, 335))
        self._seed_partition_ranges = (
            (48, 175), (176, 207), (208, 335))
        self._partition_mode = "linked"
        self._lattice_constant = 2.8
        self._layer_lattices = ("square", "square")
        self._layer_constants = (2.8, 2.8)
        self._symmetry_label = "Square–Square"
        self._seed_cells = _default_seed_cells()
        self._seed_helix_spacing = 2.8
        self._world_scale = 5.0
        self._last_scale_bar_px = 0.0
        self.setToolTip("左键拖动旋转；右键拖动平移；滚轮缩放；双击复位")

    def set_preview_background(self, color):
        self._preview_background = QColor(color)
        self.update()

    def set_design(self, project):
        self._angle = float(project.prediction["reported_angle_deg"])
        self._period_nm = project.prediction.get(
            "predicted_moire_period_nm")
        self._sst_z1 = int(getattr(
            project.settings, "sst_growth_bp_z1",
            project.settings.growth_bp_z1))
        self._nominal_spacing = int(project.settings.spacer_bp_z2)
        self._spacing = float(project.prediction.get(
            "actual_z2_spacing_bp", project.settings.spacer_bp_z2))
        self._sst_z3 = int(getattr(
            project.settings, "sst_growth_bp_z3",
            project.settings.growth_bp_z3))
        partition = project.prediction.get("preview_seed_partition", {})
        self._seed_z1 = int(partition.get("z1_bp", 128))
        self._seed_z3 = int(partition.get("z3_bp", 128))
        self._seed_overlap_z1 = int(partition.get(
            "sst_overlap_z1_bp", min(self._seed_z1, self._sst_z1)))
        self._seed_overlap_z3 = int(partition.get(
            "sst_overlap_z3_bp", min(self._seed_z3, self._sst_z3)))
        self._seed_total = int(partition.get("total_bp", 288))
        self._sst_layer_ranges = tuple(
            tuple(map(int, item)) for item in partition.get(
                "sst_layer_ranges", ((48, 175), (208, 335))))
        self._seed_partition_ranges = tuple(
            tuple(map(int, item)) for item in partition.get(
                "seed_partition_ranges",
                ((48, 175), (176, 207), (208, 335))))
        self._partition_mode = (
            "linked" if partition.get("linked", True) else "independent")
        self._lattice_constant = float(project.settings.lattice_constant_nm)
        self._layer_lattices = tuple(project.prediction.get(
            "layer_lattice_types", ("square", "square")))
        self._layer_constants = tuple(project.prediction.get(
            "layer_lattice_constants_nm", (
                self._lattice_constant, self._lattice_constant)))
        self._symmetry_label = {
            "square_square_c4": "Square–Square",
            "kagome_kagome": "Kagome–Kagome",
            "square_kagome": "Square–Kagome",
        }.get(project.settings.lattice_symmetry,
              project.settings.lattice_symmetry)
        self._seed_cells = {tuple(map(int, cell)) for cell in
                            project.seed_plan.get(
                                "cross_section_cells",
                                _default_seed_cells())}
        self.update()

    def mousePressEvent(self, event):
        if event.button() in (Qt.MouseButton.LeftButton,
                              Qt.MouseButton.RightButton):
            self._last = event.position()
            self._drag_button = event.button()

    def mouseMoveEvent(self, event):
        if self._last is None:
            return
        delta = event.position()-self._last
        self._last = event.position()
        if self._drag_button == Qt.MouseButton.RightButton:
            self._pan += delta
        else:
            # Grab-and-drag interaction: the model follows the cursor rather
            # than orbiting in the opposite camera direction.
            self._yaw -= delta.x()*0.008
            # ±89° permits a true cross-section/top-down view.
            self._pitch = max(-1.553, min(1.553,
                self._pitch-delta.y()*0.008))
        self.update()

    def mouseReleaseEvent(self, event):
        self._last = None
        self._drag_button = None

    def mouseDoubleClickEvent(self, event):
        self._yaw = 0.0
        self._pitch = 0.0
        self._zoom = 1.0
        self._pan = QPointF(0.0, 0.0)
        self.update()

    def wheelEvent(self, event):
        self._zoom *= 1.12 if event.angleDelta().y() > 0 else 1/1.12
        self._zoom = max(.45, min(2.8, self._zoom))
        self.update()

    def _project(self, point):
        x, y, z = _rotate(point, self._yaw, self._pitch)
        scale = self._world_scale*self._zoom
        content_top = PREVIEW_CONTENT_TOP
        content_bottom = PREVIEW_CONTENT_BOTTOM
        content_height = max(1.0, self.height()-content_top-content_bottom)
        # Orthographic projection: depth controls draw order/fading only and
        # never changes the apparent size of an object.
        # Legends now live below the canvas, so the side-view model remains
        # centred in its full black preview field.
        return QPointF(self.width()*.50+self._pan.x()+x*scale,
                       content_top+content_height*.50+
                       self._pan.y()-z*scale), y

    def _line3(self, painter, first, second, pen):
        a, da = self._project(first)
        b, db = self._project(second)
        faded = QPen(pen)
        depth = (da+db)/2.0
        color = faded.color()
        color.setAlpha(max(70, min(255, int(215-depth*2.0))))
        faded.setColor(color)
        painter.setPen(faded)
        painter.drawLine(a, b)

    def _point3(self, painter, point, color, radius=2.6):
        mapped, depth = self._project(point)
        faded = QColor(color)
        faded.setAlpha(max(85, min(255, int(225-depth*2.0))))
        painter.setPen(QPen(QColor(7, 11, 16, faded.alpha()), 1.0))
        painter.setBrush(faded)
        painter.drawEllipse(mapped, radius, radius)

    def _draw_grid(self, painter, z, rotation_deg, color,
                   lattice_type="square", lattice_constant=None):
        step = float(lattice_constant or self._lattice_constant)
        points, edges = lattice_graph(
            lattice_type, step, self.LATTICE_HALF_EXTENT_NM)
        points, edges = rotated_graph_outside_square(
            points, edges, rotation_deg,
            4.0*self._seed_helix_spacing+self.HELIX_DIAMETER_NM/2.0)
        pen_color = QColor(color)
        pen_color.setAlpha(190)
        pen = QPen(pen_color, 1.35)
        pen.setCosmetic(True)
        rotated = [(x, y, z) for x, y in points]
        for first, second in edges:
            self._line3(painter, rotated[first], rotated[second], pen)
        for point in rotated:
            self._point3(painter, point, color)

    def _draw_lattice_volume(self, painter, low, high, rotation_deg, color,
                             lattice_type="square", lattice_constant=None):
        """Draw cross-section networks enclosed by four translucent faces."""
        lattice_constant = float(lattice_constant or self._lattice_constant)
        self._draw_lattice_side_faces(
            painter, low, high, rotation_deg, color)
        self._draw_grid(
            painter, low, rotation_deg, color,
            lattice_type, lattice_constant)
        self._draw_grid(
            painter, high, rotation_deg, color,
            lattice_type, lattice_constant)

    def _draw_lattice_side_faces(self, painter, low, high, rotation_deg,
                                 color):
        """Close the outer SST volume with faces instead of many Z lines."""
        half = self.LATTICE_HALF_EXTENT_NM
        angle = math.radians(rotation_deg)
        cosine, sine = math.cos(angle), math.sin(angle)
        corners = []
        for x, y in ((-half, -half), (half, -half),
                     (half, half), (-half, half)):
            corners.append((cosine*x-sine*y, sine*x+cosine*y))
        faces = []
        for index in range(4):
            first = corners[index]
            second = corners[(index+1) % 4]
            world = [
                (first[0], first[1], low),
                (second[0], second[1], low),
                (second[0], second[1], high),
                (first[0], first[1], high),
            ]
            projected = [self._project(point) for point in world]
            faces.append((
                sum(depth for unused_point, depth in projected)/4.0,
                QPolygonF([point for point, unused_depth in projected])))
        painter.save()
        outline = QColor(color)
        outline.setAlpha(135)
        fill = QColor(color)
        fill.setAlpha(34)
        pen = QPen(outline, 1.0)
        pen.setCosmetic(True)
        painter.setPen(pen)
        painter.setBrush(fill)
        for unused_depth, polygon in sorted(faces, key=lambda item: item[0]):
            painter.drawPolygon(polygon)
        painter.restore()

    def _draw_rod_path(self, painter, world_points, scale,
                       body_color, highlight_color):
        """Match the cadnano Twist and Bend cylinder rendering exactly."""
        projected = [self._project(point)[0] for point in world_points]
        if len(projected) < 2:
            return
        path = QPainterPath()
        path.moveTo(projected[0])
        for point in projected[1:]:
            path.lineTo(point)
        rod_width = max(2.0, 2.0*scale)
        outline_width = rod_width+max(1.2, .18*scale)

        def draw(color, width):
            pen = QPen(color)
            pen.setWidthF(width)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(path)

        draw(QColor(8, 12, 16, 245), outline_width)
        draw(QColor(body_color), rod_width)
        draw(QColor(highlight_color), max(1.0, rod_width*.52))

    def _extreme_anchor(self, z, half_extent, rotation_deg, side):
        """Return the left/right-most projected corner at one axial boundary."""
        angle = math.radians(rotation_deg)
        cosine, sine = math.cos(angle), math.sin(angle)
        projected = []
        for x, y in ((-half_extent, -half_extent),
                     (half_extent, -half_extent),
                     (half_extent, half_extent),
                     (-half_extent, half_extent)):
            world = (cosine*x-sine*y, sine*x+cosine*y, z)
            projected.append(self._project(world)[0])
        selector = min if side == "left" else max
        return selector(projected, key=lambda point: point.x())

    def _draw_dimension(self, painter, first_anchor, second_anchor,
                        dimension_x, text, color, side):
        """Draw projected extension lines and a horizontal dimension label."""
        first_y = first_anchor.y()
        second_y = second_anchor.y()
        extension_pen = QPen(QColor(196, 211, 223, 205), 1.0,
                             Qt.PenStyle.DashLine)
        extension_pen.setCosmetic(True)
        painter.save()
        painter.setPen(extension_pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawLine(first_anchor, QPointF(dimension_x, first_y))
        painter.drawLine(second_anchor, QPointF(dimension_x, second_y))
        _draw_double_arrow(
            painter, QPointF(dimension_x, first_y),
            QPointF(dimension_x, second_y), "#f3f6f9")

        painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
        painter.setPen(QColor(color))
        middle_y = (first_y+second_y)/2.0
        if side == "left":
            text_rect = QRectF(
                7.0, middle_y-18.0,
                max(24.0, dimension_x-16.0), 36.0)
            alignment = (Qt.AlignmentFlag.AlignRight |
                         Qt.AlignmentFlag.AlignVCenter)
        else:
            text_rect = QRectF(
                dimension_x+9.0, middle_y-18.0,
                max(24.0, self.width()-dimension_x-16.0), 36.0)
            alignment = (Qt.AlignmentFlag.AlignLeft |
                         Qt.AlignmentFlag.AlignVCenter)
        painter.drawText(text_rect, alignment, text)
        painter.restore()

    def _draw_axial_dimensions(
            self, painter, sst_first, sst_second,
            seed_first, seed_spacing, seed_second,
            to_world, visual_twist):
        """Project model-anchored SST/Seed dimensions into screen space."""
        sst_segments = (
            (sst_first[0], sst_first[1], 0.0,
             "1st layer\n%s" % _format_preview_bp(self._sst_z1),
             "#2a78d1"),
            (sst_first[1], sst_second[0], visual_twist*.5,
             "Spacing\n%s" % _format_preview_bp(self._spacing),
             "#ffffff"),
            (sst_second[0], sst_second[1], visual_twist,
             "2nd layer\n%s" % _format_preview_bp(self._sst_z3),
             "#d65b74"),
        )
        seed_first_overlap = (
            max(sst_first[0], seed_first[0]),
            min(sst_first[1], seed_first[1]))
        seed_second_overlap = (
            max(sst_second[0], seed_second[0]),
            min(sst_second[1], seed_second[1]))
        seed_segments = (
            (seed_first_overlap[0], seed_first_overlap[1], 0.0,
             "Z1\n%s" % _format_preview_bp(self._seed_overlap_z1),
             "#2a78d1"),
            (seed_spacing[0], seed_spacing[1], visual_twist*.5,
             "Z2\n%s" % _format_preview_bp(self._spacing),
             "#ffffff"),
            (seed_second_overlap[0], seed_second_overlap[1], visual_twist,
             "Z3\n%s" % _format_preview_bp(self._seed_overlap_z3),
             "#d65b74"),
        )

        def anchors(segments, half_extent, side):
            rows = []
            for low, high, rotation, text, color in segments:
                low_anchor = self._extreme_anchor(
                    to_world(low), half_extent, rotation, side)
                high_anchor = self._extreme_anchor(
                    to_world(high), half_extent, rotation, side)
                rows.append((low_anchor, high_anchor, text, color))
            return rows

        sst_rows = anchors(
            sst_segments, self.LATTICE_HALF_EXTENT_NM, "left")
        # Seed values are referenced on the right, but the guide lines begin
        # at the complete structure silhouette.  They therefore never cross
        # the SST volume or the central Seed bundle.
        seed_rows = anchors(
            seed_segments, self.LATTICE_HALF_EXTENT_NM, "right")
        left_edge = min(
            min(first.x(), second.x())
            for first, second, unused_text, unused_color in sst_rows)
        right_edge = max(
            max(first.x(), second.x())
            for first, second, unused_text, unused_color in seed_rows)
        # Keep the measurement rails decisively outside the rendered
        # silhouette.  A compact 72-pixel edge reserve still fits the two-line
        # labels at the default pane width, while allowing more outward travel
        # than the previous conservative clamp.
        dimension_clearance = 90.0
        label_edge_reserve = 72.0
        left_x = max(
            label_edge_reserve,
            min(self.width()*.30, left_edge-dimension_clearance))
        right_x = min(
            self.width()-label_edge_reserve,
            max(self.width()*.70, right_edge+dimension_clearance))
        for first, second, text, color in sst_rows:
            self._draw_dimension(
                painter, first, second, left_x, text, color, "left")
        for first, second, text, color in seed_rows:
            self._draw_dimension(
                painter, first, second, right_x, text, color, "right")

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), self._preview_background)

        # These coordinates are the same absolute centred ranges serialized
        # into the Square caDNAno JSON.  Preview and output therefore cannot
        # drift when Z1/Z2/Z3 changes.
        sst_first = tuple(map(float, self._sst_layer_ranges[0]))
        sst_second = tuple(map(float, self._sst_layer_ranges[1]))
        seed_first = tuple(map(float, self._seed_partition_ranges[0]))
        seed_spacing = tuple(map(float, self._seed_partition_ranges[1]))
        seed_second = tuple(map(float, self._seed_partition_ranges[2]))
        axial_low = min(sst_first[0], seed_first[0])
        axial_high = max(sst_second[1], seed_second[1])
        total = max(1.0, axial_high - axial_low)
        height = total * 0.34
        width = max(7*self._seed_helix_spacing + 2.0,
                    2*self.LATTICE_HALF_EXTENT_NM + 2.0)
        content_height = max(
            1.0, self.height()-PREVIEW_CONTENT_TOP-PREVIEW_CONTENT_BOTTOM)
        self._world_scale = min(
            self.width()/max(30.0, width*1.7),
            content_height/max(
                30.0, height*PREVIEW_MODEL_HEIGHT_FACTOR))
        center_bp = (axial_low + axial_high) / 2.0
        to_world = lambda value: (value - center_bp) * 0.34
        sst_z0, sst_z1_end = map(to_world, sst_first)
        sst_z2_start, sst_z3 = map(to_world, sst_second)
        seed_z0, seed_z1_end = map(to_world, seed_first)
        seed_z2_start, seed_z3 = map(to_world, seed_second)
        visual_twist = _right_handed_preview_rotation(self._angle)
        # SST volumes and Seed rods are deliberately independent so overhang
        # and complete-support margins are visible immediately.
        self._draw_lattice_volume(
            painter, sst_z0, sst_z1_end, 0.0, QColor("#2a78d1"),
            self._layer_lattices[0], self._layer_constants[0])
        self._draw_lattice_volume(
            painter, sst_z2_start, sst_z3,
            visual_twist, QColor("#d65b74"),
            self._layer_lattices[1], self._layer_constants[1])

        coordinates = [
            ((col-3.5)*self._seed_helix_spacing,
             (row-3.5)*self._seed_helix_spacing)
            for row, col in sorted(self._seed_cells)]
        rods = []
        for x, y in coordinates:
            lower_points = [(x, y, seed_z0), (x, y, seed_z1_end)]
            spacer_points = [(x, y, seed_z1_end)]
            spacer_steps = max(
                1, min(16, int(self._spacing // 8)
                       if self._spacing else 1))
            for index in range(1, spacer_steps+1):
                fraction = index/spacer_steps
                rotation = math.radians(visual_twist*fraction)
                current = (
                    math.cos(rotation)*x-math.sin(rotation)*y,
                    math.sin(rotation)*x+math.cos(rotation)*y,
                    seed_z1_end +
                    (to_world(seed_spacing[1])-seed_z1_end)*fraction)
                spacer_points.append(current)
            rotation = math.radians(visual_twist)
            upper_x = math.cos(rotation)*x-math.sin(rotation)*y
            upper_y = math.sin(rotation)*x+math.cos(rotation)*y
            upper_points = [(upper_x, upper_y, seed_z2_start),
                            (upper_x, upper_y, seed_z3)]
            all_points = lower_points+spacer_points+upper_points
            mean_depth = sum(_rotate(point, self._yaw, self._pitch)[1]
                             for point in all_points)/len(all_points)
            rods.append((mean_depth, lower_points, spacer_points, upper_points))
        scale = self._world_scale*self._zoom
        for unused_depth, lower, spacer, upper in sorted(
                rods, key=lambda item: item[0]):
            self._draw_rod_path(
                painter, lower, scale, "#4f8fce", "#c8e2fa")
            self._draw_rod_path(
                painter, spacer, scale, "#d9dee3", "#f4f6f8")
            self._draw_rod_path(
                painter, upper, scale, "#d96a82", "#f8ced8")

        # Dimension anchors are projected from the current model pose on
        # every frame.  Labels are then painted in screen space, so they stay
        # horizontal while following rotation, pan and zoom.
        self._draw_axial_dimensions(
            painter, sst_first, sst_second,
            seed_first, seed_spacing, seed_second,
            to_world, visual_twist)

        _draw_model_coverage_frame(painter, self.width(), self.height())

        scale_pixels_per_nm = self._world_scale*self._zoom
        self._last_scale_bar_px = _draw_scale_bar(
            painter, self.width(), self.height(),
            scale_pixels_per_nm)


class MoireTopViewPreview(QWidget):
    """Interactive top view of Square/Kagome bilayer lattices."""

    GRID_COUNT = 100
    HELIX_DIAMETER_NM = 2.0
    PIXELS_PER_CELL = 10
    FIELD_SPAN_NM = GRID_COUNT*2.8

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(360, 260)
        self._preview_background = QColor(PREVIEW_BACKGROUND)
        self._twist_angle = 3.2967555036483183
        self._period_nm = None
        self._lattice_constant = 2.8
        self._layer_lattices = ("square", "square")
        self._layer_constants = (2.8, 2.8)
        self._symmetry_label = "Square–Square"
        self._seed_cells = _default_seed_cells()
        self._seed_helix_spacing = 2.8
        self._view_rotation = 0.0
        self._zoom = 1.0
        self._pan = QPointF(0.0, 0.0)
        self._last = None
        self._drag_button = None
        self._grid_cache_key = None
        self._blue_lattice = QPixmap()
        self._red_lattice = QPixmap()
        self._last_scale_bar_px = 0.0
        self.setToolTip(
            "100×100根helix俯视图；左键拖动旋转；右键拖动平移；"
            "滚轮缩放；双击复位")

    def set_preview_background(self, color):
        self._preview_background = QColor(color)
        self._grid_cache_key = None
        self.update()

    def set_design(self, project):
        self._twist_angle = float(project.prediction["reported_angle_deg"])
        self._period_nm = project.prediction.get(
            "predicted_moire_period_nm")
        self._lattice_constant = float(project.settings.lattice_constant_nm)
        self._layer_lattices = tuple(project.prediction.get(
            "layer_lattice_types", ("square", "square")))
        self._layer_constants = tuple(project.prediction.get(
            "layer_lattice_constants_nm", (
                self._lattice_constant, self._lattice_constant)))
        self._symmetry_label = {
            "square_square_c4": "Square–Square",
            "kagome_kagome": "Kagome–Kagome",
            "square_kagome": "Square–Kagome",
        }.get(project.settings.lattice_symmetry,
              project.settings.lattice_symmetry)
        self._seed_cells = {tuple(map(int, cell)) for cell in
                            project.seed_plan.get(
                                "cross_section_cells",
                                _default_seed_cells())}
        self._grid_cache_key = None
        self.update()

    def set_configuration(self, symmetry, seed_cells, layer_constants,
                          twist_angle=3.2967555036483183, period_nm=None):
        self._layer_lattices = {
            "square_square_c4": ("square", "square"),
            "kagome_kagome": ("kagome", "kagome"),
            "square_kagome": ("square", "kagome"),
        }.get(symmetry, ("square", "square"))
        self._symmetry_label = {
            "square_square_c4": "Square–Square",
            "kagome_kagome": "Kagome–Kagome",
            "square_kagome": "Square–Kagome",
        }.get(symmetry, str(symmetry))
        self._layer_constants = tuple(map(float, layer_constants))
        self._lattice_constant = self._layer_constants[0]
        self._twist_angle = float(twist_angle)
        self._period_nm = period_nm
        self._seed_cells = {tuple(map(int, cell)) for cell in seed_cells}
        self._grid_cache_key = None
        self.update()

    def mousePressEvent(self, event):
        if event.button() in (Qt.MouseButton.LeftButton,
                              Qt.MouseButton.RightButton):
            self._last = event.position()
            self._drag_button = event.button()

    def mouseMoveEvent(self, event):
        if self._last is None:
            return
        delta = event.position()-self._last
        self._last = event.position()
        if self._drag_button == Qt.MouseButton.RightButton:
            self._pan += delta
        else:
            self._view_rotation += (delta.x()+delta.y()*.25)*.28
        self.update()

    def mouseReleaseEvent(self, event):
        self._last = None
        self._drag_button = None

    def mouseDoubleClickEvent(self, event):
        self._view_rotation = 0.0
        self._zoom = 1.0
        self._pan = QPointF(0.0, 0.0)
        self.update()

    def wheelEvent(self, event):
        self._zoom *= 1.12 if event.angleDelta().y() > 0 else 1/1.12
        self._zoom = max(.35, min(8.0, self._zoom))
        self.update()

    def _make_lattice_pixmap(self, color, lattice_type, lattice_constant):
        size = self.GRID_COUNT*self.PIXELS_PER_CELL
        lattice = QPixmap(size, size)
        lattice.fill(Qt.GlobalColor.transparent)
        painter = QPainter(lattice)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        points, edges = lattice_graph(
            lattice_type, lattice_constant, self.FIELD_SPAN_NM/2.0)
        pixels_per_nm = size/self.FIELD_SPAN_NM

        def mapped(point):
            return QPointF(
                (point[0]+self.FIELD_SPAN_NM/2.0)*pixels_per_nm,
                (point[1]+self.FIELD_SPAN_NM/2.0)*pixels_per_nm)

        line_color = QColor(color)
        line_color.setAlpha(max(105, min(175, line_color.alpha())))
        line_pen = QPen(line_color, 1.15)
        line_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(line_pen)
        for first, second in edges:
            painter.drawLine(mapped(points[first]), mapped(points[second]))
        node_color = QColor(color)
        node_color.setAlpha(max(175, node_color.alpha()))
        radius = self.HELIX_DIAMETER_NM*pixels_per_nm/2.0
        painter.setPen(QPen(QColor(5, 10, 15, 205), .75))
        painter.setBrush(node_color)
        for point in points:
            painter.drawEllipse(mapped(point), radius, radius)
        painter.end()
        return lattice

    def _ensure_grid_pixmaps(self):
        key = (self._layer_lattices,
               tuple(round(value, 6) for value in self._layer_constants),
               self.GRID_COUNT, self.HELIX_DIAMETER_NM)
        if key == self._grid_cache_key:
            return
        self._grid_cache_key = key
        self._blue_lattice = self._make_lattice_pixmap(
            QColor(42, 120, 209, 145), self._layer_lattices[0],
            self._layer_constants[0])
        self._red_lattice = self._make_lattice_pixmap(
            QColor(214, 91, 116, 145), self._layer_lattices[1],
            self._layer_constants[1])

    def _draw_lattice(self, painter, relative_angle, lattice,
                      lattice_constant):
        painter.save()
        painter.rotate(relative_angle)
        width_nm = self.FIELD_SPAN_NM
        target = QRectF(-width_nm/2.0, -width_nm/2.0,
                        width_nm, width_nm)
        painter.drawPixmap(target, lattice, QRectF(lattice.rect()))
        painter.restore()

    def _draw_seed(self, painter):
        """Overlay the selected Square-grid Seed cross-section in white."""
        radius = self.HELIX_DIAMETER_NM/2.0
        outline = QPen(QColor(28, 38, 48, 230))
        outline.setWidthF(.8)
        outline.setCosmetic(True)
        painter.setPen(outline)
        painter.setBrush(QColor(255, 255, 255, 245))
        for row, column in sorted(self._seed_cells):
            x = (column-3.5)*self._seed_helix_spacing
            y = (row-3.5)*self._seed_helix_spacing
            painter.drawEllipse(QPointF(x, y), radius, radius)

    def _clear_seed_footprint(self, painter):
        """Remove SST dots from the complete 8 x 8 Seed-S footprint."""
        # Include one SST helix radius so a node centred just outside the
        # nominal 8 x 8 cell boundary cannot remain underneath an outer Seed
        # helix after antialiasing.
        side = (8.0*self._seed_helix_spacing+self.HELIX_DIAMETER_NM)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(self._preview_background)
        painter.drawRect(QRectF(-side/2.0, -side/2.0, side, side))

    def _map_world_to_screen(self, point, scale):
        """Apply the live top-view pan/rotation/scale to one world point."""
        angle = math.radians(self._view_rotation)
        cosine, sine = math.cos(angle), math.sin(angle)
        if isinstance(point, QPointF):
            x, y = point.x(), point.y()
        else:
            x, y = point
        x, y = cosine*x-sine*y, sine*x+cosine*y
        top_margin = PREVIEW_CONTENT_TOP
        usable_height = max(
            1.0, self.height()-top_margin-PREVIEW_CONTENT_BOTTOM)
        return QPointF(
            self.width()/2.0+self._pan.x()+x*scale,
            top_margin+usable_height/2.0+self._pan.y()+y*scale)

    def _draw_twist_annotation(self, painter, scale):
        """Trace both complete top edges and label their included angle."""
        visual_angle = _right_handed_preview_rotation(self._twist_angle)
        angle = math.radians(visual_angle)
        half = self.FIELD_SPAN_NM/2.0
        # Intersection of the first layer's top edge and the rotated second
        # layer's top edge.  -h*tan(a/2) is stable at small angles.
        intersection_x = -half*math.tan(angle/2.0)
        intersection_x = max(-half*.78, min(half*.78, intersection_x))
        origin_world = QPointF(intersection_x, -half)
        origin = self._map_world_to_screen(origin_world, scale)

        # The guides use the actual complete top boundary of each finite
        # lattice, rather than two short decorative rays.
        first_start = self._map_world_to_screen(
            QPointF(-half, -half), scale)
        first_end = self._map_world_to_screen(
            QPointF(half, -half), scale)
        local_start = QPointF(-half, -half)
        local_end = QPointF(half, -half)
        second_start = self._map_world_to_screen(QPointF(
            math.cos(angle)*local_start.x()-math.sin(angle)*local_start.y(),
            math.sin(angle)*local_start.x()+math.cos(angle)*local_start.y()),
            scale)
        second_end = self._map_world_to_screen(QPointF(
            math.cos(angle)*local_end.x()-math.sin(angle)*local_end.y(),
            math.sin(angle)*local_end.x()+math.cos(angle)*local_end.y()),
            scale)

        painter.save()
        ray_pen = QPen(QColor(222, 231, 238, 220), 1.15,
                       Qt.PenStyle.DashLine)
        ray_pen.setCosmetic(True)
        painter.setPen(ray_pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawLine(first_start, first_end)
        painter.drawLine(second_start, second_end)

        # A short angular arc degenerates into a spike at experimentally
        # common small twists.  The two complete lattice-edge guides already
        # show the physical angle, so label their true intersection directly.
        painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
        painter.setPen(QColor("#ffffff"))
        if self._twist_angle > 1e-9:
            handedness = "right-handed"
        elif self._twist_angle < -1e-9:
            handedness = "left-handed"
        else:
            handedness = "untwisted"
        label_text = "Twist: %+.1f° (%s)" % (
            self._twist_angle, handedness)
        label_width = painter.fontMetrics().horizontalAdvance(
            label_text)+8.0
        label_rect = QRectF(
            max(4.0, min(self.width()-label_width-4.0,
                         origin.x()-label_width/2.0)),
            max(3.0, min(self.height()-52.0, origin.y()-28.0)),
            label_width, 20.0)
        painter.drawText(
            label_rect, Qt.AlignmentFlag.AlignCenter,
            label_text)
        painter.restore()

    def _draw_period_annotation(self, painter, scale):
        """Outline one physical Moiré unit and label its period."""
        if self._period_nm is None:
            return
        period = float(self._period_nm)
        if not math.isfinite(period) or period <= 0.0:
            return
        # Square uses orthogonal basis vectors; Kagome uses its 60-degree
        # triangular Bravais basis.  The displayed unit must start at a real
        # repeated coincidence point, but it must not use the Seed-centre
        # coincidence point as one of its corners.
        base_angle = math.radians(
            _right_handed_preview_rotation(self._twist_angle)*.5)
        if self._layer_lattices == ("kagome", "kagome"):
            base_angle += math.radians(30.0)
            included_angle = math.radians(60.0)
        else:
            included_angle = math.radians(90.0)
        first_vector = QPointF(
            period*math.cos(base_angle), period*math.sin(base_angle))
        second_angle = base_angle+included_angle
        second_vector = QPointF(
            period*math.cos(second_angle), period*math.sin(second_angle))
        # Select the nearest translated primitive cell whose bounding box is
        # fully outside the Seed footprint.  Translation is restricted to
        # integer combinations of the two Moiré basis vectors, so the outline
        # remains a genuine repeat unit rather than an arbitrarily moved box.
        seed_half = (
            8.0*self._seed_helix_spacing+self.HELIX_DIAMETER_NM)/2.0
        search_radius = max(3, int(math.ceil(
            (2.0*seed_half)/max(period, 1e-6)))+3)
        candidates = []
        for first_step in range(-search_radius, search_radius+1):
            for second_step in range(-search_radius, search_radius+1):
                anchor = QPointF(
                    first_step*first_vector.x()+
                    second_step*second_vector.x(),
                    first_step*first_vector.y()+
                    second_step*second_vector.y())
                vertices = (
                    anchor,
                    anchor+first_vector,
                    anchor+first_vector+second_vector,
                    anchor+second_vector)
                xs = [point.x() for point in vertices]
                ys = [point.y() for point in vertices]
                outside_seed = (
                    min(xs) > seed_half or max(xs) < -seed_half or
                    min(ys) > seed_half or max(ys) < -seed_half)
                if not outside_seed:
                    continue
                center_x = sum(xs)/4.0
                center_y = sum(ys)/4.0
                # Prefer the lower-right equivalent only as a deterministic
                # tie-breaker between equally near physical repeat cells.
                quadrant_penalty = 0 if center_x >= 0 and center_y >= 0 else 1
                candidates.append((
                    round(math.hypot(center_x, center_y), 9),
                    quadrant_penalty, abs(center_y), abs(center_x), anchor))
        coincidence = min(candidates, key=lambda item: item[:4])[-1] \
            if candidates else QPointF(first_vector+second_vector)
        polygon_world = [
            coincidence,
            coincidence+first_vector,
            coincidence+first_vector+second_vector,
            coincidence+second_vector,
        ]
        polygon = QPolygonF([
            self._map_world_to_screen(point, scale)
            for point in polygon_world])
        painter.save()
        period_pen = QPen(QColor("#ffffff"), 1.25,
                          Qt.PenStyle.DashLine)
        period_pen.setCosmetic(True)
        painter.setPen(period_pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPolygon(polygon)
        bounds = polygon.boundingRect()
        painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
        painter.setPen(QColor("#ffffff"))
        label_text = "Moiré period: %.1f nm" % period
        label_width = painter.fontMetrics().horizontalAdvance(label_text)+4.0
        label_rect = QRectF(
            max(4.0, min(self.width()-label_width-4.0,
                         bounds.center().x()-label_width/2.0)),
            max(4.0, min(self.height()-56.0, bounds.top()-24.0)),
            label_width, 20.0)
        painter.drawText(
            label_rect,
            Qt.AlignmentFlag.AlignCenter,
            label_text)
        painter.restore()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        painter.fillRect(self.rect(), self._preview_background)
        self._ensure_grid_pixmaps()
        top_margin = PREVIEW_CONTENT_TOP
        bottom_margin = PREVIEW_CONTENT_BOTTOM
        usable_width = max(1.0, self.width()-34.0)
        usable_height = max(1.0, self.height()-top_margin-bottom_margin)
        size_nm = self.FIELD_SPAN_NM
        # Match the default apparent model height of BilayerPreview instead
        # of filling the entire 2D canvas. Zoom remains user-controlled.
        target_height = usable_height/PREVIEW_MODEL_HEIGHT_FACTOR
        scale = min(usable_width, target_height)/size_nm*self._zoom
        painter.save()
        painter.translate(self.width()/2+self._pan.x(),
                          top_margin+usable_height/2.0+self._pan.y())
        painter.rotate(self._view_rotation)
        painter.scale(scale, scale)
        self._draw_lattice(
            painter, 0.0, self._blue_lattice, self._layer_constants[0])
        self._draw_lattice(
            painter, _right_handed_preview_rotation(self._twist_angle),
            self._red_lattice,
            self._layer_constants[1])
        self._clear_seed_footprint(painter)
        self._draw_seed(painter)
        painter.restore()

        self._draw_twist_annotation(painter, scale)
        self._draw_period_annotation(painter, scale)

        _draw_model_coverage_frame(painter, self.width(), self.height())

        self._last_scale_bar_px = _draw_scale_bar(
            painter, self.width(), self.height(), scale)


class CaptureMapPreview(QWidget):
    """Schematic fully cooperative north/south capture map."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(300)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#fbfcfe"))
        size = min(self.width()-100, self.height()-70)
        cell = size/8.0
        left = (self.width()-size)/2.0
        top = 38.0
        active = QColor("#6b61bd")
        inactive = QColor("#d6dce5")
        edge = QColor("#536273")
        for row in range(8):
            for col in range(8):
                if 2 <= row <= 5 and 2 <= col <= 5:
                    continue
                rect = QRectF(left+col*cell+2, top+row*cell+2,
                              cell-4, cell-4)
                color = active if row in (0, 7) else inactive
                painter.setBrush(color)
                painter.setPen(QPen(edge, 1))
                painter.drawEllipse(rect)
                if row in (0, 7) and col % 2 == 0:
                    painter.setBrush(QColor("#f4c95d"))
                    painter.drawEllipse(rect.center(), 3.5, 3.5)
        painter.setPen(QColor("#233143"))
        painter.setFont(QFont("Arial", 11, QFont.Weight.DemiBold))
        painter.drawText(18, 24, "Fully cooperative capture map")
        painter.setFont(QFont("Arial", 11))
        painter.drawText(18, self.height()-18,
            "Purple: active north/south surfaces · Gold: capture-0/1 pair groups")


class StructureDesignPreview(QWidget):
    """Interactive cadnano-like Slice/Path preview for structure workflow."""

    EMPTY = [-1, -1, -1, -1]
    PATH_IMAGE_WIDTH = 3200
    PATH_IMAGE_HEIGHT = 3600
    SEQUENCE_SCAFFOLD_COLORS = (
        "#1769aa", "#d1495b", "#2a9d8f", "#7b61b8",
        "#e17c05", "#348aa7", "#8f5d2f", "#5c946e")

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(640, 480)
        self.setMouseTracking(True)
        self._payload = None
        self._rows = {}
        self._filename = None
        self._channel = None
        self._channel_title = "No structure design available"
        self._sst_helices = set()
        self._seed_helices = set()
        self._capture_seed = {}
        self._capture_seed_colors = {}
        self._capture_sst = set()
        self._capture_seed_ends = {}
        self._capture_faces = []
        self._staple_colors = {}
        self._sequence_scaffold_colors = {}
        self._sequence_scaffold_nodes = set()
        self._sequence_base_labels = {}
        self._sst_input_nodes = set()
        self._slice_zoom = 1.0
        self._slice_pan = QPointF()
        self._face_zoom = 1.0
        self._face_pan = QPointF()
        self._path_zoom = 1.0
        self._path_pan = QPointF()
        self._drag_target = None
        self._last_position = None
        self._slice_rect = QRectF()
        self._capture_rect = QRectF()
        self._path_rect = QRectF()
        self._path_outer_rect = QRectF()
        self._path_report_rect = QRectF()
        self._path_report_handle = QRectF()
        self._path_report_text = ""
        self._path_report_ratio = .28
        self._path_pixmap = None
        self._panel_handles = []
        # Slice and Capture share the upper row; Path uses the full lower row.
        self._complete_panel_boundaries = [.36, .54]
        self._simple_panel_boundary = .38
        self.setToolTip(
            "Pan and zoom the Slice, Capture side, and Path views "
            "independently. The Path view uses a pre-rendered high-resolution "
            "image. Squares and triangles indicate 5′ and 3′ ends, "
            "respectively. Drag the panel dividers vertically; drag the Path "
            "report divider horizontally to collapse it. Double-click to "
            "reset a view.")

    def clear(self, message="No structure design is available for display."):
        self._payload = None
        self._rows = {}
        self._filename = None
        self._channel = None
        self._channel_title = message
        self._capture_seed = {}
        self._capture_seed_colors = {}
        self._capture_sst = set()
        self._capture_seed_ends = {}
        self._capture_faces = []
        self._staple_colors = {}
        self._sequence_scaffold_colors = {}
        self._sequence_scaffold_nodes = set()
        self._sequence_base_labels = {}
        self._sst_input_nodes = set()
        self._path_report_text = ""
        self._invalidate_path_pixmap()
        self.setMinimumHeight(480)
        self.update()

    def set_path_report(self, text):
        """Show a collapsible report immediately to the left of Path."""
        self._path_report_text = str(text or "").strip()
        self._layout_panels()
        self.update()

    def set_source(self, filename, channel, title=None):
        path = Path(filename).expanduser().resolve()
        payload = json.loads(path.read_text(encoding="utf-8"))
        self._payload = payload
        self._rows = {
            int(row["num"]): row for row in payload.get("vstrands", [])}
        metadata = payload.get("moire_structure_metadata", {})
        sst = metadata.get("sst_helix_numbers")
        seed = metadata.get("seed_helix_numbers")
        if sst is None or seed is None:
            sst_first = metadata.get("helix_numbering") == "sst_first"
            if sst_first:
                sst, seed = range(16), range(16, 64)
            else:
                seed, sst = range(48), range(48, 64)
        self._sst_helices = {int(number) for number in sst}
        self._seed_helices = {int(number) for number in seed}
        sst_first = metadata.get("helix_numbering") == "sst_first"
        face_definitions = metadata.get("capture_face_definitions") or [
            {
                "id": "face1",
                "label": "Face 1 · upper Seed edge",
                "internal_seed_helices": list(range(0, 8)),
                "sst_first_seed_helices": list(range(16, 24)),
                "color": "#7b61b8",
            },
            {
                "id": "face2",
                "label": "Face 2 · lower Seed edge",
                "internal_seed_helices": list(range(24, 32)),
                "sst_first_seed_helices": list(range(40, 48)),
                "color": "#2a9d8f",
            },
        ]
        helix_key = ("sst_first_seed_helices" if sst_first else
                     "internal_seed_helices")
        self._capture_faces = [{
            "id": item.get("id", "face%d" % (index + 1)),
            "label": item.get(
                "label", "Face %d" % (index + 1)),
            "helices": [int(number) for number in item.get(helix_key, ())],
            "color": item.get(
                "color", "#7b61b8" if index == 0 else "#2a9d8f"),
        } for index, item in enumerate(face_definitions[:2])]
        self._filename = str(path)
        self._channel = str(channel)
        self._channel_title = title or path.name
        self._capture_seed, self._capture_sst = self._capture_sites()
        self._capture_seed_ends = self._capture_end_labels()
        self._staple_colors = self._staple_component_colors()
        self._capture_seed, self._capture_seed_colors = \
            self._capture_preview_sites(
                metadata, face_definitions, self._capture_seed)
        self._sequence_scaffold_colors = self._scaffold_component_colors()
        self._sequence_scaffold_nodes = set()
        self._sequence_base_labels = {}
        self._sst_input_nodes = set()
        if self._channel == "complete":
            self.setMinimumHeight(560)
        elif self._channel in ("sequence_scaffold", "sst_input"):
            self.setMinimumHeight(440)
        else:
            self.setMinimumHeight(480)
        self._slice_zoom = 1.0
        self._slice_pan = QPointF()
        self._face_zoom = 1.0
        self._face_pan = QPointF()
        self._path_zoom = 1.0
        self._path_pan = QPointF()
        self._invalidate_path_pixmap()
        self.update()

    def _invalidate_path_pixmap(self):
        self._path_pixmap = None

    def _layout_panels(self):
        """Lay out resizable Slice/Capture/Path panels and drag handles."""
        left = 14.0
        width = max(1.0, self.width() - 28.0)
        # The outer group box already identifies the preview.  Begin the
        # actual panels immediately below it instead of reserving a second,
        # redundant channel-title row.
        top = 10.0
        # The former mouse-operation footer has been removed.  Let the
        # preview panels use the reclaimed vertical space.
        bottom = max(top + 1.0, self.height() - 14.0)
        span = max(1.0, bottom - top)
        gap = 12.0
        self._panel_handles = []
        if self._channel in ("sequence_scaffold", "sst_input"):
            self._slice_rect = QRectF()
            self._capture_rect = QRectF()
            self._path_rect = QRectF(left, top, width, span)
            self._layout_path_report()
            return
        if self._channel == "complete":
            first = top + span * self._complete_panel_boundaries[0]
            column_gap = 12.0
            column_width = max(1.0, (width - column_gap) / 2.0)
            self._slice_rect = QRectF(
                left, top, column_width,
                max(1.0, first - top - gap / 2.0))
            self._capture_rect = QRectF(
                left + column_width + column_gap, top, column_width,
                max(1.0, first - top - gap / 2.0))
            self._path_rect = QRectF(
                left, first + gap / 2.0, width,
                max(1.0, bottom - first - gap / 2.0))
            self._panel_handles = [
                QRectF(left, first - gap / 2.0, width, gap),
            ]
        else:
            boundary = top + span * self._simple_panel_boundary
            self._slice_rect = QRectF(
                left, top, width,
                max(1.0, boundary - top - gap / 2.0))
            self._capture_rect = QRectF()
            self._path_rect = QRectF(
                left, boundary + gap / 2.0, width,
                max(1.0, bottom - boundary - gap / 2.0))
            self._panel_handles = [
                QRectF(left, boundary - gap / 2.0, width, gap)]
        self._layout_path_report()

    def _layout_path_report(self):
        """Split the Path row into a report column and the Path viewport."""
        self._path_outer_rect = QRectF(self._path_rect)
        self._path_report_rect = QRectF()
        self._path_report_handle = QRectF()
        if not self._path_report_text or self._path_outer_rect.isEmpty():
            return
        outer = self._path_outer_rect
        handle_width = 12.0
        ratio = max(0.0, min(.56, float(self._path_report_ratio)))
        report_width = outer.width() * ratio
        handle_x = (outer.left() + handle_width / 2.0
                    if report_width < 2.0 else
                    outer.left() + report_width)
        if report_width >= 2.0:
            self._path_report_rect = QRectF(
                outer.left(), outer.top(),
                max(1.0, report_width - handle_width / 2.0),
                outer.height())
        self._path_report_handle = QRectF(
            handle_x - handle_width / 2.0, outer.top(),
            handle_width, outer.height())
        path_left = handle_x + handle_width / 2.0
        self._path_rect = QRectF(
            path_left, outer.top(),
            max(1.0, outer.right() - path_left), outer.height())

    def _draw_panel_handles(self, painter):
        painter.save()
        for handle in self._panel_handles:
            painter.setPen(QPen(QColor("#c0cad5"), 1.0))
            painter.setBrush(QColor("#e5ebf1"))
            grip = QRectF(handle.center().x() - 24.0,
                          handle.center().y() - 2.0, 48.0, 4.0)
            painter.drawRoundedRect(grip, 2.0, 2.0)
        if not self._path_report_handle.isEmpty():
            handle = self._path_report_handle
            painter.setPen(QPen(QColor("#b9c4d0"), 1.0))
            painter.setBrush(QColor("#e2e8ef"))
            painter.drawRoundedRect(handle, 3.0, 3.0)
            grip = QRectF(handle.center().x() - 2.0,
                          handle.center().y() - 24.0, 4.0, 48.0)
            painter.setBrush(QColor("#9eabb9"))
            painter.drawRoundedRect(grip, 2.0, 2.0)
        painter.restore()

    def _draw_path_report(self, painter):
        rect = self._path_report_rect
        if rect.isEmpty() or not self._path_report_text:
            return
        painter.save()
        painter.setClipRect(rect)
        painter.setPen(QPen(QColor("#e5c769"), 1.1))
        painter.setBrush(QColor("#fff7d8"))
        painter.drawRoundedRect(rect.adjusted(.5, .5, -.5, -.5), 8, 8)
        painter.setPen(QColor("#705820"))
        painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
        painter.drawText(
            rect.adjusted(12, 10, -10, -8),
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop |
            Qt.TextFlag.TextWordWrap,
            self._path_report_text)
        painter.restore()

    def _panel_handle_at(self, position):
        for index, handle in enumerate(self._panel_handles):
            if handle.contains(position):
                return index
        return None

    @staticmethod
    def _capture_protrusion_side(number):
        """Return the actual cadnano staple side used by a capture bridge."""
        # In cadnano's Path view staple runs on the lower side of an even
        # helix and the upper side of an odd helix.  The capture bridge exits
        # from that same side; adjacent helices therefore form the expected
        # periodic up/down pattern without inventing idealized coordinates.
        return 1.0 if int(number) % 2 == 0 else -1.0

    @staticmethod
    def _occupied(record):
        return list(record) != StructureDesignPreview.EMPTY

    @classmethod
    def _records_linked(cls, number, records, first_index, second_index):
        """True only when adjacent occupied bases share a backbone edge."""
        if (first_index < 0 or second_index < 0 or
                first_index >= len(records) or second_index >= len(records) or
                not cls._occupied(records[first_index]) or
                not cls._occupied(records[second_index])):
            return False
        first = records[first_index]
        second = records[second_index]
        first_neighbors = {(int(first[0]), int(first[1])),
                           (int(first[2]), int(first[3]))}
        second_neighbors = {(int(second[0]), int(second[1])),
                            (int(second[2]), int(second[3]))}
        return ((int(number), int(second_index)) in first_neighbors or
                (int(number), int(first_index)) in second_neighbors)

    @classmethod
    def _segments(cls, records, number=None):
        """Return occupied runs without visually bridging a real nick."""
        segments = []
        start = None
        for index, record in enumerate(records):
            occupied = cls._occupied(record)
            if occupied and start is None:
                start = index
            elif (occupied and start is not None and number is not None and
                  not cls._records_linked(number, records, index - 1, index)):
                segments.append((start, index - 1))
                start = index
            elif not occupied and start is not None:
                segments.append((start, index - 1))
                start = None
        if start is not None:
            segments.append((start, len(records) - 1))
        return segments

    @classmethod
    def _strand_special_nodes(cls, number, records):
        """Bases carrying an endpoint/nick or leaving via a crossover."""
        special = set()
        for index, record in enumerate(records):
            if not cls._occupied(record):
                continue
            for offset in (0, 2):
                partner = int(record[offset])
                if partner < 0 or partner != int(number):
                    special.add(index)
        return special

    def _capture_sites(self):
        seed_sites = {}
        sst_sites = set()
        for number, row in self._rows.items():
            for index, record in enumerate(row.get("stap", [])):
                for offset in (0, 2):
                    partner, partner_index = map(
                        int, record[offset:offset + 2])
                    if partner < 0:
                        continue
                    if number in self._seed_helices and \
                            partner in self._sst_helices:
                        seed_sites.setdefault(number, set()).add(index)
                        sst_sites.add((partner, partner_index))
                    elif number in self._sst_helices and \
                            partner in self._seed_helices:
                        sst_sites.add((number, index))
                        seed_sites.setdefault(partner, set()).add(
                            partner_index)
        return seed_sites, sst_sites

    def _capture_preview_sites(self, metadata, face_definitions,
                               actual_sites):
        """Include every physical and export-only capture-helix site."""
        sites = {int(number): set(positions)
                 for number, positions in actual_sites.items()}
        colors = {
            (int(number), int(position)): self._staple_colors.get(
                (int(number), int(position)), 0x000000)
            for number, positions in sites.items()
            for position in positions}
        assignments = (
            metadata.get("capture_export_site_assignments_internal") or
            metadata.get("capture_site_assignments_internal") or [])
        if not assignments:
            return sites, colors

        number_map = {}
        sst_first = metadata.get("helix_numbering") == "sst_first"
        for face in face_definitions:
            internal = [int(value) for value in
                        face.get("internal_seed_helices", ())]
            displayed = [int(value) for value in face.get(
                "sst_first_seed_helices" if sst_first else
                "internal_seed_helices", ())]
            number_map.update(zip(internal, displayed))

        pair_by_position = {}
        pair_offset = 0
        for layer_positions in metadata.get(
                "capture_positions_by_layer", ()):
            layer_positions = [int(value) for value in layer_positions]
            for position_index, position in enumerate(layer_positions):
                pair_by_position[position] = pair_offset + position_index // 2
            pair_offset += (len(layer_positions) + 1) // 2
        pair_colors = [int(value) for value in
                       metadata.get("capture_pair_colors", ())]

        for assignment in assignments:
            position = int(assignment["position"])
            pair_index = pair_by_position.get(position)
            fallback = (pair_colors[pair_index]
                        if pair_index is not None and
                        pair_index < len(pair_colors) else 0x000000)
            for bridge in assignment.get("bridges", ()):
                internal_number = int(bridge["seed_helix"])
                number = number_map.get(internal_number, internal_number)
                if number not in self._seed_helices:
                    continue
                sites.setdefault(number, set()).add(position)
                colors[(number, position)] = self._staple_colors.get(
                    (number, position), fallback)
        return sites, colors

    def _capture_end_labels(self):
        """Return the actual staple side that exits Seed toward SST."""
        result = {}
        for number in self._seed_helices:
            row = self._rows.get(number, {})
            for index, record in enumerate(row.get("stap", [])):
                for offset, label in ((0, "5′"), (2, "3′")):
                    if int(record[offset]) in self._sst_helices:
                        result.setdefault(number, {})[index] = label
        return result

    def _staple_component_colors(self):
        nodes = {
            (number, index)
            for number, row in self._rows.items()
            for index, record in enumerate(row.get("stap", []))
            if self._occupied(record)}
        neighbors = {node: [] for node in nodes}
        markers = {}
        for number, row in self._rows.items():
            for index, color in row.get("stap_colors", []):
                markers[(number, int(index))] = int(color)
        for number, index in nodes:
            record = self._rows[number]["stap"][index]
            for offset in (0, 2):
                partner = (int(record[offset]), int(record[offset + 1]))
                if partner in nodes:
                    neighbors[(number, index)].append(partner)
        colors = {}
        visited = set()
        for first in sorted(nodes):
            if first in visited:
                continue
            component = set([first])
            stack = [first]
            visited.add(first)
            while stack:
                node = stack.pop()
                for other in neighbors[node]:
                    if other not in visited:
                        visited.add(other)
                        component.add(other)
                        stack.append(other)
            color = next(
                (markers[node] for node in sorted(component)
                 if node in markers), 0x000000)
            for node in component:
                colors[node] = color
        return colors

    def _scaffold_component_colors(self, markers=None):
        markers = markers or {}
        nodes = {
            (number, index)
            for number, row in self._rows.items()
            for index, record in enumerate(row.get("scaf", []))
            if self._occupied(record)}
        neighbors = {node: [] for node in nodes}
        for number, index in nodes:
            record = self._rows[number]["scaf"][index]
            for offset in (0, 2):
                partner = (int(record[offset]), int(record[offset + 1]))
                if partner in nodes:
                    neighbors[(number, index)].append(partner)
        colors = {}
        visited = set()
        component_index = 0
        for first in sorted(nodes):
            if first in visited:
                continue
            component = {first}
            stack = [first]
            visited.add(first)
            while stack:
                node = stack.pop()
                for other in neighbors[node]:
                    if other not in visited:
                        visited.add(other)
                        component.add(other)
                        stack.append(other)
            color = next(
                (markers[node] for node in sorted(component)
                 if node in markers),
                self.SEQUENCE_SCAFFOLD_COLORS[
                    component_index % len(self.SEQUENCE_SCAFFOLD_COLORS)])
            component_index += 1
            for node in component:
                colors[node] = color
        return colors

    def set_sequence_scaffold_targets(self, targets):
        """Use analysis-card colors for their corresponding Path components."""
        markers = {
            (int(item["start_vh"]), int(item["start_idx"])):
            item.get("color", "#1769aa")
            for item in targets
        }
        self._sequence_scaffold_nodes = self._connected_scaffold_nodes(
            set(markers))
        self._sequence_scaffold_colors = self._scaffold_component_colors(
            markers)
        self._invalidate_path_pixmap()
        self.update()

    def _connected_scaffold_nodes(self, starts):
        nodes = {
            (number, index)
            for number, row in self._rows.items()
            for index, record in enumerate(row.get("scaf", []))
            if self._occupied(record)
        }
        selected = set()
        stack = [node for node in starts if node in nodes]
        while stack:
            node = stack.pop()
            if node in selected:
                continue
            selected.add(node)
            number, index = node
            record = self._rows[number]["scaf"][index]
            for offset in (0, 2):
                partner = (int(record[offset]), int(record[offset + 1]))
                if partner in nodes and partner not in selected:
                    stack.append(partner)
        return selected

    def set_sequence_scaffold_assignments(self, assignments):
        """Overlay accepted scaffold bases along each 5'-to-3' routing."""
        labels = {}
        for assignment in assignments:
            if assignment.get("category") != "seed_scaffold":
                continue
            node = (int(assignment.get("start_vh", -1)),
                    int(assignment.get("start_idx", -1)))
            for base in str(assignment.get("sequence", "")).upper():
                if node not in self._sequence_scaffold_nodes:
                    break
                labels[node] = base
                number, index = node
                record = self._rows[number]["scaf"][index]
                node = (int(record[2]), int(record[3]))
        self._sequence_base_labels = labels
        self._invalidate_path_pixmap()
        self.update()

    def set_sst_input_assignments(self, assignments):
        """Overlay added SST input bases along each 5'-to-3' input path."""
        labels = {}
        for assignment in assignments:
            if not str(assignment.get("category", "")).startswith(
                    "sst_input_layer_"):
                continue
            node = (int(assignment.get("start_vh", -1)),
                    int(assignment.get("start_idx", -1)))
            for base in str(assignment.get("sequence", "")).upper():
                if node not in self._sst_input_nodes:
                    break
                labels[node] = base
                number, index = node
                record = self._rows[number]["scaf"][index]
                node = (int(record[2]), int(record[3]))
        self._sequence_base_labels = labels
        self._invalidate_path_pixmap()
        self.update()

    def set_sst_input_targets(self, targets):
        """Restrict the Path preview to the selected SST input oligos."""
        starts = {
            (int(item["start_vh"]), int(item["start_idx"]))
            for item in targets
        }
        self._sst_input_nodes = self._connected_scaffold_nodes(starts)
        self._invalidate_path_pixmap()
        self.update()

    def _path_draws_record(self, number, index, field):
        if self._channel == "sequence_scaffold":
            return (field == "scaf" and
                    (number, index) in self._sequence_scaffold_nodes)
        if self._channel == "sst_input":
            return field == "scaf" and (number, index) in self._sst_input_nodes
        if field == "stap":
            return self._path_shows_staple(number)
        return True

    def _path_records(self, number, field):
        records = self._rows[number].get(field, [])
        if self._channel not in ("sequence_scaffold", "sst_input"):
            return records
        return [
            record if self._path_draws_record(number, index, field)
            else self.EMPTY
            for index, record in enumerate(records)
        ]

    def _path_scaffold_color(self, number, index):
        if self._channel in ("sequence_scaffold", "scaffold", "complete"):
            return QColor(self._sequence_scaffold_colors.get(
                (number, index), "#1769aa"))
        return QColor("#0066cc")

    def _path_shows_staple(self, number):
        if self._channel in ("sequence_scaffold", "sst_input"):
            return False
        return (self._channel in ("sst", "complete") or
                number in self._sst_helices)

    def _visible_rows(self):
        if not self._rows:
            return []
        if self._channel in ("sst", "sst_input"):
            candidates = self._sst_helices
        else:
            candidates = set(self._rows)
        visible = []
        for number in sorted(candidates):
            row = self._rows.get(number)
            if row is None:
                continue
            fields = [self._path_records(number, "scaf")]
            if self._path_shows_staple(number):
                fields.append(self._path_records(number, "stap"))
            if any(any(self._occupied(record) for record in records)
                   for records in fields):
                visible.append(number)
        return visible

    def mousePressEvent(self, event):
        if event.button() not in (Qt.MouseButton.LeftButton,
                                  Qt.MouseButton.RightButton):
            return
        self._last_position = event.position()
        position = event.position()
        handle_index = self._panel_handle_at(position)
        if self._path_report_handle.contains(position):
            self._drag_target = "path_report_handle"
        elif handle_index is not None:
            self._drag_target = "panel_handle_%d" % handle_index
        elif self._path_report_rect.contains(position):
            self._drag_target = "path_report"
        elif self._slice_rect.contains(position):
            self._drag_target = "slice"
        elif self._capture_rect.contains(position):
            self._drag_target = "capture"
        else:
            self._drag_target = "path"

    def mouseMoveEvent(self, event):
        if self._last_position is None:
            return
        delta = event.position() - self._last_position
        self._last_position = event.position()
        path_only = False
        if self._drag_target == "path_report_handle":
            outer = self._path_outer_rect
            ratio = ((event.position().x() - outer.left()) /
                     max(1.0, outer.width()))
            # Snap fully shut near the left edge, while leaving the narrow
            # handle available so the report can always be reopened.
            self._path_report_ratio = 0.0 if ratio < .035 else max(
                .08, min(.56, ratio))
            self._layout_panels()
        elif str(self._drag_target).startswith("panel_handle_"):
            index = int(self._drag_target.rsplit("_", 1)[1])
            top = 42.0
            span = max(1.0, self.height() - 70.0)
            ratio = (event.position().y() - top) / span
            if self._channel == "complete":
                self._complete_panel_boundaries[0] = max(
                    .22, min(.64, ratio))
            else:
                self._simple_panel_boundary = max(.18, min(.82, ratio))
            self._layout_panels()
        elif self._drag_target == "slice":
            self._slice_pan += delta
        elif self._drag_target == "capture":
            self._face_pan += delta
        elif self._drag_target == "path_report":
            pass
        else:
            self._path_pan += delta
            path_only = True
        if path_only:
            self.update(self._path_rect.toAlignedRect())
        else:
            self.update()

    def mouseReleaseEvent(self, event):
        self._last_position = None
        self._drag_target = None

    def wheelEvent(self, event):
        if self._path_report_rect.contains(event.position()):
            event.accept()
            return
        factor = 1.12 if event.angleDelta().y() > 0 else 1 / 1.12
        path_only = False
        if self._slice_rect.contains(event.position()):
            self._slice_zoom = max(.5, min(4.0,
                self._slice_zoom * factor))
        elif self._capture_rect.contains(event.position()):
            self._face_zoom = max(.5, min(5.0,
                self._face_zoom * factor))
        else:
            self._path_zoom = max(.45, min(5.0,
                self._path_zoom * factor))
            path_only = True
        if path_only:
            self.update(self._path_rect.toAlignedRect())
        else:
            self.update()

    def mouseDoubleClickEvent(self, event):
        handle_index = self._panel_handle_at(event.position())
        path_only = False
        if self._path_report_handle.contains(event.position()):
            self._path_report_ratio = .28
            self._layout_panels()
        elif handle_index is not None:
            self._complete_panel_boundaries = [.36, .54]
            self._simple_panel_boundary = .38
            self._layout_panels()
        elif self._slice_rect.contains(event.position()):
            self._slice_zoom = 1.0
            self._slice_pan = QPointF()
        elif self._capture_rect.contains(event.position()):
            self._face_zoom = 1.0
            self._face_pan = QPointF()
        else:
            self._path_zoom = 1.0
            self._path_pan = QPointF()
            path_only = True
        if path_only:
            self.update(self._path_rect.toAlignedRect())
        else:
            self.update()

    def _draw_panel(self, painter, rect, title):
        painter.setPen(QPen(QColor("#c8d1dc"), 1.0))
        painter.setBrush(QColor(255, 255, 255, 242))
        painter.drawRoundedRect(rect, 8, 8)
        painter.setPen(QColor("#243548"))
        painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
        painter.drawText(
            rect.adjusted(9, 5, -8, -4),
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop |
            Qt.TextFlag.TextWordWrap,
            title)

    def _draw_slice(self, painter, visible):
        rect = self._slice_rect
        self._draw_panel(
            painter, rect,
            "caDNAno slice view")
        active = [self._rows[number] for number in visible]
        if not active:
            return
        active_coords = {(int(row["row"]), int(row["col"]))
                         for row in active}
        min_row = min(item[0] for item in active_coords) - 1
        max_row = max(item[0] for item in active_coords) + 1
        min_col = min(item[1] for item in active_coords) - 1
        max_col = max(item[1] for item in active_coords) + 1
        side = max(max_row - min_row + 1, max_col - min_col + 1)
        center_row = (min_row + max_row) / 2.0
        center_col = (min_col + max_col) / 2.0
        min_row = math.floor(center_row - (side - 1) / 2.0)
        min_col = math.floor(center_col - (side - 1) / 2.0)
        body = rect.adjusted(12, 29, -12, -10)
        # cadnano Slice uses a 15 px helix radius on a 30 px square lattice.
        fit = min(body.width(), body.height()) / max(30.0, side * 30.0)
        cell = 30.0 * fit * self._slice_zoom
        origin = QPointF(
            body.center().x() - ((side - 1) * cell) / 2.0 +
            self._slice_pan.x(),
            body.center().y() - ((side - 1) * cell) / 2.0 +
            self._slice_pan.y())
        coord_to_number = {
            (int(row["row"]), int(row["col"])): number
            for number, row in self._rows.items()}
        painter.save()
        painter.setClipRect(body)
        radius = 15.0 * fit * self._slice_zoom
        label_size = max(5, int(10.0 * fit * self._slice_zoom))
        for row_offset in range(side):
            for col_offset in range(side):
                coord = (min_row + row_offset, min_col + col_offset)
                center = QPointF(origin.x() + col_offset * cell,
                                 origin.y() + row_offset * cell)
                number = coord_to_number.get(coord)
                if coord in active_coords and number is not None:
                    face = next((item for item in self._capture_faces
                                 if number in item["helices"]), None)
                    if face is not None:
                        face_color = QColor(face["color"])
                        fill = QColor(face_color)
                        fill.setAlpha(210)
                        painter.setBrush(fill)
                        painter.setPen(QPen(
                            face_color.darker(125),
                            max(.7, .8 * fit * self._slice_zoom)))
                    else:
                        painter.setBrush(QColor("#ffeab7"))
                        painter.setPen(QPen(QColor("#ea8451"),
                            max(.5, .5 * fit * self._slice_zoom)))
                else:
                    painter.setBrush(QColor("#eeeeee"))
                    painter.setPen(QPen(QColor("#666666"),
                                        max(.5, .5 * fit)))
                painter.drawEllipse(center, radius, radius)
                if coord in active_coords and number is not None:
                    painter.setFont(
                        _preview_text_font(QFont.Weight.DemiBold))
                    painter.setPen(QColor("#222222"))
                    painter.drawText(
                        QRectF(center.x() - radius, center.y() - radius,
                               radius * 2, radius * 2),
                        Qt.AlignmentFlag.AlignCenter, str(number))
        # Label each capture face beside the corresponding Seed helix row.
        # Both labels share one x origin and remain horizontal/parallel as the
        # Slice view is zoomed or panned.
        painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
        active_centers = {
            number: QPointF(
                origin.x() + (int(row["col"]) - min_col) * cell,
                origin.y() + (int(row["row"]) - min_row) * cell)
            for number, row in self._rows.items()
            if number in visible}
        face_seed_numbers = {
            number for face in self._capture_faces
            for number in face["helices"]}
        rightmost_x = max(
            (center.x() for number, center in active_centers.items()
             if number in face_seed_numbers),
            default=body.center().x())
        face_label_x = rightmost_x + radius + max(8.0, cell * .28)
        for face_index, face in enumerate(self._capture_faces):
            centers = [active_centers[number]
                       for number in face["helices"]
                       if number in active_centers]
            if not centers:
                continue
            y = sum(center.y() for center in centers) / len(centers)
            painter.setBrush(QColor(face["color"]))
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawEllipse(QPointF(face_label_x, y), 4, 4)
            painter.setPen(QColor("#526375"))
            painter.drawText(QPointF(face_label_x + 7, y + 4),
                             "Face %d" % (face_index + 1))
        painter.restore()

    def _draw_capture_faces(self, painter):
        if self._channel != "complete" or not self._capture_seed:
            self._capture_rect = QRectF()
            return
        faces = self._capture_faces[:2]
        section = self._capture_rect
        self._draw_panel(
            painter, section,
            "Seed capture side view (line endpoints = staple extension sites)")
        gap = 8.0
        body_width = section.width() - 24.0
        panel_width = (body_width - gap) / 2.0
        panel_height = section.height() - 43.0
        left = section.left() + 12.0
        for face_index, face in enumerate(faces):
            numbers = face["helices"]
            rect = QRectF(
                left + face_index * (panel_width + gap),
                section.top() + 30.0, panel_width, panel_height)
            self._draw_panel(painter, rect, "Face %d" % (face_index + 1))
            occupied_by_helix = {}
            for number in numbers:
                row = self._rows.get(number, {})
                occupied = [
                    index for field in ("scaf", "stap")
                    for index, record in enumerate(row.get(field, []))
                    if self._occupied(record)]
                if occupied:
                    occupied_by_helix[number] = (min(occupied), max(occupied))
            if not occupied_by_helix:
                continue
            low = min(value[0] for value in occupied_by_helix.values())
            high = max(value[1] for value in occupied_by_helix.values())
            body = rect.adjusted(15, 47, -12, -28)
            line_body = body.adjusted(0, 18, 0, 0)
            column_fit = body.width() / max(12.0, len(numbers) * 22.0)
            column_step = 22.0 * column_fit * self._face_zoom
            base_fit = line_body.height() / max(
                20.0, (high - low + 1) * 12.0)
            base_scale = 12.0 * base_fit * self._face_zoom
            x_center = body.center().x() + self._face_pan.x()
            x_start = x_center - (len(numbers) - 1) * column_step / 2.0
            y_center = (line_body.center().y() + self._face_pan.y())

            def mapped_y(position):
                return y_center + (float(position) -
                                   (low + high) / 2.0) * base_scale

            painter.save()
            painter.setClipRect(body)
            for column_index, number in enumerate(numbers):
                if number not in occupied_by_helix:
                    continue
                x = x_start + column_index * column_step
                helix_low, helix_high = occupied_by_helix[number]
                y_low, y_high = mapped_y(helix_low), mapped_y(helix_high)
                painter.setFont(
                    _preview_text_font(QFont.Weight.DemiBold))
                painter.setPen(QColor("#35485b"))
                painter.drawText(
                    QRectF(x - column_step * .48, body.top(),
                           column_step * .96, 18.0),
                    Qt.AlignmentFlag.AlignHCenter |
                    Qt.AlignmentFlag.AlignTop, str(number))
                painter.setPen(QPen(QColor(face["color"]), 1.4))
                painter.drawLine(QPointF(x, y_low), QPointF(x, y_high))
                for position in sorted(self._capture_seed.get(number, ())):
                    if not (helix_low <= position <= helix_high):
                        continue
                    y = mapped_y(position)
                    protrusion = max(5.0, min(13.0, column_step * .36))
                    tip_x = x + self._capture_protrusion_side(
                        number) * protrusion
                    color_value = self._capture_seed_colors.get(
                        (number, position),
                        self._staple_colors.get(
                            (number, position), 0x000000))
                    site_color = QColor(
                        (color_value >> 16) & 255,
                        (color_value >> 8) & 255,
                        color_value & 255)
                    painter.setPen(QPen(
                        site_color.darker(135), 1.2))
                    painter.drawLine(QPointF(x, y), QPointF(tip_x, y))
                    painter.setBrush(site_color)
                    painter.drawEllipse(QPointF(tip_x, y), 3.2, 3.2)
                    painter.setFont(
                        _preview_text_font(QFont.Weight.DemiBold))
                    painter.setPen(site_color.darker(155))
                    label = self._capture_seed_ends.get(
                        number, {}).get(position, "")
                    text_x = tip_x + 4.0 if tip_x >= x else tip_x - 20.0
                    painter.drawText(QPointF(text_x, y - 4.0), label)
            painter.restore()

    def _draw_path(self, painter, visible, render_area=None,
                   draw_panel=True):
        """Rasterize Path content; interaction never calls this repeatedly."""
        area = render_area or self._path_rect
        path_title = (
            "caDNAno path view · scaffold routing only"
            if self._channel == "sequence_scaffold" else
            "caDNAno path view · SST sublattice input only"
            if self._channel == "sst_input" else
            "caDNAno path view")
        if draw_panel:
            self._draw_panel(
                painter, area, path_title)
            body = area.adjusted(54, 34, -14, -25)
        else:
            body = area.adjusted(42, 28, -22, -28)
        occupied_indices = []
        for number in visible:
            row = self._rows[number]
            fields = (("scaf",) if self._channel in
                      ("sequence_scaffold", "sst_input")
                      else ("scaf", "stap"))
            for field in fields:
                occupied_indices.extend(
                    index for index, record in enumerate(
                        self._path_records(number, field))
                    if self._occupied(record))
        if not occupied_indices:
            return
        base_low, base_high = min(occupied_indices), max(occupied_indices)
        base_span = max(1, base_high - base_low)
        # Preserve cadnano's 20-px base, 40-px helix and 50-px padding
        # proportions.  The initial fit is a view transform, like zoomToFit.
        fit_x = body.width() / max(20.0, (base_span + 1) * 20.0)
        fit_y = body.height() / max(40.0, len(visible) * 90.0)
        fit = max(.02, min(fit_x, fit_y))
        x_scale = 20.0 * fit
        helix_height = 40.0 * fit
        row_step = 90.0 * fit
        strand_width = max(.65, 3.0 * fit)
        scaffold_width = max(.9, 5.0 * fit)
        grid_width = max(.35, .5 * fit)
        content_width = (base_span + 1) * x_scale
        content_height = max(
            helix_height, (len(visible) - 1) * row_step + helix_height)
        y_start = body.center().y() - content_height / 2.0
        x_start = body.center().x() - content_width / 2.0
        helix_top = {number: y_start + index * row_step
                     for index, number in enumerate(visible)}

        def mapped_x(index):
            return x_start + (index - base_low + .5) * x_scale

        def is_drawn_5_to_3(number, field):
            even = number % 2 == 0
            return ((field == "scaf" and even) or
                    (field == "stap" and not even))

        def strand_y(number, field):
            on_top = is_drawn_5_to_3(number, field)
            return helix_top[number] + (
                .25 * helix_height if on_top else .75 * helix_height)

        painter.save()
        painter.setClipRect(body)

        for number in visible:
            top_y = helix_top[number]
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor("#999999"), grid_width))
            painter.drawRect(QRectF(
                mapped_x(base_low) - x_scale / 2.0, top_y,
                (base_span + 1) * x_scale, helix_height))
            painter.setPen(QPen(QColor("#cccccc"), grid_width))
            painter.drawLine(
                QPointF(mapped_x(base_low) - x_scale / 2.0,
                        top_y + helix_height / 2.0),
                QPointF(mapped_x(base_high) + x_scale / 2.0,
                        top_y + helix_height / 2.0))
            grid_step = 1 if x_scale >= 5.0 else 8
            first_grid = ((base_low + grid_step - 1) // grid_step) * grid_step
            for base in range(first_grid, base_high + 1, grid_step):
                x = mapped_x(base) - x_scale / 2.0
                color = "#999999" if base % 8 == 0 else "#cccccc"
                painter.setPen(QPen(QColor(color), grid_width))
                painter.drawLine(QPointF(x, top_y),
                                 QPointF(x, top_y + helix_height))
            handle_radius = max(4.0, 15.0 * fit)
            handle_center = QPointF(
                body.left() - 25.0, top_y + helix_height / 2.0)
            first_scaffold = next((index for index, record in enumerate(
                self._path_records(number, "scaf"))
                if self._occupied(record)), 0)
            handle_color = self._path_scaffold_color(number, first_scaffold)
            painter.setBrush(handle_color)
            painter.setPen(QPen(handle_color,
                                max(.7, 2.0 * fit)))
            painter.drawEllipse(handle_center, handle_radius, handle_radius)
            painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
            painter.setPen(QColor("#ffffff"))
            painter.drawText(
                QRectF(handle_center.x() - handle_radius,
                       handle_center.y() - handle_radius,
                       handle_radius * 2, handle_radius * 2),
                Qt.AlignmentFlag.AlignCenter, str(number))
            row = self._rows[number]
            scaffold_records = self._path_records(number, "scaf")
            scaffold_special = self._strand_special_nodes(
                number, scaffold_records)
            node_gap = max(.18, min(x_scale * .20, 4.0 * fit))
            y = strand_y(number, "scaf")
            for index in range(len(scaffold_records) - 1):
                if not self._records_linked(
                        number, scaffold_records, index, index + 1):
                    continue
                pen = QPen(self._path_scaffold_color(number, index),
                           scaffold_width)
                pen.setCapStyle(Qt.PenCapStyle.FlatCap)
                painter.setPen(pen)
                x1, x2 = mapped_x(index), mapped_x(index + 1)
                if index in scaffold_special:
                    x1 += node_gap
                if index + 1 in scaffold_special:
                    x2 -= node_gap
                if x2 > x1:
                    painter.drawLine(QPointF(x1, y), QPointF(x2, y))
            show_staple = self._path_shows_staple(number)
            if show_staple:
                staple_records = self._path_records(number, "stap")
                staple_special = self._strand_special_nodes(
                    number, staple_records)
                y = strand_y(number, "stap")
                for index in range(len(staple_records) - 1):
                    if not self._records_linked(
                            number, staple_records, index, index + 1):
                        continue
                    color_value = self._staple_colors.get(
                        (number, index), 0x000000)
                    color = QColor(
                        (color_value >> 16) & 255,
                        (color_value >> 8) & 255,
                        color_value & 255)
                    pen = QPen(color, strand_width)
                    pen.setCapStyle(Qt.PenCapStyle.FlatCap)
                    painter.setPen(pen)
                    x1, x2 = mapped_x(index), mapped_x(index + 1)
                    if index in staple_special:
                        x1 += node_gap
                    if index + 1 in staple_special:
                        x2 -= node_gap
                    if x2 > x1:
                        painter.drawLine(QPointF(x1, y), QPointF(x2, y))

        # Indels are part of the final design geometry.  Previously this
        # preview rendered connectivity only, so its Z2 looked nominal even
        # though the exported JSON contained the correct loop/skip arrays.
        for number in visible:
            row = self._rows[number]
            loops = row.get("loop", [])
            skips = row.get("skip", [])
            limit = min(len(loops), len(skips))
            for index in range(max(0, base_low), min(base_high + 1, limit)):
                loop_value = int(loops[index])
                skip_value = int(skips[index])
                if not loop_value and not skip_value:
                    continue
                occupied_fields = []
                for field in ("scaf", "stap"):
                    if field == "stap" and not self._path_shows_staple(number):
                        continue
                    records = row.get(field, [])
                    if index >= len(records):
                        continue
                    if self._path_draws_record(number, index, field) and \
                            self._occupied(records[index]):
                        occupied_fields.append(field)
                if not occupied_fields:
                    continue
                x = mapped_x(index)
                if loop_value:
                    field = ("scaf" if "scaf" in occupied_fields else
                             occupied_fields[0])
                    y = strand_y(number, field)
                    direction = (-1.0 if is_drawn_5_to_3(number, field)
                                 else 1.0)
                    half = max(2.5, min(7.0, x_scale * .42))
                    height = max(5.0, min(16.0, helix_height * .42))
                    loop_path = QPainterPath(QPointF(x - half, y))
                    loop_path.quadTo(
                        QPointF(x, y + direction * height),
                        QPointF(x + half, y))
                    painter.setBrush(Qt.BrushStyle.NoBrush)
                    painter.setPen(QPen(
                        QColor("#27a9e1"), max(1.1, 2.3 * fit)))
                    painter.drawPath(loop_path)
                    painter.setFont(_preview_text_font(
                        QFont.Weight.DemiBold))
                    painter.setPen(QColor("#1769aa"))
                    painter.drawText(
                        QPointF(x + half + 2.0,
                                y + direction * max(5.0, height * .55)),
                        "+%d" % loop_value)
                if skip_value:
                    y = helix_top[number] + helix_height / 2.0
                    half_x = max(2.8, min(7.0, x_scale * .38))
                    half_y = max(4.0, min(11.0, helix_height * .32))
                    painter.setPen(QPen(
                        QColor("#e33b3b"), max(1.2, 2.5 * fit)))
                    painter.drawLine(QPointF(x - half_x, y - half_y),
                                     QPointF(x + half_x, y + half_y))
                    painter.drawLine(QPointF(x - half_x, y + half_y),
                                     QPointF(x + half_x, y - half_y))

        fields = ["scaf", "stap"]
        for field in fields:
            for number in visible:
                if field == "stap" and not self._path_shows_staple(number):
                    continue
                row = self._rows[number]
                for index, record in enumerate(row.get(field, [])):
                    if not self._path_draws_record(number, index, field):
                        continue
                    partner, partner_index = map(int, record[2:4])
                    if partner == number or partner not in helix_top:
                        continue
                    if partner < 0:
                        continue
                    if field == "scaf":
                        color = self._path_scaffold_color(number, index)
                        width = scaffold_width
                    else:
                        color_value = self._staple_colors.get(
                            (number, index), 0x000000)
                        color = QColor(
                            (color_value >> 16) & 255,
                            (color_value >> 8) & 255,
                            color_value & 255)
                        width = strand_width
                    start = QPointF(
                        mapped_x(index), strand_y(number, field))
                    end = QPointF(
                        mapped_x(partner_index), strand_y(partner, field))
                    five_is_top = is_drawn_5_to_3(number, field)
                    five_is_5_to_3 = five_is_top
                    three_is_top = is_drawn_5_to_3(partner, field)
                    three_is_5_to_3 = three_is_top
                    half_base = x_scale / 2.0
                    five_enter = QPointF(
                        start.x() + (-half_base if five_is_5_to_3
                                     else half_base), start.y())
                    five_exit = QPointF(
                        start.x(), helix_top[number] +
                        (0 if five_is_top else helix_height))
                    three_enter = QPointF(
                        end.x(), helix_top[partner] +
                        (0 if three_is_top else helix_height))
                    three_exit = QPointF(
                        end.x() + (half_base if three_is_5_to_3
                                   else -half_base), end.y())
                    vertical_distance = abs(
                        three_enter.y() - five_exit.y())
                    if five_is_5_to_3 == three_is_5_to_3:
                        control = QPointF(
                            five_exit.x() + .035 * vertical_distance,
                            (five_exit.y() + three_enter.y()) / 2.0)
                    else:
                        sign = (-1.0 if five_is_top and five_is_5_to_3
                                else 1.0)
                        control = QPointF(
                            five_exit.x() + sign * .035 * vertical_distance,
                            (five_exit.y() + three_enter.y()) / 2.0)
                    # Match cadnano XoverItem: horizontal entry, vertical
                    # exit, one quadratic bridge, then the reciprocal node.
                    path = QPainterPath(five_enter)
                    path.lineTo(start)
                    path.lineTo(five_exit)
                    path.quadTo(control, three_enter)
                    path.lineTo(end)
                    path.lineTo(three_exit)
                    pen = QPen(color, width)
                    pen.setCapStyle(Qt.PenCapStyle.FlatCap)
                    painter.setPen(pen)
                    painter.setBrush(Qt.BrushStyle.NoBrush)
                    painter.drawPath(path)

        # cadnano endpoint grammar: a square is 5′ and a triangle is 3′.
        for field in fields:
            for number in visible:
                if field == "stap" and not self._path_shows_staple(number):
                    continue
                row = self._rows[number]
                drawn_5_to_3 = is_drawn_5_to_3(number, field)
                for index, record in enumerate(row.get(field, [])):
                    if not self._path_draws_record(number, index, field):
                        continue
                    if not self._occupied(record):
                        continue
                    if field == "scaf":
                        color = self._path_scaffold_color(number, index)
                    else:
                        color_value = self._staple_colors.get(
                            (number, index), 0x000000)
                        color = QColor(
                            (color_value >> 16) & 255,
                            (color_value >> 8) & 255,
                            color_value & 255)
                    center = QPointF(
                        mapped_x(index), strand_y(number, field))
                    cap = max(3.2, min(7.5, x_scale * .34))
                    painter.setPen(QPen(color, max(.6, strand_width * .55)))
                    painter.setBrush(color)
                    for offset, label in ((0, "5′"), (2, "3′")):
                        if int(record[offset]) >= 0:
                            continue
                        is_left = ((offset == 0 and drawn_5_to_3) or
                                   (offset == 2 and not drawn_5_to_3))
                        if offset == 0:
                            painter.drawRect(QRectF(
                                center.x() - cap, center.y() - cap,
                                cap * 2.0, cap * 2.0))
                        else:
                            point_x = (center.x() - cap if is_left else
                                       center.x() + cap)
                            back_x = (center.x() + cap if is_left else
                                      center.x() - cap)
                            triangle = QPolygonF([
                                QPointF(point_x, center.y()),
                                QPointF(back_x, center.y() - cap),
                                QPointF(back_x, center.y() + cap)])
                            painter.drawPolygon(triangle)
                        if x_scale >= 2.0:
                            painter.setFont(_preview_text_font(
                                QFont.Weight.DemiBold))
                            text_x = (center.x() - cap - 27 if is_left else
                                      center.x() + cap + 4)
                            painter.drawText(
                                QPointF(text_x, center.y() - cap - 3), label)

        if self._channel in ("sequence_scaffold", "sst_input") and \
                x_scale >= 3.0:
            painter.setFont(
                _preview_text_font(QFont.Weight.DemiBold))
            for (number, index), base in self._sequence_base_labels.items():
                if number not in helix_top or not (
                        base_low <= index <= base_high):
                    continue
                y = strand_y(number, "scaf")
                offset = (-max(3.0, 5.0 * fit)
                          if is_drawn_5_to_3(number, "scaf") else
                          max(7.0, 10.0 * fit))
                cell = QRectF(
                    mapped_x(index) - x_scale / 2.0, y + offset - 7.0,
                    x_scale, 12.0)
                painter.setPen(self._path_scaffold_color(number, index))
                painter.drawText(cell, Qt.AlignmentFlag.AlignCenter, base)

        painter.restore()

    def _draw_path_pixmap(self, painter, visible):
        """Transform one fixed HD raster; never repaint Path while interacting."""
        area = self._path_rect
        path_title = (
            "caDNAno path view · scaffold routing only"
            if self._channel == "sequence_scaffold" else
            "caDNAno path view · SST sublattice input only"
            if self._channel == "sst_input" else
            "caDNAno path view")
        self._draw_panel(
            painter, area, path_title)
        body = area.adjusted(12, 30, -12, -12)
        if self._path_pixmap is None:
            pixmap = QPixmap(
                self.PATH_IMAGE_WIDTH, self.PATH_IMAGE_HEIGHT)
            pixmap.fill(QColor("#ffffff"))
            image_painter = QPainter(pixmap)
            image_painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            self._draw_path(
                image_painter, visible,
                QRectF(0.0, 0.0, float(self.PATH_IMAGE_WIDTH),
                       float(self.PATH_IMAGE_HEIGHT)),
                draw_panel=False)
            image_painter.end()
            self._path_pixmap = pixmap
        base_scale = min(
            body.width() / max(1.0, self._path_pixmap.width()),
            body.height() / max(1.0, self._path_pixmap.height()))
        draw_width = self._path_pixmap.width() * base_scale * self._path_zoom
        draw_height = self._path_pixmap.height() * base_scale * self._path_zoom
        target = QRectF(
            body.center().x() - draw_width / 2.0 + self._path_pan.x(),
            body.center().y() - draw_height / 2.0 + self._path_pan.y(),
            draw_width, draw_height)
        visible_target = target.intersected(body)
        if visible_target.isEmpty():
            return
        source = QRectF(
            (visible_target.left() - target.left()) /
            max(1e-9, target.width()) * self._path_pixmap.width(),
            (visible_target.top() - target.top()) /
            max(1e-9, target.height()) * self._path_pixmap.height(),
            visible_target.width() / max(1e-9, target.width()) *
            self._path_pixmap.width(),
            visible_target.height() / max(1e-9, target.height()) *
            self._path_pixmap.height())
        painter.save()
        painter.setClipRect(body)
        painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        painter.drawPixmap(visible_target, self._path_pixmap, source)
        painter.restore()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#f3f6f9"))
        if not self._payload:
            painter.setFont(_preview_text_font(QFont.Weight.DemiBold))
            painter.setPen(QColor("#728193"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "Generate or import a design to display it here.")
            return
        self._layout_panels()
        visible = self._visible_rows()
        self._draw_path_report(painter)
        self._draw_path_pixmap(painter, visible)
        if self._channel not in ("sequence_scaffold", "sst_input"):
            self._draw_slice(painter, visible)
            self._draw_capture_faces(painter)
            self._draw_panel_handles(painter)

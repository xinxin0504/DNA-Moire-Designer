"""Bulk TEM/Moiré analysis UI and original-resolution exporters."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import tempfile
import time
from datetime import datetime
from pathlib import Path

from PyQt6.QtCore import QPointF, QRect, QRectF, QSize, QSignalBlocker, Qt
from PyQt6.QtGui import (
    QColor, QFont, QImage, QPainter, QPainterPath, QPen, QPixmap, QPolygonF)
from PyQt6.QtSvg import QSvgGenerator
from PyQt6.QtWidgets import (
    QApplication, QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QFrame, QGraphicsPixmapItem, QGraphicsScene, QGraphicsView,
    QGridLayout, QGroupBox, QHBoxLayout, QLabel,
    QInputDialog, QMessageBox, QPushButton, QScrollArea, QSizePolicy, QSlider, QSplitter,
    QStackedWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from .i18n import localize_csv, localize_svg, translate
from .scale_metadata import raw_pixel_size_nm, scale_nm_from_filename
from moire_runtime import worker_command


class PanZoomImageView(QGraphicsView):
    """Aspect-safe image viewer: wheel zoom, left-drag pan, double-click reset."""

    def __init__(self, placeholder="", parent=None):
        super().__init__(parent)
        self._scene = QGraphicsScene(self)
        self._item = QGraphicsPixmapItem()
        self._item.setTransformationMode(
            Qt.TransformationMode.SmoothTransformation)
        self._scene.addItem(self._item)
        self.setScene(self._scene)
        self._has_image = False
        self._user_view = False
        self.setMinimumSize(420, 680)
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setBackgroundBrush(QColor("#0b1014"))
        self.setRenderHints(QPainter.RenderHint.Antialiasing |
                            QPainter.RenderHint.TextAntialiasing |
                            QPainter.RenderHint.SmoothPixmapTransform)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(
            QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self.setToolTip((placeholder + "\n" if placeholder else "") +
                        "滚轮缩放 · 左键拖动平移 · 双击复位")

    def set_image(self, image):
        pixmap = QPixmap.fromImage(image)
        self._item.setPixmap(pixmap)
        self._scene.setSceneRect(QRectF(pixmap.rect()))
        self._has_image = not pixmap.isNull()
        self._user_view = False
        self.reset_view()

    def reset_view(self):
        if not self._has_image:
            return
        self.resetTransform()
        self.fitInView(self._item, Qt.AspectRatioMode.KeepAspectRatio)
        self._user_view = False

    def wheelEvent(self, event):
        if not self._has_image:
            return super().wheelEvent(event)
        factor = 1.18 if event.angleDelta().y() > 0 else 1.0 / 1.18
        current = self.transform().m11()
        target = current * factor
        if 0.03 <= target <= 80.0:
            self.scale(factor, factor)
            self._user_view = True
        event.accept()

    def mousePressEvent(self, event):
        if self._has_image and event.button() == Qt.MouseButton.LeftButton:
            self._user_view = True
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event):
        self.reset_view()
        event.accept()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._has_image and not self._user_view:
            self.fitInView(self._item, Qt.AspectRatioMode.KeepAspectRatio)


class AnalysisBulkMixin:
    """Methods mixed into ``MoireDesignerWindow``.

    All raster exports use the exact TEM source width and height and SVG view
    boxes use those same dimensions.  Embedded FFTs additionally use isotropic
    reciprocal-space axes, so a rectangular TEM cannot stretch FFT angles or
    turn a Square first-order peak set into a rectangle.
    """

    ANALYSIS_BUTTON_HEIGHT = 34
    ANNOTATION_DOTS_PER_METER = 3780  # 96 DPI

    @classmethod
    def _normalize_annotation_dpi(cls, image):
        """Prevent source-image DPI metadata from resizing annotations."""
        image.setDotsPerMeterX(cls.ANNOTATION_DOTS_PER_METER)
        image.setDotsPerMeterY(cls.ANNOTATION_DOTS_PER_METER)
        return image

    @staticmethod
    def _annotation_font(pixel_size=14, weight=QFont.Weight.Normal):
        """Return a true pixel-sized annotation font independent of DPI."""
        font = QFont("Arial")
        font.setPixelSize(int(pixel_size))
        font.setWeight(weight)
        return font

    def _standardize_analysis_button_heights(self, root):
        """Match Design action geometry and state colors in Analysis."""
        for button in root.findChildren(QPushButton):
            button.setFixedHeight(self.ANALYSIS_BUTTON_HEIGHT)
            if button.objectName() not in {
                    "optionalButton", "acceptedButton", "parameterStepButton"}:
                button.setObjectName("primaryButton")
            button.style().unpolish(button)
            button.style().polish(button)

    @staticmethod
    def _analysis_sidebar_width(widget, minimum=390, maximum=470):
        """Return a bounded width that accommodates the sidebar controls."""
        layout = widget.layout()
        if layout is not None:
            layout.activate()
        hinted_width = widget.sizeHint().width()
        if hinted_width <= 0:
            hinted_width = minimum
        return max(int(minimum), min(int(maximum), int(hinted_width)))

    def _build_analysis_tab(self):
        """Host the standalone Moiré analysis workspace."""
        tab = QWidget()
        outer = QVBoxLayout(tab)
        # Match the Design pages: keep the first analysis action 8 px below
        # the workflow bar and use the same lower inset.
        outer.setContentsMargins(0, 8, 0, 8)
        outer.setSpacing(0)
        self.analysis_module_stack = QStackedWidget()
        self.analysis_module_stack.addWidget(self._build_crystal_analysis_panel())
        outer.addWidget(self.analysis_module_stack, 1)
        self._standardize_analysis_button_heights(tab)
        return tab

    def _analysis_module_changed(self, index):
        self._open_analysis_module(index)

    def _open_analysis_module(self, index):
        index = max(0, min(int(index), self.analysis_module_stack.count() - 1))
        self.analysis_module_stack.setCurrentIndex(index)
        self._go_to_step(4 + index)

    def _build_crystal_analysis_panel(self):
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter, 1)

        # Left controls deliberately have no scroll area. They fill the full
        # column while result figures and detailed metrics scroll on the right.
        left_panel = QWidget()
        left = QVBoxLayout(left_panel)
        left.setContentsMargins(3, 2, 9, 4)
        left.setSpacing(4)
        left_panel.setMinimumWidth(390)
        left_panel.setMaximumWidth(470)
        splitter.addWidget(left_panel)

        # Compatibility-only summary target retained for older signal paths;
        # steps 1–3 are collapsible and step 4 remains visible at all times.
        self.analysis_intro = QLabel()
        self.analysis_intro.hide()

        upload_box = QGroupBox()
        upload_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        upload_layout = QVBoxLayout(upload_box)
        upload_layout.setSpacing(6)
        self.bulk_analysis_checkbox = QCheckBox("Bulk analysis（多图批量）")
        self.bulk_analysis_checkbox.setToolTip(
            "开启后可选择多个TEM样品；运行前选择文件夹，随后逐图自动导出并生成CSV汇总。")
        upload_layout.addWidget(self.bulk_analysis_checkbox)
        self.bulk_scale_widget = QWidget()
        self.bulk_scale_widget.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        bulk_scale_layout = QVBoxLayout(self.bulk_scale_widget)
        bulk_scale_layout.setContentsMargins(0, 2, 0, 2)
        bulk_scale_layout.setSpacing(4)
        bulk_scale_layout.addWidget(QLabel("Bulk scale 来源"))
        self.bulk_scale_mode = QComboBox()
        self.bulk_scale_mode.addItem("所有图片的 scale bar 数值相同", "same")
        self.bulk_scale_mode.addItem("各图数值不同（文件名开头标注）", "filename")
        self.bulk_scale_mode.addItem("原始 Raw data（读取内置像素尺度）", "raw")
        bulk_scale_layout.addWidget(self.bulk_scale_mode)
        self.bulk_common_scale_widget = QWidget()
        common_row = QHBoxLayout(self.bulk_common_scale_widget)
        common_row.setContentsMargins(0, 0, 0, 0)
        common_row.addWidget(QLabel("统一 scale bar 标注值"))
        self.bulk_common_scale_nm = QDoubleSpinBox()
        self.bulk_common_scale_nm.setRange(0, 100000)
        self.bulk_common_scale_nm.setDecimals(2)
        self.bulk_common_scale_nm.setSuffix(" nm")
        self.bulk_common_scale_nm.setSpecialValueText("请输入")
        common_row.addWidget(self.bulk_common_scale_nm, 1)
        bulk_scale_layout.addWidget(self.bulk_common_scale_widget)
        self.bulk_scale_help = QLabel()
        self.bulk_scale_help.setObjectName("subtitle")
        self.bulk_scale_help.setWordWrap(True)
        self.bulk_scale_help.setMinimumHeight(44)
        self.bulk_scale_help.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        bulk_scale_layout.addWidget(self.bulk_scale_help)
        self.bulk_scale_widget.hide()
        upload_layout.addWidget(self.bulk_scale_widget)
        tem_row = QHBoxLayout()
        self.select_tem_button = QPushButton("上传TEM图")
        self.tem_path_label = QLabel("尚未选择")
        self.tem_path_label.setWordWrap(True)
        self.tem_path_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        tem_row.addWidget(self.select_tem_button)
        tem_row.addWidget(self.tem_path_label, 1)
        upload_layout.addLayout(tem_row)

        # A single image is scale-checked immediately after selection.  Keep
        # these values next to the upload action because they gate every later
        # measurement; the result summary no longer doubles as an input form.
        self.single_scale_controls = QWidget()
        scale_controls_layout = QVBoxLayout(self.single_scale_controls)
        scale_controls_layout.setContentsMargins(0, 0, 0, 0)
        scale_form = QFormLayout()
        self.scale_bar_pixels = QDoubleSpinBox()
        self.scale_bar_pixels.setRange(0, 100000)
        self.scale_bar_pixels.setDecimals(1)
        self.scale_bar_pixels.setSuffix(" px")
        self.scale_bar_nm = QDoubleSpinBox()
        self.scale_bar_nm.setRange(0, 100000)
        self.scale_bar_nm.setDecimals(2)
        self.scale_bar_nm.setSuffix(" nm")
        scale_form.addRow("Scale bar pixel length", self.scale_bar_pixels)
        scale_form.addRow("Scale bar label value", self.scale_bar_nm)
        scale_controls_layout.addLayout(scale_form)
        self.scale_detection_status = QLabel(
            "Select one TEM image to detect its scale bar.")
        self.scale_detection_status.setWordWrap(True)
        scale_controls_layout.addWidget(self.scale_detection_status)
        upload_layout.addWidget(self.single_scale_controls)

        mode_box = QGroupBox()
        mode_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        mode_layout = QVBoxLayout(mode_box)
        mode_layout.setSpacing(6)
        self.image_analysis_mode = QComboBox()
        self.image_analysis_mode.addItem("双层分析", "bilayer")
        self.image_analysis_mode.addItem("单层分析", "single")
        mode_layout.addWidget(self.image_analysis_mode)
        self.analysis_selector_widget = QWidget()
        selector_row = QHBoxLayout(self.analysis_selector_widget)
        selector_row.setContentsMargins(0, 0, 0, 0)
        selector_row.addWidget(QLabel("当前样品"))
        self.analysis_file_selector = QComboBox()
        self.analysis_file_selector.setEnabled(False)
        selector_row.addWidget(self.analysis_file_selector, 1)
        self.analysis_selector_widget.hide()
        mode_layout.addWidget(self.analysis_selector_widget)

        self.bulk_overlay_widget = QWidget()
        bulk_overlay = QGridLayout(self.bulk_overlay_widget)
        bulk_overlay.setContentsMargins(0, 2, 0, 0)
        bulk_view_label = QLabel("Bulk Selected-spot view")
        bulk_view_label.setWordWrap(True)
        bulk_overlay.addWidget(bulk_view_label, 0, 0, 1, 2)
        self.bulk_ifft_view_mode = QComboBox()
        self.bulk_ifft_view_mode.addItem("Pure inverse FFT", "pure")
        self.bulk_ifft_view_mode.addItem("Overlay", "overlay")
        bulk_overlay.addWidget(self.bulk_ifft_view_mode, 1, 0, 1, 2)
        bulk_strength_label = QLabel("统一 Blend strength")
        bulk_strength_label.setWordWrap(True)
        bulk_overlay.addWidget(bulk_strength_label, 2, 0)
        self.bulk_ifft_strength_value = QLabel("100%")
        bulk_overlay.addWidget(self.bulk_ifft_strength_value, 2, 1)
        self.bulk_ifft_strength = QSlider(Qt.Orientation.Horizontal)
        self.bulk_ifft_strength.setRange(0, 100)
        self.bulk_ifft_strength.setValue(100)
        bulk_overlay.addWidget(self.bulk_ifft_strength, 3, 0, 1, 2)
        self.bulk_overlay_widget.hide()
        self.bulk_overlay_widget.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        calibration_box = QGroupBox()
        self.analysis_calibration_box = calibration_box
        calibration_layout = QVBoxLayout(calibration_box)
        self.single_result_summary = QWidget()
        summary_layout = QVBoxLayout(self.single_result_summary)
        summary_layout.setContentsMargins(0, 0, 0, 0)
        result_grid = QGridLayout()
        result_grid.setHorizontalSpacing(12)
        result_grid.setVerticalSpacing(6)
        result_grid.addWidget(QLabel(""), 0, 0)
        tem_header = QLabel("TEM")
        fft_header = QLabel("FFT")
        tem_header.setObjectName("analysisColumnHeader")
        fft_header.setObjectName("analysisColumnHeader")
        result_grid.addWidget(tem_header, 0, 1)
        result_grid.addWidget(fft_header, 0, 2)
        result_grid.addWidget(QLabel("a"), 1, 0)
        self.tem_a_value = QLabel("—")
        self.fft_a_value = QLabel("—")
        self.tem_a_value.setWordWrap(True)
        self.fft_a_value.setWordWrap(True)
        result_grid.addWidget(self.tem_a_value, 1, 1)
        result_grid.addWidget(self.fft_a_value, 1, 2)

        self.analysis_twist_label = QLabel("Twist")
        self.tem_twist_value = QLabel("—")
        self.fft_twist_value = QLabel("—")
        self.tem_twist_value.setWordWrap(True)
        self.fft_twist_value.setWordWrap(True)
        result_grid.addWidget(self.analysis_twist_label, 2, 0)
        result_grid.addWidget(self.tem_twist_value, 2, 1)
        result_grid.addWidget(self.fft_twist_value, 2, 2)
        self.analysis_period_label = QLabel("Period")
        self.tem_period_value = QLabel("—")
        self.fft_period_value = QLabel("—")
        self.tem_period_value.setWordWrap(True)
        self.fft_period_value.setWordWrap(True)
        result_grid.addWidget(self.analysis_period_label, 3, 0)
        result_grid.addWidget(self.tem_period_value, 3, 1)
        result_grid.addWidget(self.fft_period_value, 3, 2)
        summary_layout.addLayout(result_grid)
        calibration_layout.addWidget(self.single_result_summary)

        # Compatibility attributes retained for existing signal wiring and
        # project serialization; detailed values are displayed in labels.
        self.measured_lattice = QDoubleSpinBox()
        self.measured_moire = QDoubleSpinBox()
        for hidden in (self.measured_lattice, self.measured_moire):
            hidden.setRange(0, 100000)
            hidden.setDecimals(4)
            hidden.hide()
        self.tem_twist_result = self.tem_twist_value
        self.fft_twist_result = self.fft_twist_value
        self.final_twist_result = QLabel()
        self.final_twist_result.hide()

        self.image_analysis_status = QLabel("等待TEM图像。")
        self.image_analysis_status.setWordWrap(True)
        self.image_analysis_status.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
        self.image_analysis_status.setMinimumHeight(42)
        calibration_layout.addWidget(self.image_analysis_status)

        self.crystal_section_group = QButtonGroup(self)
        self.crystal_section_group.setExclusive(False)
        self.crystal_step_buttons = []

        def crystal_step_button(text, content):
            button = QPushButton(text)
            button.setObjectName("parameterStepButton")
            button.setCheckable(True)
            button.setFixedHeight(self.ANALYSIS_BUTTON_HEIGHT)
            self.crystal_section_group.addButton(button)
            button.toggled.connect(content.setVisible)
            button.setChecked(False)
            content.hide()
            left.addWidget(button)
            left.addWidget(content)
            self.crystal_step_buttons.append(button)
            return button

        self.crystal_upload_step_button = crystal_step_button(
            "1. 上传 TEM 图像", upload_box)
        self.crystal_mode_step_button = crystal_step_button(
            "2. 选择单层或双层分析", mode_box)
        # Step 3 owns both the optional FFT references and the analysis
        # action. Keep its controls hidden until the numbered step is opened.
        self.crystal_run_section = QWidget()
        run_section_layout = QVBoxLayout(self.crystal_run_section)
        run_section_layout.setContentsMargins(8, 6, 8, 8)
        run_section_layout.setSpacing(5)
        run_section_layout.addWidget(self._build_theoretical_reference_widget())
        self.run_image_analysis_button = QPushButton("Run Analysis")
        self.run_image_analysis_button.setObjectName("primaryButton")
        self.run_image_analysis_button.setFixedHeight(
            self.ANALYSIS_BUTTON_HEIGHT)
        run_section_layout.addWidget(self.run_image_analysis_button)
        run_section_layout.addWidget(self.bulk_overlay_widget)
        self.crystal_run_step_button = crystal_step_button(
            "3. Automatically Identify and Analyze",
            self.crystal_run_section)
        self.analysis_calibration_box.setTitle(translate(
            "4. 结果摘要"))
        left.addWidget(calibration_box)
        # Single-image numerical results live in the right-hand report.  The
        # left status box is retained only for the always-visible Bulk run
        # status requested by the batch workflow.
        self.analysis_calibration_box.hide()
        left.addStretch(1)
        crystal_sidebar_width = self._analysis_sidebar_width(
            left_panel, minimum=390, maximum=430)
        left_panel.setMinimumWidth(crystal_sidebar_width)
        left_panel.setMaximumWidth(crystal_sidebar_width + 80)

        # The complete optimized result is one vertically scrollable document.
        right_scroll = QScrollArea()
        right_scroll.setWidgetResizable(True)
        right_scroll.setFrameShape(QFrame.Shape.NoFrame)
        right_content = QWidget()
        right = QVBoxLayout(right_content)
        right.setContentsMargins(4, 2, 5, 8)
        right.setSpacing(8)
        right_scroll.setWidget(right_content)
        splitter.addWidget(right_scroll)

        self.analysis_preview_title = QLabel("TEM / FFT analysis")
        self.analysis_preview_title.setObjectName("title")
        self.analysis_preview_title.hide()
        right.addWidget(self.analysis_preview_title)
        self.analysis_figures = QWidget()
        figures = self.analysis_figures
        # Full-width figures remain stacked vertically; the outer right-hand
        # scroll area performs the page scrolling.
        figures_layout = QVBoxLayout(figures)
        figures_layout.setContentsMargins(0, 0, 0, 0)
        figures_layout.setSpacing(2)

        def figure_panel(title, subtitle, placeholder):
            # Use the same titled outer frame as Particle Analysis so all
            # analysis previews share one visual hierarchy.
            panel = QGroupBox(title)
            panel.setObjectName("previewGroupBox")
            layout = QVBoxLayout(panel)
            layout.setContentsMargins(10, 10, 8, 8)
            layout.setSpacing(4)
            heading = QLabel(subtitle)
            heading.setAlignment(Qt.AlignmentFlag.AlignCenter)
            heading.setWordWrap(True)
            heading.setObjectName("analysisFigureTitle")
            image = PanZoomImageView(placeholder)
            image.setSizePolicy(QSizePolicy.Policy.Expanding,
                                QSizePolicy.Policy.Fixed)
            image.setMinimumHeight(680)
            layout.addWidget(heading)
            layout.addWidget(image, 1)
            figures_layout.addWidget(panel)
            return image, heading, panel

        (self.tem_analysis_image, self.tem_analysis_heading,
         self.tem_analysis_box) = figure_panel(
            "Original TEM",
            "Real-space image and FFT",
            "Original TEM appears here")
        (self.reconstructed_analysis_image,
         self.reconstructed_analysis_heading,
         self.reconstructed_analysis_box) = figure_panel(
            "Reconstructed lattice",
            "Measured lattice model",
            "Reconstructed lattice appears here")
        (self.ifft_analysis_image, self.ifft_analysis_heading,
         self.ifft_analysis_box) = figure_panel(
            "Inverse FFT",
            "Phase-preserving selected-spot reconstruction",
            "Inverse FFT appears here")
        right.addWidget(figures)

        self.analysis_single_controls = QFrame()
        controls = self.analysis_single_controls
        controls.setObjectName("analysisControls")
        controls_layout = QHBoxLayout(controls)
        controls_layout.setContentsMargins(9, 7, 9, 7)
        self.analysis_overlay_label = QLabel("Selected-spot view")
        controls_layout.addWidget(self.analysis_overlay_label)
        self.ifft_view_mode = QComboBox()
        self.ifft_view_mode.addItem("Pure inverse FFT", "pure")
        self.ifft_view_mode.addItem("Overlay", "overlay")
        controls_layout.addWidget(self.ifft_view_mode)
        controls_layout.addWidget(QLabel("Blend strength"))
        self.ifft_strength = QSlider(Qt.Orientation.Horizontal)
        self.ifft_strength.setRange(0, 100)
        self.ifft_strength.setValue(100)
        self.ifft_strength_value = QLabel("100%")
        controls_layout.addWidget(self.ifft_strength, 1)
        controls_layout.addWidget(self.ifft_strength_value)
        self.save_image_analysis_button = QPushButton("Export current SVG")
        self.save_image_analysis_button.setEnabled(False)
        controls_layout.addWidget(self.save_image_analysis_button)
        right.addWidget(controls)

        self.analysis_cards = QWidget()
        analysis_cards = self.analysis_cards
        cards_layout = QHBoxLayout(analysis_cards)
        cards_layout.setContentsMargins(0, 0, 0, 0)
        cards_layout.setSpacing(2)

        def analysis_card(title):
            card = QGroupBox(title)
            layout = QVBoxLayout(card)
            text = QLabel("等待分析结果。")
            text.setWordWrap(True)
            text.setTextFormat(Qt.TextFormat.RichText)
            text.setAlignment(Qt.AlignmentFlag.AlignTop)
            layout.addWidget(text)
            cards_layout.addWidget(card, 1)
            return card, text

        self.analysis_lattice_card, self.analysis_lattice_text = analysis_card(
            "Bilayer lattice and Moiré analysis")
        self.analysis_comparison_card, self.analysis_comparison_text = analysis_card(
            "Comparison")
        self.analysis_ifft_card, self.analysis_ifft_text = analysis_card(
            "Selected-spot inverse FFT")
        right.addWidget(analysis_cards)
        self.analysis_preview_note = QLabel(
            "图片和具体分析数值使用与确认版相同的布局；所有图保持原始宽高比。"
            "三张图均支持滚轮缩放、左键拖动平移、双击复位。")
        self.analysis_preview_note.setWordWrap(True)
        self.analysis_preview_note.hide()
        right.addWidget(self.analysis_preview_note)
        self.bulk_no_preview_label = QLabel(
            "Bulk analysis在左侧逐张运行并自动导出；右侧不加载批量图片。")
        self.bulk_no_preview_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.bulk_no_preview_label.setWordWrap(True)
        self.bulk_no_preview_label.hide()
        right.addWidget(self.bulk_no_preview_label, 1)
        right.addStretch(1)
        splitter.setSizes([crystal_sidebar_width, 1100])

        self.measurements_table = QTableWidget(0, 5)
        self.measurements_table.setHorizontalHeaderLabels([
            "来源", "Twist angle", "Moiré period", "预测角度", "角度误差"])
        self._tem_image_path = None
        self._tem_image_paths = []
        self._fft_image_path = None
        self._tem_analysis = None
        self._fft_analysis = None
        self._analysis_records = []
        self._analysis_output_paths = {}
        self._final_analysis_angle = None
        self._tem_analysis_angle = None
        self._fft_analysis_angle = None
        self._preflight_scale_result = None
        self._analysis_bulk_mode_active = False
        self._single_analysis_state = None
        self._bulk_analysis_state = None
        self._update_single_scale_gate()
        return tab

    def _build_theoretical_reference_widget(self):
        widget = QGroupBox("FFT Recognition References (Optional)")
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(8, 6, 8, 6)
        layout.setSpacing(4)
        layout.addWidget(QLabel("Expected symmetries (maximum 2)"))
        self.theoretical_symmetry_rows = []
        self.theoretical_symmetry_rows_layout = QVBoxLayout()
        self.theoretical_symmetry_rows_layout.setSpacing(3)
        layout.addLayout(self.theoretical_symmetry_rows_layout)
        self._add_theoretical_symmetry_row(primary=True)
        layout.addWidget(QLabel("Expected lattice constants a (multiple values)"))
        self.theoretical_a_rows = []
        self.theoretical_a_rows_layout = QVBoxLayout()
        self.theoretical_a_rows_layout.setSpacing(3)
        layout.addLayout(self.theoretical_a_rows_layout)
        self._add_theoretical_a_row(primary=True)
        return widget

    def _add_theoretical_symmetry_row(self, primary=False):
        if len(getattr(self, "theoretical_symmetry_rows", [])) >= 2:
            return
        row_widget = QWidget()
        row = QHBoxLayout(row_widget)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)
        symmetry = QComboBox()
        symmetry.addItem("Automatic symmetry", "")
        symmetry.addItem("Square", "Square")
        symmetry.addItem("Honeycomb", "Honeycomb")
        symmetry.addItem("Kagome", "Kagome")
        action = QPushButton("+" if primary else "−")
        action.setFixedWidth(34)
        row.addWidget(symmetry, 1)
        row.addWidget(action)
        entry = {"widget": row_widget, "value": symmetry,
                 "action": action, "primary": primary}
        self.theoretical_symmetry_rows.append(entry)
        self.theoretical_symmetry_rows_layout.addWidget(row_widget)
        if primary:
            action.clicked.connect(self._add_theoretical_symmetry_row)
        else:
            action.clicked.connect(
                lambda unused=False, item=entry:
                self._remove_theoretical_symmetry_row(item))
        self._update_theoretical_symmetry_buttons()

    def _remove_theoretical_symmetry_row(self, entry):
        if (entry not in self.theoretical_symmetry_rows or
                entry["primary"]):
            return
        self.theoretical_symmetry_rows.remove(entry)
        entry["widget"].deleteLater()
        self._update_theoretical_symmetry_buttons()

    def _update_theoretical_symmetry_buttons(self):
        full = len(getattr(self, "theoretical_symmetry_rows", [])) >= 2
        for entry in getattr(self, "theoretical_symmetry_rows", []):
            if entry["primary"]:
                entry["action"].setEnabled(not full)
                entry["action"].setToolTip(
                    "Maximum of two symmetries" if full else
                    "Add another expected symmetry")

    def _add_theoretical_a_row(self, primary=False):
        row_widget = QWidget()
        row = QHBoxLayout(row_widget)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)
        value = QDoubleSpinBox()
        value.setRange(0.0, 100000.0)
        value.setDecimals(3)
        value.setSuffix(" nm")
        value.setSpecialValueText("Automatic a")
        action = QPushButton("+" if primary else "−")
        action.setFixedWidth(34)
        row.addWidget(value, 1)
        row.addWidget(action)
        entry = {"widget": row_widget, "value": value,
                 "action": action, "primary": primary}
        self.theoretical_a_rows.append(entry)
        self.theoretical_a_rows_layout.addWidget(row_widget)
        if primary:
            action.clicked.connect(self._add_theoretical_a_row)
        else:
            action.clicked.connect(
                lambda unused=False, item=entry:
                self._remove_theoretical_a_row(item))

    def _remove_theoretical_a_row(self, entry):
        if entry not in self.theoretical_a_rows or entry["primary"]:
            return
        self.theoretical_a_rows.remove(entry)
        entry["widget"].deleteLater()

    def _theoretical_reference_values(self):
        return {
            "symmetries": [
                entry["value"].currentData()
                for entry in getattr(self, "theoretical_symmetry_rows", [])
                if entry["value"].currentData()],
            "a_nm": [
                entry["value"].value()
                for entry in getattr(self, "theoretical_a_rows", [])
                if entry["value"].value() > 0],
        }

    def _set_theoretical_reference_values(self, references):
        references = references or {}
        # Accept the earlier paired-row state format when reopening an in-memory
        # mode snapshot created before this UI was rebuilt.
        if isinstance(references, list):
            symmetries = [item.get("symmetry") for item in references
                          if item.get("symmetry")][:2]
            a_values = [item.get("a_nm") for item in references
                        if item.get("a_nm")]
        else:
            symmetries = list(references.get("symmetries") or [])[:2]
            a_values = list(references.get("a_nm") or [])
        while len(self.theoretical_symmetry_rows) > 1:
            self._remove_theoretical_symmetry_row(
                self.theoretical_symmetry_rows[-1])
        if len(symmetries) > 1:
            self._add_theoretical_symmetry_row()
        for index, entry in enumerate(self.theoretical_symmetry_rows):
            value = symmetries[index] if index < len(symmetries) else ""
            symmetry_index = entry["value"].findData(value)
            entry["value"].setCurrentIndex(max(0, symmetry_index))
        while len(self.theoretical_a_rows) > 1:
            self._remove_theoretical_a_row(self.theoretical_a_rows[-1])
        for unused in range(max(0, len(a_values) - 1)):
            self._add_theoretical_a_row()
        for index, entry in enumerate(self.theoretical_a_rows):
            entry["value"].setValue(float(
                a_values[index] if index < len(a_values) else 0.0))

    def _bulk_enabled(self):
        return bool(self.bulk_analysis_checkbox.isChecked())

    def _bulk_analysis_changed(self, enabled):
        enabled = bool(enabled)
        previous_enabled = bool(getattr(
            self, "_analysis_bulk_mode_active", False))
        mode_changed = enabled != previous_enabled
        if mode_changed:
            current_state = self._capture_analysis_mode_state()
            if previous_enabled:
                self._bulk_analysis_state = current_state
            else:
                self._single_analysis_state = current_state
            self._analysis_bulk_mode_active = enabled
            target_state = (self._bulk_analysis_state if enabled else
                            self._single_analysis_state)
        else:
            target_state = None
        # Bulk results are intentionally never browsed in the right-hand
        # document.  They are rendered directly into the selected folder.
        self.analysis_selector_widget.setVisible(False)
        self.bulk_overlay_widget.setVisible(bool(enabled))
        self.bulk_scale_widget.setVisible(bool(enabled))
        self._bulk_scale_mode_changed()
        self.single_scale_controls.setVisible(not enabled)
        self.single_result_summary.hide()
        self.analysis_calibration_box.setVisible(bool(enabled))
        self.analysis_calibration_box.setTitle(translate(
            "4. 批量运行状态" if enabled else
            "4. 结果摘要"))
        for widget in (self.analysis_figures, self.analysis_single_controls,
                       self.analysis_cards):
            widget.setVisible(not enabled)
        self.analysis_preview_title.hide()
        self.analysis_preview_note.hide()
        self.bulk_no_preview_label.setVisible(bool(enabled))
        self.select_tem_button.setText(
            "上传多个TEM图" if enabled else "上传TEM图")
        self.crystal_run_step_button.setText(
            "3. Run Batch Detection and Analysis" if enabled else
            "3. Automatically Identify and Analyze")
        self.run_image_analysis_button.setText(
            "Run Batch Analysis" if enabled else "Run Analysis")
        self.save_image_analysis_button.setText(
            "Export all results + CSV" if enabled else
            "Export current SVG")
        self._update_crystal_intro()
        self.analysis_preview_note.setText(
            ("所有样品使用真实FFT和相位保持的Selected-spot孔径，并自动生成"
             "原始分辨率PNG、保持比例的SVG、JSON和TEM/FFT分列的CSV统计。"
             if enabled else
             "Original TEM、Reconstructed和Selected-spot IFFT保持确认版布局。"
             "分析完成后不会自动写出SVG；点击图片下方Export current result"
             "才按当前显示比例导出。"))
        if mode_changed:
            self._restore_analysis_mode_state(target_state, bulk=enabled)
        self.image_analysis_status.setMinimumHeight(70 if enabled else 42)
        self._update_single_scale_gate()

    def _capture_analysis_mode_state(self):
        return {
            "paths": list(getattr(self, "_tem_image_paths", [])),
            "path": getattr(self, "_tem_image_path", None),
            "records": list(getattr(self, "_analysis_records", [])),
            "record_index": self.analysis_file_selector.currentIndex(),
            "preflight": getattr(self, "_preflight_scale_result", None),
            "scale_pixels": self.scale_bar_pixels.value(),
            "scale_nm": self.scale_bar_nm.value(),
            "scale_status": self.scale_detection_status.text(),
            "analysis_status": self.image_analysis_status.text(),
            "tem_analysis": getattr(self, "_tem_analysis", None),
            "final_angle": getattr(self, "_final_analysis_angle", None),
            "tem_angle": getattr(self, "_tem_analysis_angle", None),
            "fft_angle": getattr(self, "_fft_analysis_angle", None),
            "theoretical_references": self._theoretical_reference_values(),
        }

    def _restore_analysis_mode_state(self, state, bulk=False):
        state = state or {}
        self._tem_image_paths = list(state.get("paths") or [])
        self._tem_image_path = state.get("path")
        self._analysis_records = list(state.get("records") or [])
        self._preflight_scale_result = state.get("preflight")
        self._tem_analysis = state.get("tem_analysis")
        self._final_analysis_angle = state.get("final_angle")
        self._tem_analysis_angle = state.get("tem_angle")
        self._fft_analysis_angle = state.get("fft_angle")
        self._set_theoretical_reference_values(
            state.get("theoretical_references"))
        blockers = [QSignalBlocker(self.scale_bar_pixels),
                    QSignalBlocker(self.scale_bar_nm)]
        try:
            self.scale_bar_pixels.setValue(float(
                state.get("scale_pixels") or 0.0))
            self.scale_bar_nm.setValue(float(state.get("scale_nm") or 0.0))
        finally:
            del blockers
        self.scale_detection_status.setText(state.get("scale_status") or
            "Select one TEM image to detect its scale bar.")
        self.image_analysis_status.setText(state.get("analysis_status") or
                                           "等待TEM图像。")
        names = [Path(path).name for path in self._tem_image_paths]
        self.tem_path_label.setText(
            (("%d个样品：%s" % (len(names), "、".join(names[:4]) +
                                ("…" if len(names) > 4 else "")))
             if bulk else (names[0] if names else "尚未选择")))
        self.analysis_file_selector.blockSignals(True)
        self.analysis_file_selector.clear()
        for record in self._analysis_records:
            self.analysis_file_selector.addItem(Path(record["source"]).name)
        index = int(state.get("record_index", 0))
        if self._analysis_records:
            index = max(0, min(index, len(self._analysis_records) - 1))
            self.analysis_file_selector.setCurrentIndex(index)
        self.analysis_file_selector.blockSignals(False)
        self.analysis_file_selector.setEnabled(
            bool(self._analysis_records) and bulk)
        self.save_image_analysis_button.setEnabled(bool(self._analysis_records))
        if not bulk:
            if self._analysis_records:
                self._activate_analysis_record(index)
            elif self._tem_image_path:
                self._show_scale_preflight_preview()

    def _bulk_scale_mode_changed(self):
        if not hasattr(self, "bulk_scale_mode"):
            return
        mode = self.bulk_scale_mode.currentData()
        self.bulk_common_scale_widget.setVisible(mode == "same")
        messages = {
            "same": (
                "Enter the shared scale-bar value before selecting images. "
                "Bar pixel lengths are detected per image."),
            "filename": (
                "Prefix each filename with '<value> <unit>_', for example "
                "'5 nm_sample.tif'."),
            "raw": (
                "Use raw TIFF files with embedded physical pixel size. "
                "Image scale bars are ignored."),
        }
        self.bulk_scale_help.setText(messages.get(mode, ""))
        self._update_crystal_intro()

    def _bulk_worker_scale_inputs(self, source):
        """Resolve Bulk scale before FFT/theoretical-a candidate filtering."""
        mode = (self.bulk_scale_mode.currentData()
                if hasattr(self, "bulk_scale_mode") else None)
        if mode == "raw":
            return raw_pixel_size_nm(source), None
        if mode == "filename":
            return None, scale_nm_from_filename(source)
        if mode == "same":
            value_nm = float(self.bulk_common_scale_nm.value())
            return None, value_nm if value_nm > 0 else None
        return None, None

    def _update_crystal_intro(self):
        if not hasattr(self, "analysis_intro"):
            return
        if self._bulk_enabled():
            mode = self.bulk_scale_mode.currentData()
            scale_step = {
                "same": "输入所有图共用的 scale bar 标注值，再批量选图",
                "filename": "按‘<数值> <单位>_文件名’为每张图标注尺度",
                "raw": "选择含内置物理像素尺度的原始 TIFF",
            }.get(mode, "选择 Bulk scale 模式")
            self.analysis_intro.setText(
                "1. %s\n2. 选择单层或双层分析\n"
                "3. 开始 Bulk 分析，逐图识别 TEM / FFT\n"
                "4. 自动导出每图结果与 CSV 汇总" % scale_step)
        else:
            self.analysis_intro.setText(
                "1. 上传 TEM 图；可以‘<数值> <单位>_’开头命名（如 20 nm_xx）\n"
                "2. 选择单层或双层分析\n"
                "3. 点击分析，自动识别 scale bar、TEM 与 FFT\n"
                "4. 核对尺度与结果；无 OCR/文件名尺度时才手动输入")

    def _overlay_mode(self):
        return (self.bulk_ifft_view_mode.currentData() if self._bulk_enabled()
                else self.ifft_view_mode.currentData())

    def _overlay_strength(self):
        return (self.bulk_ifft_strength.value() if self._bulk_enabled()
                else self.ifft_strength.value())

    def _analysis_mode_changed(self):
        bilayer = self.image_analysis_mode.currentData() == "bilayer"
        self.measured_moire.setEnabled(bilayer)
        for widget in (self.analysis_twist_label, self.tem_twist_value,
                       self.fft_twist_value, self.analysis_period_label,
                       self.tem_period_value, self.fft_period_value):
            widget.setVisible(bilayer)
        self.analysis_comparison_card.setVisible(bilayer)
        self.analysis_ifft_card.setVisible(bilayer)
        self.analysis_lattice_card.setTitle(
            "Bilayer lattice and Moiré analysis" if bilayer else
            "Lattice analysis")
        self.reconstructed_analysis_box.setTitle(
            "Reconstructed bilayer" if bilayer else
            "Reconstructed lattice")
        self.reconstructed_analysis_heading.setText(
            "Measured a and relative twist" if bilayer else
            "Measured FFT lattice constant")
        if self._analysis_records:
            self._analysis_records = []
            self.analysis_file_selector.clear()
            self.analysis_file_selector.setEnabled(False)
            self.save_image_analysis_button.setEnabled(False)
            self.image_analysis_status.setText("分析类型已改变，请重新运行分析。")
        self._analysis_values_changed()

    def _analysis_overlay_changed(self):
        self.ifft_strength_value.setText("%d%%" % self.ifft_strength.value())
        self.bulk_ifft_strength_value.setText(
            "%d%%" % self.bulk_ifft_strength.value())
        if not self._bulk_enabled():
            self._refresh_analysis_previews()

    def select_tem_image(self):
        if (self._bulk_enabled() and
                self.bulk_scale_mode.currentData() == "same" and
                self.bulk_common_scale_nm.value() <= 0):
            QMessageBox.information(
                self, "请先输入统一尺度",
                "请先输入这批图共用的 scale bar 标注值（nm），再选择图片。")
            self.bulk_common_scale_nm.setFocus()
            return
        file_filter = ("Raw TIFF (*.tif *.tiff)" if
                       self._bulk_enabled() and
                       self.bulk_scale_mode.currentData() == "raw" else
                       "Image (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if self._bulk_enabled():
            filenames, unused = QFileDialog.getOpenFileNames(
                self, "选择多个TEM图", str(Path.home() / "Desktop"),
                file_filter)
        else:
            filename, unused = QFileDialog.getOpenFileName(
                self, "选择TEM图", str(Path.home() / "Desktop"),
                "Image (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
            filenames = [filename] if filename else []
        if not filenames:
            return
        if (self._bulk_enabled() and
                self.bulk_scale_mode.currentData() == "filename"):
            unnamed = [Path(name).name for name in filenames
                       if scale_nm_from_filename(name) is None]
            if unnamed:
                QMessageBox.warning(
                    self, "文件名缺少尺度",
                    "当前 Bulk 模式要求每个文件以‘<数值> <单位>_’开头。\n"
                    "例：5 nm_sample.tif、0.2 µm_test.tif\n\n"
                    "请重命名后重新选择：\n" +
                    "\n".join(unnamed[:12]))
                return
        self._tem_image_paths = [str(Path(name).resolve()) for name in filenames]
        self._tem_image_path = self._tem_image_paths[0]
        names = [Path(name).name for name in self._tem_image_paths]
        self.tem_path_label.setText(
            ("%d个样品：%s" % (len(names), "、".join(names[:4]) +
                              ("…" if len(names) > 4 else "")))
            if self._bulk_enabled() else names[0])
        self._analysis_records = []
        self.analysis_file_selector.clear()
        self.analysis_file_selector.setEnabled(False)
        self.save_image_analysis_button.setEnabled(False)
        if not self._bulk_enabled():
            # Show the source immediately, then add the detected scale overlay
            # as soon as the lightweight preflight returns.
            source_image = QImage(self._tem_image_path)
            if not source_image.isNull():
                self.tem_analysis_image.set_image(source_image)
            self._detect_single_image_scale()
        else:
            self._update_single_scale_gate()

    @staticmethod
    def _run_scale_detection(source):
        with tempfile.TemporaryDirectory(prefix="moire-scale-detection-") as folder:
            pgm = Path(folder) / "image.pgm"
            AnalysisBulkMixin._host_image_to_pgm(source, pgm)
            result = subprocess.run(
                worker_command("scale-detection", str(pgm),
                               "--original", str(source)),
                check=False, text=True, capture_output=True, timeout=45)
            if result.returncode:
                raise RuntimeError(result.stderr.strip() or
                                   "Scale detection failed.")
            try:
                return json.loads(result.stdout)
            except Exception as error:
                raise RuntimeError(
                    "Scale detection returned an invalid result.") from error

    def _detect_single_image_scale(self):
        source = self._tem_image_path
        if not source:
            self._update_single_scale_gate()
            return
        self.scale_detection_status.setText("Detecting the scale bar…")
        QApplication.processEvents()
        self._busy(True)
        try:
            detected = self._run_scale_detection(source)
        except Exception as error:
            detected = {"scale_bar": None, "scale_value_nm": None,
                        "error": str(error)}
        finally:
            self._busy(False)
        filename_nm = scale_nm_from_filename(source)
        image_nm = float(detected.get("scale_value_nm") or 0.0)
        resolved_nm = image_nm or float(filename_nm or 0.0)
        detected["resolved_scale_nm"] = resolved_nm
        detected["scale_source"] = (
            "image OCR" if image_nm else
            "filename" if filename_nm else None)
        self._preflight_scale_result = detected
        bar = detected.get("scale_bar") or {}
        blockers = [QSignalBlocker(self.scale_bar_pixels),
                    QSignalBlocker(self.scale_bar_nm)]
        try:
            self.scale_bar_pixels.setValue(float(
                bar.get("pixel_length") or 0.0))
            self.scale_bar_nm.setValue(resolved_nm)
        finally:
            del blockers
        if bar and resolved_nm > 0:
            self.scale_detection_status.setText(
                "Scale detected: %.1f px = %g nm (%s)." % (
                    float(bar.get("pixel_length") or 0.0), resolved_nm,
                    detected.get("scale_source")))
        elif bar:
            self.scale_detection_status.setText(
                "Scale bar detected, but its label value was not recognized. "
                "Enter the value in nm to continue.")
        else:
            self.scale_detection_status.setText(
                "Scale bar not detected. Enter its pixel length and label "
                "value to continue.")
        self._show_scale_preflight_preview()
        self._update_single_scale_gate()

    def _show_scale_preflight_preview(self):
        if self._bulk_enabled() or not self._tem_image_path:
            return
        image = QImage(self._tem_image_path).convertToFormat(
            QImage.Format.Format_RGB32)
        if image.isNull():
            return
        self._normalize_annotation_dpi(image)
        detected = self._preflight_scale_result or {}
        bar = detected.get("scale_bar") or {}
        if bar:
            painter = QPainter(image)
            scale = max(1.0, min(image.width(), image.height()) / 1100.0)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            color = QColor("#ffd447")
            painter.setPen(QPen(color, max(2, int(3 * scale))))
            rect = QRectF(float(bar.get("x0", 0)), float(bar.get("y0", 0)),
                          max(1.0, float(bar.get("x1", 0)) -
                              float(bar.get("x0", 0))),
                          max(2.0, float(bar.get("y1", 0)) -
                              float(bar.get("y0", 0)) + 1.0))
            painter.drawRect(rect)
            value_nm = self.scale_bar_nm.value()
            text = ("Scale bar: %g nm" % value_nm if value_nm > 0 else
                    "Scale bar detected · value not recognized")
            # Keep this upload-time status annotation compact.  It is UI
            # guidance rather than a publication annotation, so high-resolution
            # TEM pixels must not scale it into an oversized banner.
            font = QFont("Arial")
            font.setPixelSize(14)
            font.setWeight(QFont.Weight.Bold)
            painter.setFont(font)
            metrics = painter.fontMetrics()
            text_width = metrics.horizontalAdvance(text) + 12
            text_height = metrics.height() + 6
            gap = 6
            # Keep the original scale label unobstructed: place our detection
            # label beside the bar, preferring the right side and falling back
            # to the left when the bar is close to the image edge.
            if rect.right() + gap + text_width <= image.width() - 2:
                text_x = rect.right() + gap
            else:
                text_x = max(2.0, rect.left() - gap - text_width)
            text_y = max(2.0, min(rect.center().y() - text_height / 2.0,
                                  image.height() - text_height - 2.0))
            text_rect = QRectF(text_x, text_y, text_width, text_height)
            painter.fillRect(text_rect, QColor(8, 15, 20, 205))
            painter.setPen(color)
            painter.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, text)
            painter.end()
        self.tem_analysis_image.set_image(image)

    def _update_single_scale_gate(self):
        if not hasattr(self, "crystal_mode_step_button"):
            return
        if self._bulk_enabled():
            # Numbered sections are navigation and remain available at every
            # stage. Required inputs gate only the final Run action.
            self.crystal_mode_step_button.setEnabled(True)
            ready = bool(self._tem_image_paths)
            # The numbered section is navigation and must remain openable;
            # only its Run action is gated by required input.
            self.crystal_run_step_button.setEnabled(True)
            self.run_image_analysis_button.setEnabled(ready)
            self.scale_bar_nm.setStyleSheet("")
            self.scale_bar_pixels.setStyleSheet("")
            return
        has_image = bool(self._tem_image_path)
        pixels_valid = self.scale_bar_pixels.value() > 0
        value_valid = self.scale_bar_nm.value() > 0
        ready = has_image and pixels_valid and value_valid
        self.crystal_mode_step_button.setEnabled(True)
        self.crystal_run_step_button.setEnabled(True)
        self.run_image_analysis_button.setEnabled(ready)
        invalid_style = (
            "QDoubleSpinBox { color: #a40000; background: #ffe6e6; "
            "border: 2px solid #d63c3c; }")
        self.scale_bar_pixels.setStyleSheet(
            invalid_style if has_image and not pixels_valid else "")
        self.scale_bar_nm.setStyleSheet(
            invalid_style if has_image and not value_valid else "")

    @staticmethod
    def _run_image_worker(source, analysis_kind, theoretical_references=None,
                          pixel_size_nm=None, scale_value_nm=None):
        folder = Path(tempfile.mkdtemp(prefix="moire-image-analysis-"))
        pgm = folder / "image.pgm"
        # This method is supplied by the host window.
        AnalysisBulkMixin._host_image_to_pgm(source, pgm)
        command = worker_command(
            "image-analysis", "tem", str(pgm),
            "--original", str(source),
            "--analysis-kind", str(analysis_kind),
            "--output-dir", str(folder))
        theoretical_references = theoretical_references or {}
        for symmetry in theoretical_references.get("symmetries") or []:
            command.extend(["--theoretical-symmetry", str(symmetry)])
        for a_nm in theoretical_references.get("a_nm") or []:
            command.extend(["--theoretical-a-nm", str(float(a_nm))])
        if pixel_size_nm and float(pixel_size_nm) > 0:
            command.extend(["--pixel-size-nm", str(float(pixel_size_nm))])
        if scale_value_nm and float(scale_value_nm) > 0:
            command.extend(["--scale-value-nm", str(float(scale_value_nm))])
        # The worker can legitimately spend tens of seconds on a 4k TEM and
        # returns a multi-megabyte JSON document.  ``subprocess.run`` blocked
        # the Qt event loop for the whole interval, making a successful run
        # look like a permanent hang.  Stream its output to files (avoiding a
        # full pipe deadlock), poll while servicing Qt events, and retain the
        # same hard timeout/error semantics.
        stdout_path = folder / "analysis_result.json"
        stderr_path = folder / "analysis_worker.log"
        started = time.monotonic()
        with stdout_path.open("w", encoding="utf-8") as stdout_handle, \
                stderr_path.open("w", encoding="utf-8") as stderr_handle:
            process = subprocess.Popen(
                command, stdout=stdout_handle, stderr=stderr_handle,
                text=True)
            while process.poll() is None:
                QApplication.processEvents()
                if time.monotonic() - started > 240.0:
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    raise RuntimeError("图像分析超过240秒，已安全停止。")
                time.sleep(.05)
        stderr_text = stderr_path.read_text(
            encoding="utf-8", errors="replace").strip()
        if process.returncode:
            raise RuntimeError(stderr_text or "图像分析进程失败。")
        try:
            payload = json.loads(stdout_path.read_text(encoding="utf-8"))
        except Exception as error:
            raise RuntimeError("图像分析没有返回有效结果。") from error
        payload["working_directory"] = str(folder)
        return payload

    @staticmethod
    def _host_image_to_pgm(source, target):
        image = QImage(str(source))
        if image.isNull():
            raise ValueError("Qt无法读取图像：%s" % source)
        gray = image.convertToFormat(QImage.Format.Format_Grayscale8)
        pointer = gray.constBits()
        pointer.setsize(gray.sizeInBytes())
        raw = bytes(pointer)
        rows = b"".join(
            raw[y * gray.bytesPerLine():
                y * gray.bytesPerLine() + gray.width()]
            for y in range(gray.height()))
        Path(target).write_bytes(
            ("P5\n%d %d\n255\n" % (gray.width(), gray.height())).encode(
                "ascii") + rows)

    def run_image_analysis(self):
        if not self._tem_image_paths:
            QMessageBox.information(
                self, "缺少TEM图",
                "请先上传多个TEM图。" if self._bulk_enabled() else
                "请先上传TEM图。")
            return
        if self._bulk_enabled():
            self._run_bulk_analysis()
            return
        if (self.scale_bar_pixels.value() <= 0 or
                self.scale_bar_nm.value() <= 0):
            self._update_single_scale_gate()
            return
        analysis_kind = self.image_analysis_mode.currentData()
        theoretical_references = self._theoretical_reference_values()
        pixel_size_nm = (self.scale_bar_nm.value() /
                         self.scale_bar_pixels.value())
        records, failures = [], []
        sources = (self._tem_image_paths if self._bulk_enabled() else
                   self._tem_image_paths[:1])
        self._busy(True)
        try:
            for index, source in enumerate(sources, 1):
                self.image_analysis_status.setText(
                    ("正在分析 %d/%d：%s" % (
                        index, len(sources), Path(source).name)
                     if self._bulk_enabled() else
                     "正在分析：%s" % Path(source).name))
                QApplication.processEvents()
                try:
                    result = self._run_image_worker(
                        source, analysis_kind, theoretical_references,
                        pixel_size_nm)
                    records.append({"source": source, "result": result})
                except Exception as error:
                    failures.append("%s：%s" % (Path(source).name, error))
        finally:
            self._busy(False)
        if records:
            # The single-image scale was already reviewed in Step 1.  Preserve
            # a manual correction instead of prompting again after the full
            # FFT/Moiré analysis.
            records[0]["scale_bar_pixels"] = self.scale_bar_pixels.value()
            records[0]["scale_bar_nm"] = self.scale_bar_nm.value()
            records[0]["scale_source"] = (
                (self._preflight_scale_result or {}).get("scale_source") or
                "manual")
        self._analysis_records = records
        self.analysis_file_selector.blockSignals(True)
        self.analysis_file_selector.clear()
        for record in records:
            self.analysis_file_selector.addItem(Path(record["source"]).name)
        self.analysis_file_selector.blockSignals(False)
        self.analysis_file_selector.setEnabled(
            bool(records) and self._bulk_enabled())
        self.save_image_analysis_button.setEnabled(bool(records))
        if records:
            self.analysis_file_selector.setCurrentIndex(0)
            self._activate_analysis_record(0)
        message = ("完成%d个样品" % len(records)
                   if self._bulk_enabled() else
                   ("分析完成：%s" % Path(records[0]["source"]).name
                    if records else "分析失败"))
        if failures:
            message += "；%d个失败：\n%s" % (len(failures), "\n".join(failures))
        warnings = [record.get("scale_warning") for record in records
                    if record.get("scale_warning")]
        if warnings:
            message += "\n尺度比对提示：" + "\n".join(warnings)
        self.image_analysis_status.setText(message)
        if not records:
            QMessageBox.critical(
                self, "批量分析失败" if self._bulk_enabled() else "图像分析失败",
                message)

    def _request_scale_components(self, record, bulk=True):
        """Resolve scale using OCR, filename/common fallback, or raw metadata."""
        result = record["result"]
        sample = Path(record["source"]).name
        bar = result.get("scale_bar") or {}
        pixels = float(bar.get("pixel_length") or 0.0)
        ocr_nm = float(result.get("scale_value_nm") or 0.0)
        filename_nm = scale_nm_from_filename(record["source"])
        mode = (self.bulk_scale_mode.currentData()
                if bulk and self._bulk_enabled() else "auto")
        if mode == "raw":
            pixel_nm = raw_pixel_size_nm(record["source"])
            if not pixel_nm:
                pixel_nm, accepted = QInputDialog.getDouble(
                    self, "Raw data 像素尺度：%s" % sample,
                    "文件中未读到物理像素尺度。请输入每像素对应的 nm：",
                    0.1, 0.000001, 100000.0, 6)
                if not accepted:
                    return False, "用户取消了 Raw data 像素尺度输入。"
                source_name = "manual_pixel_size"
            else:
                source_name = "raw_metadata"
            record["pixel_size_nm"] = float(pixel_nm)
            record["scale_source"] = source_name
            return True, ""

        fallback_nm = None
        fallback_source = None
        if mode == "same":
            fallback_nm = float(self.bulk_common_scale_nm.value())
            fallback_source = "bulk_common"
        elif filename_nm:
            fallback_nm = float(filename_nm)
            fallback_source = "filename"
        value_nm = ocr_nm or fallback_nm or 0.0
        source_name = "image_ocr" if ocr_nm else fallback_source
        comparison_nm = fallback_nm or filename_nm
        comparison_label = ("统一输入" if mode == "same" else "文件名")
        if ocr_nm and comparison_nm:
            tolerance = max(0.05, abs(comparison_nm) * 0.02)
            if abs(ocr_nm - comparison_nm) > tolerance:
                record["scale_warning"] = (
                    "OCR=%g nm，%s=%g nm；已优先使用 OCR，请核对。" %
                    (ocr_nm, comparison_label, comparison_nm))
        if pixels <= 0:
            pixels, accepted = QInputDialog.getDouble(
                self, "%sScale bar像素长度：%s" % (
                    "Bulk · " if bulk else "", sample),
                "无法识别scale bar横线。请输入该图scale bar的像素长度：",
                100.0, 0.1, 100000.0, 1)
            if not accepted:
                return False, "用户取消了scale bar像素长度输入。"
        if value_nm <= 0:
            value_nm, accepted = QInputDialog.getDouble(
                self, "%sScale bar标注值：%s" % (
                    "Bulk · " if bulk else "", sample),
                "无法读取scale bar标注数值。请输入该图对应的实际nm值：",
                100.0, 0.01, 100000.0, 2)
            if not accepted:
                return False, "用户取消了scale bar nm数值输入。"
            source_name = "manual"
        record["scale_bar_pixels"] = pixels
        record["scale_bar_nm"] = value_nm
        record["scale_source"] = source_name or "manual"
        record["filename_scale_nm"] = filename_nm
        return True, ""

    def _run_bulk_analysis(self):
        root = (self._project_output_dir("analysis/moire_twist")
                if hasattr(self, "_project_output_dir") else None)
        if root is None:
            root = QFileDialog.getExistingDirectory(
                self, "选择Bulk analysis自动导出文件夹",
                str(Path.home() / "Desktop"))
        if not root:
            return
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output = Path(root) / ("TEM_bulk_analysis_" + stamp)
        counter = 1
        while output.exists():
            output = Path(root) / (
                "TEM_bulk_analysis_%s_%d" % (stamp, counter))
            counter += 1
        output.mkdir(parents=True)
        rows, reports, failures = [], [], []
        analysis_kind = self.image_analysis_mode.currentData()
        theoretical_references = self._theoretical_reference_values()
        for index, source in enumerate(self._tem_image_paths, 1):
            self.image_analysis_status.setText(
                "Bulk %d/%d：正在分析并导出 %s" % (
                    index, len(self._tem_image_paths), Path(source).name))
            QApplication.processEvents()
            self._busy(True)
            try:
                pixel_size_nm, scale_value_nm = (
                    self._bulk_worker_scale_inputs(source))
                result = self._run_image_worker(
                    source, analysis_kind, theoretical_references,
                    pixel_size_nm=pixel_size_nm,
                    scale_value_nm=scale_value_nm)
            except Exception as error:
                failures.append("%s：%s" % (Path(source).name, error))
                continue
            finally:
                self._busy(False)
            record = {"source": source, "result": result}
            accepted, error = self._request_scale_components(
                record, bulk=True)
            if not accepted:
                failures.append("%s：%s" % (Path(source).name, error))
                continue
            sample = self._safe_stem(source)
            folder = output / sample
            suffix = 1
            while folder.exists():
                folder = output / (sample + "_%d" % suffix)
                suffix += 1
            try:
                report = self._export_analysis_record(record, folder)
                reports.append(report)
                rows.append(self._summary_row(report))
            except Exception as error:
                failures.append("%s：导出失败：%s" % (
                    Path(source).name, error))
        if rows:
            with (output / "bulk_analysis_summary.csv").open(
                    "w", newline="", encoding="utf-8-sig") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            localize_csv(output / "bulk_analysis_summary.csv")
            (output / "bulk_analysis_summary.json").write_text(
                json.dumps(reports, ensure_ascii=False, indent=2),
                encoding="utf-8")
        if failures:
            (output / "failed_samples.txt").write_text(
                "\n".join(failures), encoding="utf-8")
        # Bulk never populates the interactive result document on the right.
        self._analysis_records = []
        self.analysis_file_selector.clear()
        self.save_image_analysis_button.setEnabled(False)
        self.image_analysis_status.setText(
            "Bulk完成：%d个成功，%d个失败。\n自动导出：%s" %
            (len(rows), len(failures), output))
        QMessageBox.information(
            self, "Bulk analysis完成",
            "已逐张分析并自动导出。\n成功：%d\n失败：%d\n\n%s" %
            (len(rows), len(failures), output))

    def _activate_analysis_record(self, index):
        if index < 0 or index >= len(self._analysis_records):
            return
        record = self._analysis_records[index]
        self._tem_image_path = record["source"]
        self._tem_analysis = record["result"]
        result = record["result"]
        bar = result.get("scale_bar") or {}
        blockers = [QSignalBlocker(widget) for widget in (
            self.scale_bar_pixels, self.scale_bar_nm,
            self.measured_lattice, self.measured_moire)]
        try:
            self.scale_bar_pixels.setValue(float(
                record.get("scale_bar_pixels", bar.get("pixel_length") or 0)))
            self.scale_bar_nm.setValue(float(
                record.get("scale_bar_nm", result.get("scale_value_nm") or 0)))
            metrics = self._record_metrics(record)
            self.measured_lattice.setValue(float(
                record.get("lattice_nm", metrics.get("mean_a_nm") or 0)))
            self.measured_moire.setValue(float(
                record.get("moire_nm", metrics.get("tem_period_nm") or 0)))
        finally:
            del blockers
        self._analysis_values_changed()

    def _current_record(self):
        index = self.analysis_file_selector.currentIndex()
        return (self._analysis_records[index]
                if 0 <= index < len(self._analysis_records) else None)

    @staticmethod
    def _fft_lattice_entries(result, pixel_nm):
        layers = (result.get("lattice_fft") or {}).get("layers") or []
        if result.get("analysis_kind", "bilayer") == "bilayer":
            layers = layers[:2]
        raw = []
        for layer in layers:
            spacing_px = layer.get("lattice_constant_px")
            if not spacing_px:
                continue
            raw.append({
                "symmetry": str(layer.get("symmetry") or "Square"),
                "a_nm": float(spacing_px) * pixel_nm if pixel_nm else None,
                "orientation_deg": layer.get("orientation_deg"),
                "layer_role": layer.get("layer_role")})
        if not raw:
            spacing_px = (result.get("fft_lattice_constant_px") or
                          result.get("lattice_constant_px"))
            if spacing_px:
                raw.append({"symmetry": str(result.get("symmetry") or "Square"),
                            "a_nm": float(spacing_px) * pixel_nm
                            if pixel_nm else None,
                            "orientation_deg": None})
        # Collapse numerically identical layers, but keep genuinely different
        # lattice constants or symmetries as separate vertically listed values.
        merged = []
        for item in raw:
            match = None
            for old in merged:
                same_symmetry = old["symmetry"] == item["symmetry"]
                if old["a_nm"] is None or item["a_nm"] is None:
                    same_a = old["a_nm"] is item["a_nm"]
                else:
                    same_a = abs(old["a_nm"] - item["a_nm"]) <= max(
                        0.02, 0.02 * old["a_nm"])
                same_role = old.get("layer_role") == item.get("layer_role")
                if same_symmetry and same_a and same_role:
                    match = old
                    break
            if match is None:
                merged.append(dict(item, count=1))
            else:
                if match["a_nm"] is not None and item["a_nm"] is not None:
                    match["a_nm"] = ((match["a_nm"] * match["count"] +
                                      item["a_nm"]) / (match["count"] + 1))
                match["count"] += 1
        return merged

    @staticmethod
    def _domain_inter_axis_angles(domain):
        measured = domain.get("inter_axis_angles_deg") or []
        if measured:
            return [float(value) for value in measured]
        axes = sorted(float(value) % 180.0 for value in
                      (domain.get("reciprocal_axis_angles_deg") or []))
        symmetry = str(domain.get("symmetry") or "Square")
        if symmetry == "Square" and len(axes) >= 2:
            gap = axes[1]-axes[0]
            return [gap, 180.0-gap]
        if len(axes) >= 3:
            return [axes[1]-axes[0], axes[2]-axes[1],
                    180.0-axes[2]+axes[0]]
        return []

    @staticmethod
    def _spatial_domain_fields(result, pixel_nm):
        spatial = result.get("orientation_domains") or {}
        domains = spatial.get("domains") or []
        return {
            "Spatial_domain_symmetries": ";".join(
                str(domain.get("symmetry") or "Square")
                for domain in domains),
            "Spatial_domain_a_nm": ";".join(
                "%.6f" % (float(domain["lattice_constant_px"])*pixel_nm)
                for domain in domains
                if domain.get("lattice_constant_px") is not None and
                pixel_nm),
            "Spatial_domain_orientations_deg": ";".join(
                "%.6f" % float(domain["orientation_deg"])
                for domain in domains
                if domain.get("orientation_deg") is not None),
            "Spatial_domain_area_fraction": ";".join(
                "%.6f" % float(domain.get("area_fraction", 0.0))
                for domain in domains),
            "Spatial_domain_inter_axis_angles_deg": ";".join(
                "/".join("%.6f" % float(value) for value in
                         AnalysisBulkMixin._domain_inter_axis_angles(domain))
                for domain in domains),
            "Sample_area_fraction": spatial.get("sample_cell_fraction"),
            "Excluded_background_fraction": spatial.get(
                "background_fraction"),
            "Background_excluded": spatial.get("background_excluded"),
        }

    def _record_metrics(self, record, current_ui=False):
        result = record["result"]
        bar = result.get("scale_bar") or {}
        if current_ui:
            bar_px = self.scale_bar_pixels.value()
            bar_nm = self.scale_bar_nm.value()
        else:
            bar_px = float(record.get(
                "scale_bar_pixels", bar.get("pixel_length") or 0))
            bar_nm = float(record.get(
                "scale_bar_nm", result.get("scale_value_nm") or 0))
        pixel_nm = (float(record.get("pixel_size_nm"))
                    if record.get("pixel_size_nm") else
                    bar_nm / bar_px if bar_px > 0 and bar_nm > 0 else None)
        fft_a_px = result.get("fft_lattice_constant_px") or result.get(
            "lattice_constant_px")
        tem_period_px = result.get("moire_period_px")
        bilayer = result.get("analysis_kind", "bilayer") == "bilayer"
        mixed_multilayer = False
        assets = result.get("fft_assets") or {}
        fft_twist_reliable = bool(
            bilayer and result.get(
                "fft_twist_reliable",
                assets.get("first_order_bilayer_valid", False)))
        fft_period_px = (result.get("fft_predicted_moire_period_px")
                         if fft_twist_reliable else None)
        fft_a_nm = fft_a_px * pixel_nm if fft_a_px and pixel_nm else None
        pore = result.get("pore_lattice") or {}
        pore_a_px = (pore.get("lattice_constant_px")
                     if pore.get("valid") else None)
        pore_a_nm = pore_a_px * pixel_nm if pore_a_px and pixel_nm else None
        lattice_entries = self._fft_lattice_entries(result, pixel_nm)
        entry_values = [item["a_nm"] for item in lattice_entries
                        if item.get("a_nm") is not None]
        if entry_values and not mixed_multilayer:
            fft_a_nm = sum(entry_values) / len(entry_values)
        tem_reliable = bool(fft_twist_reliable and result.get(
            "tem_period_reliable",
            (result.get("moire_real_space") or {}).get("valid", False)))
        tem_period_nm = (tem_period_px * pixel_nm
                         if tem_reliable and tem_period_px and pixel_nm else None)
        fft_period_nm = (fft_period_px * pixel_nm
                         if fft_period_px and pixel_nm else None)
        tem_twist = None
        if bilayer and fft_a_nm and tem_period_nm and tem_period_nm >= fft_a_nm/2:
            tem_twist = math.degrees(2 * math.asin(
                min(1.0, fft_a_nm / (2 * tem_period_nm))))
        fft_twist = (result.get("fft_twist_angle_deg")
                     if fft_twist_reliable else None)
        return {
            "source_width_px": (result.get("fft_assets") or {}).get(
                "source_width_px") or QImage(record["source"]).width(),
            "source_height_px": (result.get("fft_assets") or {}).get(
                "source_height_px") or QImage(record["source"]).height(),
            "scale_bar_px": bar_px or None, "scale_bar_nm": bar_nm or None,
            "pixel_size_nm": pixel_nm,
            "scale_source": record.get("scale_source"),
            "scale_warning": record.get("scale_warning"),
            "tem_a_nm": None, "fft_a_nm": fft_a_nm,
            "pore_a_nm": pore_a_nm,
            "pore_to_helix_ratio": (pore.get("ratio_to_helix")
                                     if pore.get("valid") else None),
            "pore_repeat_multiple": (pore.get("repeat_multiple")
                                      if pore.get("valid") else None),
            "mean_a_nm": fft_a_nm,
            "layer_a_nm": entry_values,
            "fft_lattices": lattice_entries,
            "symmetries": [item["symmetry"] for item in lattice_entries],
            "tem_period_reliable": tem_reliable,
            "fft_twist_reliable": fft_twist_reliable,
            "mixed_multilayer": mixed_multilayer,
            "detected_layer_count": result.get("detected_layer_count"),
            "primary_twist_symmetry": (result.get("lattice_fft") or {}).get(
                "primary_twist_symmetry"),
            "tem_unavailable_reason": result.get(
                "tem_period_reliability_reason") or "",
            "fft_unavailable_reason": result.get(
                "fft_twist_reliability_reason") or (
                    "FFT中的两组同阶twist峰太近或证据不足，无法可靠区分。"
                    if bilayer and not fft_twist_reliable else ""),
            "tem_twist_deg": tem_twist, "fft_twist_deg": fft_twist,
            "tem_period_nm": tem_period_nm,
            "fft_period_nm": fft_period_nm,
            "final_twist_deg": tem_twist if bilayer and tem_reliable else None,
            "reconstruction_twist_deg": (
                tem_twist if tem_reliable else
                fft_twist if fft_twist_reliable else None),
            "reconstruction_period_nm": (
                tem_period_nm if tem_reliable else
                fft_period_nm if fft_twist_reliable else None),
        }

    def _scale_calibration_changed(self):
        self._update_single_scale_gate()
        if not self._bulk_enabled() and self._preflight_scale_result:
            self._show_scale_preflight_preview()
        record = self._current_record()
        if not record:
            return
        record["scale_bar_pixels"] = self.scale_bar_pixels.value()
        record["scale_bar_nm"] = self.scale_bar_nm.value()
        metrics = self._record_metrics(record, current_ui=True)
        blockers = [QSignalBlocker(self.measured_lattice),
                    QSignalBlocker(self.measured_moire)]
        try:
            self.measured_lattice.setValue(metrics.get("mean_a_nm") or 0)
            self.measured_moire.setValue(metrics.get("tem_period_nm") or 0)
        finally:
            del blockers
        self._analysis_values_changed()

    def _render_metrics(self, record):
        return self._record_metrics(
            record, current_ui=(not self._bulk_enabled() and
                                self._current_record() is record))

    def _analysis_values_changed(self):
        record = self._current_record()
        bilayer = self.image_analysis_mode.currentData() == "bilayer"
        if not record:
            self.tem_a_value.setText("—")
            self.fft_a_value.setText("—")
            self.tem_twist_value.setText("—")
            self.fft_twist_value.setText("—")
            self.tem_period_value.setText("—")
            self.fft_period_value.setText("—")
            return
        metrics = self._render_metrics(record)
        show_twist = bool(bilayer and metrics.get("fft_twist_reliable"))
        for widget in (self.analysis_twist_label, self.tem_twist_value,
                       self.fft_twist_value, self.analysis_period_label,
                       self.tem_period_value, self.fft_period_value):
            widget.setVisible(show_twist)
        self.measured_lattice.blockSignals(True)
        self.measured_moire.blockSignals(True)
        self.measured_lattice.setValue(metrics.get("fft_a_nm") or 0)
        self.measured_moire.setValue(metrics.get("tem_period_nm") or 0)
        self.measured_lattice.blockSignals(False)
        self.measured_moire.blockSignals(False)
        entries = metrics.get("fft_lattices") or []
        fft_a_lines = []
        for index, item in enumerate(entries, 1):
            prefix = ((item.get("symmetry") or "Lattice") + " "
                      if len(entries) > 1 else "")
            fft_a_lines.append(prefix + self._fmt(item.get("a_nm"), " nm"))
        if (record["result"].get("analysis_kind") == "single" and
                metrics.get("pore_a_nm") is not None):
            fft_a_lines = [
                "Helix: " + self._fmt(metrics.get("fft_a_nm"), " nm"),
                "Pore: " + self._fmt(metrics.get("pore_a_nm"), " nm")]
        self.tem_a_value.setText("—")
        self.fft_a_value.setText("<br>".join(fft_a_lines) or "—")
        self.tem_twist_value.setText(
            self._fmt(metrics.get("tem_twist_deg"), "°"))
        self.fft_twist_value.setText(
            self._fmt(metrics.get("fft_twist_deg"), "°"))
        self.tem_period_value.setText(
            self._fmt(metrics.get("tem_period_nm"), " nm"))
        self.fft_period_value.setText(
            self._fmt(metrics.get("fft_period_nm"), " nm"))
        self._tem_analysis_angle = metrics.get("tem_twist_deg")
        self._fft_analysis_angle = metrics.get("fft_twist_deg")
        self._final_analysis_angle = metrics.get("final_twist_deg")
        self._update_analysis_cards(record, metrics)
        self._update_analysis_headings(record, metrics)
        if not self._bulk_enabled():
            self._refresh_analysis_previews()

    @staticmethod
    def _metrics_html(rows):
        body = "".join(
            "<tr><td style='color:#657582;padding:3px 12px 3px 0'>%s</td>"
            "<td style='text-align:right;font-weight:600;padding:3px 0'>%s</td></tr>" %
            (label, value) for label, value in rows)
        return "<table width='100%%' cellspacing='0'>%s</table>" % body

    def _update_analysis_headings(self, record, metrics):
        reliable = metrics.get("tem_period_reliable")
        bilayer = record["result"].get("analysis_kind") == "bilayer"
        if not bilayer:
            original_subtitle = "Spatial lattice + FFT"
        elif reliable:
            original_subtitle = "Real-space Moiré + FFT"
        elif metrics.get("mixed_multilayer"):
            original_subtitle = "Mixed multilayer FFT"
        elif metrics.get("fft_twist_reliable"):
            original_subtitle = "FFT-derived twist available"
        else:
            original_subtitle = "Lattice a only"
        self.tem_analysis_box.setTitle("Original TEM")
        self.tem_analysis_heading.setText(original_subtitle)
        if bilayer:
            if metrics.get("mixed_multilayer"):
                subtitle = "%d layers · pair twist %s" % (
                    int(metrics.get("detected_layer_count") or 3),
                    self._fmt(metrics.get("fft_twist_deg"), "°"))
                title = "Mixed multilayer FFT model"
            elif metrics.get("fft_twist_reliable"):
                source_name = "TEM" if reliable else "FFT"
                subtitle = "a %s · %s twist %s" % (
                    self._fmt(metrics.get("mean_a_nm"), " nm"), source_name,
                    self._fmt(metrics.get("reconstruction_twist_deg"), "°"))
                title = "Reconstructed bilayer"
            else:
                subtitle = "a %s · twist unavailable" % (
                    self._fmt(metrics.get("mean_a_nm"), " nm"))
                title = "Lattice constant only"
        else:
            subtitle = "FFT lattice a %s" % self._fmt(
                metrics.get("mean_a_nm"), " nm")
            twin = (record["result"].get("lattice_fft") or {}).get("twin") or {}
            if twin.get("valid"):
                subtitle += " · Square twin Δ = %s" % self._fmt(
                    twin.get("relative_orientation_deg"), "°")
            if metrics.get("pore_a_nm") is not None:
                subtitle += " · pore a = %s" % self._fmt(
                    metrics.get("pore_a_nm"), " nm")
            title = "Reconstructed lattice"
        self.reconstructed_analysis_box.setTitle(title)
        self.reconstructed_analysis_heading.setText(subtitle)
        count = len((record["result"].get("fft_assets") or {}).get(
            "selected_spots", []))
        self.ifft_analysis_box.setTitle("Inverse FFT")
        self.ifft_analysis_heading.setText(
            "%d selected FFT peaks" % count)

    def _update_analysis_cards(self, record, metrics):
        bilayer = record["result"].get("analysis_kind") == "bilayer"
        entries = metrics.get("fft_lattices") or []
        lattice_rows = []
        for index, item in enumerate(entries, 1):
            suffix = " %d" % index if len(entries) > 1 else ""
            lattice_rows.append(("Symmetry%s" % suffix,
                                 item.get("symmetry") or "—"))
            label = ("Helix lattice constant, a%s" % suffix
                     if not bilayer else
                     "FFT lattice constant, a%s" % suffix)
            lattice_rows.append((label,
                                 self._fmt(item.get("a_nm"), " nm")))
            if item.get("orientation_deg") is not None:
                lattice_rows.append((
                    "Layer orientation%s" % suffix,
                    self._fmt(item.get("orientation_deg"), "°")))
        if not bilayer and metrics.get("pore_a_nm") is not None:
            lattice_rows += [
                ("Pore lattice constant, a",
                 self._fmt(metrics.get("pore_a_nm"), " nm")),
                ("Pore/helix spacing ratio",
                 self._fmt(metrics.get("pore_to_helix_ratio"), "×")),
                ("Shared lattice orientation",
                 self._fmt(((record["result"].get("pore_lattice") or {}).get(
                     "orientation_deg")), "°"))]
        references = record["result"]
        expected_symmetries = references.get("theoretical_symmetries") or []
        expected_a = references.get("theoretical_a_nm") or []
        if expected_symmetries:
            lattice_rows.append(("Reference symmetry",
                                 " / ".join(expected_symmetries)))
        if expected_a:
            lattice_rows.append((
                "Theoretical a reference",
                " / ".join("%g nm" % float(value) for value in expected_a)))
            lattice_rows.append((
                "a-reference matching",
                "Applied" if references.get("theoretical_a_filter_applied")
                else "No sufficient match; automatic result retained"))
        twin = (record["result"].get("lattice_fft") or {}).get("twin") or {}
        if not bilayer:
            orientations = record["result"].get(
                "single_layer_orientations_deg") or []
            spatial = record["result"].get("orientation_domains") or {}
            lattice_rows.extend(
                ("Lattice orientation %d" % (index+1),
                 self._fmt(angle, "°"))
                for index, angle in enumerate(orientations))
            if spatial.get("valid"):
                lattice_rows += [
                    ("Spatial domains", str(spatial.get("domain_count", 0))),
                    ("Domain boundaries", str(spatial.get(
                        "boundary_count", 0)))]
                if spatial.get("background_excluded"):
                    lattice_rows.append((
                        "Excluded non-sample background",
                        "%.1f%% of image" % (100.0*float(
                            spatial.get("background_fraction", 0.0)))))
                pixel_nm = metrics.get("pixel_size_nm")
                for domain in spatial.get("domains") or []:
                    domain_id = int(domain.get("domain_id", 0))
                    symmetry = str(domain.get("symmetry") or "Square")
                    a_px = domain.get("lattice_constant_px")
                    a_value = (float(a_px)*pixel_nm
                               if a_px is not None and pixel_nm else None)
                    lattice_rows.append((
                        "Domain %d" % domain_id,
                        "%s · a %s · %s · area %.1f%%" % (
                            symmetry,
                            self._fmt(a_value, " nm") if a_value is not None
                            else self._fmt(a_px, " px"),
                            self._fmt(domain.get("orientation_deg"), "°"),
                            100.0*float(domain.get("area_fraction", 0.0)))))
                    inter_angles = self._domain_inter_axis_angles(domain)
                    if inter_angles:
                        lattice_rows.append((
                            "Domain %d actual reciprocal %s-axis angles" %
                            (domain_id,
                             "two" if symmetry == "Square" else "tri"),
                            " / ".join("%.1f°" % float(value)
                                       for value in inter_angles)))
            if twin.get("valid") and len(orientations) == 2:
                lattice_rows.append(("Twin relative orientation", self._fmt(
                    twin.get("relative_orientation_deg"), "°")))
            lattice_rows.append((
                "Lattice-constant interpretation",
                "FFT-measured value; conventional-TEM specimen shrinkage "
                "is not silently corrected"))
        if bilayer and metrics.get("mixed_multilayer"):
            lattice_rows += [
                ("Detected model", "%d-layer mixed symmetry" % int(
                    metrics.get("detected_layer_count") or len(entries))),
                ("TEM-derived twist / moiré period",
                 "不报告：无法唯一归属于某一层对")]
        elif bilayer and metrics.get("tem_period_reliable"):
            lattice_rows += [
                ("TEM-derived moiré period",
                 self._fmt(metrics.get("tem_period_nm"), " nm")),
                ("TEM-derived twist",
                 self._fmt(metrics.get("tem_twist_deg"), "°"))]
        elif bilayer:
            lattice_rows.append((
                "TEM-derived twist",
                "无法识别：可靠Moiré单元少于2个"))
        self.analysis_lattice_text.setText(self._metrics_html(lattice_rows))
        if not bilayer:
            return
        if not metrics.get("fft_twist_reliable"):
            self.analysis_comparison_text.setText(self._metrics_html([
                ("TEM-derived twist",
                 "无法识别：可靠Moiré单元少于2个"),
                ("FFT-derived twist",
                 "无法识别：两组同阶峰太近或证据不足"),
                ("Reported measurement", "仅显示晶格常数 a")]))
            self.analysis_ifft_text.setText(self._metrics_html([
                ("Twist classification", "未报告：一阶峰不可可靠分离")]))
            return
        comparison_rows = [
            (("%s–%s twist from FFT" % (
                metrics.get("primary_twist_symmetry") or "Layer",
                metrics.get("primary_twist_symmetry") or "Layer"))
             if metrics.get("mixed_multilayer") else
             "FFT-derived twist",
             self._fmt(metrics.get("fft_twist_deg"), "°")),
            ("FFT-derived moiré period",
             self._fmt(metrics.get("fft_period_nm"), " nm"))]
        if metrics.get("mixed_multilayer"):
            comparison_rows += [
                ("Additional layer", "reported separately; no twist pair"),
                ("TEM-derived global moiré period",
                 "not uniquely assignable")]
            self.analysis_comparison_text.setText(
                self._metrics_html(comparison_rows))
            return
        if metrics.get("tem_period_reliable"):
            tem_twist = metrics.get("tem_twist_deg")
            fft_twist = metrics.get("fft_twist_deg")
            tem_period = metrics.get("tem_period_nm")
            fft_period = metrics.get("fft_period_nm")
            comparison_rows += [
                ("Mean TEM- and FFT-derived twist",
                 self._fmt((tem_twist + fft_twist) / 2
                           if tem_twist is not None and fft_twist is not None
                           else None, "°")),
                ("FFT−TEM twist deviation",
                 self._deviation_text(fft_twist, tem_twist, "°")),
                ("FFT−TEM period deviation",
                 self._deviation_text(fft_period, tem_period, " nm")),
                ("Geometric reconstruction", "a + twist only"),
                ("TEM registration", "Not phase matched")]
        else:
            comparison_rows.append(("TEM-derived twist / moiré period",
                                    "Not reported for this field of view"))
        self.analysis_comparison_text.setText(
            self._metrics_html(comparison_rows))
        assets = record["result"].get("fft_assets") or {}
        self.analysis_ifft_text.setText(self._metrics_html([
            ("Selected clear spots", "%d fitted peak centers" %
             len(assets.get("selected_spots", []))),
            ("Detection", ("Square reciprocal-lattice integer indexing"
                           if assets.get("selection_method") ==
                           "clear_square_integer_indexed_reflections" else
                           "Visible maxima + conjugate pairs")),
            ("Accepted reflection pairs", ", ".join(
                "%s %d" % (name.replace("_", " "), count)
                for name, count in ((assets.get("reflection_selection") or
                                     {}).get("pair_count_by_class") or
                                    {}).items()) or "—"),
            ("Aperture", "Local fitted ellipses"),
            ("Conjugate pairing", "Enforced"),
            ("Complex phase", "Preserved · registered"),
            ("SVG export", "Pure inverse FFT" if self._overlay_mode() == "pure"
             else "Overlay %d%%" % self._overlay_strength())]))

    @staticmethod
    def _deviation_text(value, reference, suffix):
        if value is None or reference is None:
            return "—"
        delta = float(value) - float(reference)
        percent = 100.0 * delta / float(reference) if reference else 0.0
        return "%+.2f%s (%+.1f%%)" % (delta, suffix, percent)

    @staticmethod
    def _fit_rect(bounds, image):
        if image.isNull() or bounds.width() <= 0 or bounds.height() <= 0:
            return QRectF(bounds)
        factor = min(bounds.width() / image.width(),
                     bounds.height() / image.height())
        width, height = image.width() * factor, image.height() * factor
        return QRectF(bounds.x() + (bounds.width() - width) / 2,
                      bounds.y() + (bounds.height() - height) / 2,
                      width, height)

    @staticmethod
    def _fit_fft_rect(bounds, image, source_rect):
        """Fit an FFT without distorting reciprocal-space geometry.

        A rectangular TEM produces a rectangular FFT array, but one FFT pixel
        represents ``1 / image_width`` in x and ``1 / image_height`` in y.
        Treating those array pixels as equal-sized display pixels therefore
        stretches angles and a Square reciprocal lattice.  Fit by the
        *normalised frequency span* instead: the complete FFT is square and a
        fractional centre crop keeps the same physical x/y scale.  This is the
        native, undistorted reciprocal-space aspect used by the on-screen
        preview and by both PNG and SVG exports.
        """
        if (image.isNull() or bounds.width() <= 0 or bounds.height() <= 0 or
                source_rect.width() <= 0 or source_rect.height() <= 0):
            return QRectF(bounds)
        frequency_width = source_rect.width() / float(image.width())
        frequency_height = source_rect.height() / float(image.height())
        aspect = frequency_width / max(frequency_height, 1e-12)
        if bounds.width() / bounds.height() > aspect:
            height = bounds.height()
            width = height * aspect
        else:
            width = bounds.width()
            height = width / max(aspect, 1e-12)
        return QRectF(bounds.x() + (bounds.width() - width) / 2,
                      bounds.y() + (bounds.height() - height) / 2,
                      width, height)

    @staticmethod
    def _load_gray(path):
        image = QImage(str(path))
        return image.convertToFormat(QImage.Format.Format_RGB32)

    @staticmethod
    def _periodic_mean_deg(values, period=90.0):
        """Circular mean for lattice orientations with ``period`` symmetry."""
        clean = [float(value) for value in values if value is not None]
        if not clean:
            return None
        multiplier = 2.0 * math.pi / period
        cosine = sum(math.cos(value * multiplier) for value in clean)
        sine = sum(math.sin(value * multiplier) for value in clean)
        if abs(cosine) + abs(sine) < 1e-12:
            return clean[0] % period
        return (math.atan2(sine, cosine) / multiplier) % period

    def _reconstruction_layer_angles(self, record, metrics=None):
        """Absolute lattice orientations for a TEM-registered reconstruction.

        FFT orientations are mathematical angles (positive y points upward).
        The returned values retain that convention; drawing code converts them
        to image coordinates.  When TEM real-space period supplies a refined
        relative twist, preserve the FFT-measured absolute mean orientation and
        distribute that refined twist symmetrically around the mean.
        """
        result = record["result"]
        if metrics is None:
            metrics = self._render_metrics(record)
        lattice_layers = (result.get("lattice_fft") or {}).get("layers") or []
        measured = [float(layer["orientation_deg"]) % 90.0
                    for layer in lattice_layers
                    if layer.get("orientation_deg") is not None]
        if not measured:
            single = ((result.get("fft_assets") or {}).get(
                "single_lattice") or {})
            if single.get("orientation_deg") is not None:
                measured = [float(single["orientation_deg"]) % 90.0]
        if result.get("analysis_kind") == "single":
            twin = (result.get("lattice_fft") or {}).get("twin") or {}
            if twin.get("valid") and len(measured) >= 2:
                return tuple(measured[:2])
            return (measured[0] if measured else 0.0,)
        if (result.get("mixed_multilayer") or
                (result.get("lattice_fft") or {}).get("mixed_multilayer")):
            return tuple(measured[:2])
        twist = float(metrics.get("reconstruction_twist_deg") or 0.0)
        mean = self._periodic_mean_deg(measured[:2], 90.0)
        if mean is None:
            mean = 0.0
        return ((mean - twist / 2.0) % 90.0,
                (mean + twist / 2.0) % 90.0)

    @staticmethod
    def _reconstruction_moire_basis(first_spacing, first_angle,
                                    second_spacing, second_angle,
                                    symmetry):
        """Return image-space basis vectors of the reconstructed Moiré cells.

        The coloured layer points are drawn from two rotated direct-lattice
        bases.  Their registry maxima therefore follow the inverse difference
        of the corresponding reciprocal bases; they do not, in general, lie
        on a horizontal/vertical grid.
        """
        def inverse(matrix):
            a, b = matrix[0]
            c, d = matrix[1]
            determinant = a*d-b*c
            if abs(determinant) < 1e-10:
                return None
            return ((d/determinant, -b/determinant),
                    (-c/determinant, a/determinant))

        def direct_basis(spacing, angle):
            radians = math.radians(-float(angle))
            ca, sa = math.cos(radians), math.sin(radians)
            if str(symmetry) == "Square":
                local = ((float(spacing), 0.0),
                         (0.0, float(spacing)))
            else:
                local = ((float(spacing), 0.0),
                         (float(spacing)*.5,
                          float(spacing)*math.sqrt(3.0)/2.0))
            vectors = []
            for x, y in local:
                vectors.append((ca*x-sa*y, sa*x+ca*y))
            return ((vectors[0][0], vectors[1][0]),
                    (vectors[0][1], vectors[1][1]))

        first_inverse = inverse(direct_basis(
            first_spacing, first_angle))
        second_inverse = inverse(direct_basis(
            second_spacing, second_angle))
        if first_inverse is None or second_inverse is None:
            return None
        reciprocal_difference = (
            (second_inverse[0][0]-first_inverse[0][0],
             second_inverse[0][1]-first_inverse[0][1]),
            (second_inverse[1][0]-first_inverse[1][0],
             second_inverse[1][1]-first_inverse[1][1]))
        moire_matrix = inverse(reciprocal_difference)
        if moire_matrix is None:
            return None
        first = (moire_matrix[0][0], moire_matrix[1][0])
        second = (moire_matrix[0][1], moire_matrix[1][1])
        if not all(math.isfinite(value)
                   for vector in (first, second) for value in vector):
            return None
        return first, second

    def _reconstructed_image(self, record, painter=None):
        source = QImage(record["source"])
        width, height = source.width(), source.height()
        owns_painter = painter is None
        image = None
        if owns_painter:
            image = QImage(width, height, QImage.Format.Format_RGB32)
            image.fill(QColor("black"))
            painter = QPainter(image)
        else:
            painter.fillRect(QRectF(0, 0, width, height), QColor("black"))
        result = record["result"]
        metrics = self._render_metrics(record)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        cx, cy = width / 2, height / 2
        colors = (QColor(255, 189, 63, 184), QColor(43, 212, 240, 214),
                  QColor(112, 240, 140, 205), QColor(190, 120, 245, 205))
        angles = self._reconstruction_layer_angles(record, metrics)
        raw_layers = (result.get("lattice_fft") or {}).get("layers") or []
        default_spacing = float(result.get("fft_lattice_constant_px") or 12.0)
        components = []
        split_unclassified = False
        if (result.get("analysis_kind") == "bilayer" and
                not metrics.get("fft_twist_reliable")):
            angle = angles[0] if angles else 0.0
            components.append((default_spacing, angle, colors[0],
                               False, "lattice"))
        elif result.get("analysis_kind") == "bilayer":
            for index, angle in enumerate(angles):
                spacing = (float(raw_layers[index].get("lattice_constant_px"))
                           if index < len(raw_layers) and
                           raw_layers[index].get("lattice_constant_px") else
                           default_spacing)
                symmetry = (str(raw_layers[index].get("symmetry") or "Square")
                            if index < len(raw_layers) else "Square")
                components.append((spacing, angle, colors[index % len(colors)],
                                   index > 0, "layer", symmetry))
        else:
            primary_angle = angles[0] if angles else 0.0
            components.append((default_spacing, primary_angle, colors[0],
                               False, "helix"))
            twin = (result.get("lattice_fft") or {}).get("twin") or {}
            pore = result.get("pore_lattice") or {}
            heterogeneous = bool((result.get("lattice_fft") or {}).get(
                "heterogeneous_domains"))
            if heterogeneous and raw_layers:
                components = [(
                    float(layer.get("lattice_constant_px") or
                          default_spacing),
                    float(layer.get("orientation_deg") or 0.0),
                    colors[index % len(colors)], False,
                    "heterogeneous_domain",
                    str(layer.get("symmetry") or "Square"))
                    for index, layer in enumerate(raw_layers[:4])]
                split_unclassified = True
            elif twin.get("valid") and len(raw_layers) >= 2:
                components = []
                for index, layer in enumerate(raw_layers[:2]):
                    components.append((
                        float(layer.get("lattice_constant_px") or
                              default_spacing),
                        float(layer.get("orientation_deg") or 0.0),
                        colors[index], False, "twin_domain"))
                # A global FFT identifies both orientations but not a unique
                # real-space boundary; show the two measured domains side by
                # side instead of inventing a mixed full-field lattice.
                split_unclassified = True
            elif pore.get("valid") and pore.get("lattice_constant_px"):
                components.append((float(pore["lattice_constant_px"]),
                                   float(pore.get("orientation_deg",
                                                  primary_angle)),
                                   QColor(112, 240, 140, 225), True,
                                   "supercell"))
            else:
                # If two distinct a values are present but the frequency-only
                # evidence cannot establish supercell/twin/domain topology,
                # do not invent an overlay.  Reconstruct them side by side as
                # requested, each with its measured spacing and orientation.
                distinct = []
                for layer in raw_layers:
                    spacing = layer.get("lattice_constant_px")
                    if not spacing or any(abs(float(spacing)-old[0]) /
                                          max(old[0], 1e-9) < .05
                                          for old in distinct):
                        continue
                    distinct.append((float(spacing),
                                     float(layer.get("orientation_deg") or 0.0)))
                if len(distinct) >= 2:
                    components = [
                        (distinct[0][0], distinct[0][1], colors[0], False,
                         "unclassified_a1"),
                        (distinct[1][0], distinct[1][1], colors[1], False,
                         "unclassified_a2")]
                    split_unclassified = True

        def draw_component(spacing, angle, color, outline, clip=None,
                           symmetry="Square"):
            radius = max(1.0, spacing * (0.105 if outline else 0.13))
            extent = int(math.hypot(width, height) / max(spacing, 1)) + 3
            # FFT uses y-up mathematical angles; QImage uses y-down.
            radians = math.radians(-angle)
            ca, sa = math.cos(radians), math.sin(radians)
            painter.save()
            if clip is not None:
                painter.setClipRect(clip)
            if outline:
                painter.setPen(QPen(color, max(1.0, radius * 0.34)))
                painter.setBrush(Qt.BrushStyle.NoBrush)
            else:
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(color)
            for row in range(-extent, extent + 1):
                for column in range(-extent, extent + 1):
                    if symmetry == "Square":
                        x0, y0 = column * spacing, row * spacing
                    else:
                        x0 = (column + .5*(row & 1))*spacing
                        y0 = row*spacing*math.sqrt(3.0)/2.0
                    x = cx + ca * x0 - sa * y0
                    y = cy + sa * x0 + ca * y0
                    if -radius <= x <= width + radius and -radius <= y <= height + radius:
                        painter.drawEllipse(QPointF(x, y), radius, radius)
            painter.restore()

        if split_unclassified:
            count = max(1, len(components))
            clips = [QRectF(index*width/count, 0, width/count, height)
                     for index in range(count)]
            for component, clip in zip(components, clips):
                symmetry = component[5] if len(component) > 5 else "Square"
                draw_component(*component[:4], clip=clip,
                               symmetry=symmetry)
            painter.setPen(QPen(QColor(220, 228, 233, 180), 2.0,
                                Qt.PenStyle.DashLine))
            for index in range(1, count):
                painter.drawLine(QPointF(index*width/count, 0),
                                 QPointF(index*width/count, height))
        else:
            for component in components:
                symmetry = component[5] if len(component) > 5 else "Square"
                draw_component(*component[:4], symmetry=symmetry)
        if (result.get("analysis_kind") == "bilayer" and
                not metrics.get("mixed_multilayer")):
            period = float(result.get("moire_period_px") or
                           result.get("fft_predicted_moire_period_px") or 0.0)
            if period > 0:
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(QColor(112, 240, 140, 225))
                moire_basis = None
                if len(components) >= 2:
                    first_component, second_component = components[:2]
                    first_symmetry = (first_component[5]
                                      if len(first_component) > 5
                                      else "Square")
                    second_symmetry = (second_component[5]
                                       if len(second_component) > 5
                                       else "Square")
                    if first_symmetry == second_symmetry:
                        moire_basis = self._reconstruction_moire_basis(
                            first_component[0], first_component[1],
                            second_component[0], second_component[1],
                            first_symmetry)
                if moire_basis:
                    first_basis, second_basis = moire_basis
                    shortest = min(
                        math.hypot(*first_basis),
                        math.hypot(*second_basis))
                    count = int(math.hypot(width, height) /
                                max(shortest, 1.0)) + 3
                    for row in range(-count, count + 1):
                        for column in range(-count, count + 1):
                            x = (cx + column*first_basis[0] +
                                 row*second_basis[0])
                            y = (cy + column*first_basis[1] +
                                 row*second_basis[1])
                            if 0 <= x <= width and 0 <= y <= height:
                                painter.drawEllipse(
                                    QPointF(x, y),
                                    max(2.0, default_spacing*.07),
                                    max(2.0, default_spacing*.07))
                else:
                    count = int(math.hypot(width, height) / period) + 2
                    for row in range(-count, count + 1):
                        for column in range(-count, count + 1):
                            x, y = cx + column*period, cy + row*period
                            if 0 <= x <= width and 0 <= y <= height:
                                painter.drawEllipse(
                                    QPointF(x, y),
                                    max(2.0, default_spacing*.07),
                                    max(2.0, default_spacing*.07))
        if owns_painter:
            painter.end()
        return image

    def _ifft_image(self, record):
        original = QImage(record["source"]).convertToFormat(
            QImage.Format.Format_RGB32)
        path = (record["result"].get("fft_assets") or {}).get(
            "reconstruction_path")
        reconstruction = self._load_gray(path) if path else QImage()
        if reconstruction.isNull():
            return original
        if self._overlay_mode() == "pure":
            return reconstruction.scaled(
                original.size(), Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation)
        result = original.copy()
        painter = QPainter(result)
        painter.setCompositionMode(
            QPainter.CompositionMode.CompositionMode_Overlay)
        painter.setOpacity(self._overlay_strength() / 100.0)
        target = self._fit_rect(QRectF(result.rect()), reconstruction)
        painter.drawImage(target, reconstruction, QRectF(reconstruction.rect()))
        painter.end()
        return result

    @staticmethod
    def _uploaded_annotation_units(source):
        """Return geometry and text scales for annotations on Uploaded TEM.

        Geometry remains resolution-aware so boxes, markers and line weights
        keep their relative placement on the exported original-resolution
        image.  Font point sizes deliberately do *not* scale with raster pixel
        count: a cryo/high-resolution upload must use the same annotation type
        size as an ordinary TEM rather than treating extra pixels as a reason
        to enlarge its text.
        """
        unit = max(1.0, min(source.width() / 1350.0,
                            source.height() / 1040.0))
        text_unit = 1.0
        return unit, text_unit

    def _draw_fft_inset(self, painter, record, selected=False):
        result = record["result"]
        metrics = self._render_metrics(record)
        bilayer = result.get("analysis_kind") == "bilayer"
        assets = result.get("fft_assets") or {}
        fft = self._load_gray(assets.get("fft_path", ""))
        if fft.isNull():
            return
        source = QImage(record["source"])
        unit, text_unit = self._uploaded_annotation_units(source)
        margin = 22 * unit if not selected else 28 * unit
        inset_width = (448 if not selected else 410) * unit
        spatial_domains = ((result.get("orientation_domains") or {}).get(
            "domains") or [])
        domain_rows = int(math.ceil(len(spatial_domains)/2.0))
        if (not selected and result.get("analysis_kind") == "single" and
                spatial_domains):
            inset_height = max(405.0, 330.0+38.0*domain_rows)*unit
        else:
            inset_height = (405 if not selected else 326)*unit
        box = QRectF(source.width() - inset_width - margin, margin,
                     inset_width, inset_height)
        painter.setPen(QPen(QColor("#dce4e9"), max(1.0, 2 * unit)))
        painter.setBrush(QColor(9, 12, 15, 235))
        painter.drawRoundedRect(box, 10 * unit, 10 * unit)
        # Keep measurement text outside the FFT image.  An older floating
        # black label covered the upper edge of a first-order polygon and made
        # an otherwise closed frame look incomplete.
        title_h = (64 if not selected else 44) * unit
        painter.setPen(QColor("white"))
        painter.setFont(self._annotation_font(16, QFont.Weight.Bold))
        title = ("FFT · %d precisely fitted clear spots" % len(
                    assets.get("selected_spots", [])) if selected else
                 "Measured FFT · a-matched first-order peaks")
        title_rect = QRectF(box.x() + 16 * unit,
                            box.y() + 3 * unit,
                            box.width() - 24 * unit,
                            (34 if not selected else 38) * unit)
        painter.drawText(title_rect,
                         Qt.AlignmentFlag.AlignLeft |
                         Qt.AlignmentFlag.AlignVCenter, title)
        if (not selected and bilayer and
                metrics.get("fft_twist_reliable")):
            painter.setPen(QColor("#cbd5dc"))
            painter.setFont(self._annotation_font(13, QFont.Weight.Bold))
            painter.drawText(
                QRectF(box.x() + 16 * unit, box.y() + 34 * unit,
                       box.width() - 32 * unit, 24 * unit),
                Qt.AlignmentFlag.AlignLeft |
                Qt.AlignmentFlag.AlignVCenter,
                "FFT-derived twist = %s" % self._fmt(
                    metrics.get("fft_twist_deg"), "°"))
        if selected:
            image_bounds = QRectF(box.x() + 15 * unit, box.y() + title_h,
                                  box.width() - 30 * unit,
                                  box.height() - title_h - 17 * unit)
            source_rect = QRectF(fft.rect())
        else:
            image_bounds = QRectF(box.x() + 14 * unit, box.y() + title_h,
                                  box.width() - 28 * unit, 230 * unit)
            # Keep the approved central field of view, but use equal
            # *fractional frequency* spans in x and y.  The former 550x380
            # reference crop had unequal spans and made Square FFTs appear
            # stretched for non-square source TEM images.
            fraction = max(550.0 / 1350.0, 380.0 / 1040.0)
            crop_width = fft.width() * fraction
            crop_height = fft.height() * fraction
            source_rect = QRectF((fft.width() - crop_width) / 2.0,
                                 (fft.height() - crop_height) / 2.0,
                                 crop_width, crop_height)
        inset = self._fit_fft_rect(image_bounds, fft, source_rect)
        painter.drawImage(inset, fft, source_rect)
        sx = inset.width() / source_rect.width()
        sy = inset.height() / source_rect.height()
        painter.save()
        painter.translate(inset.x() - source_rect.x() * sx,
                          inset.y() - source_rect.y() * sy)
        painter.scale(sx, sy)
        if selected:
            for spot in assets.get("selected_spots", []):
                role = spot.get("lattice_role")
                color = (QColor("#2bd4f0") if role in (
                            "twin_orientation", "square_twin_reflection",
                            "square_layer_2_reflection") else
                         QColor("#70f08c") if role in (
                            "secondary_a", "square_supercell_first_order",
                            "square_supercell_reflection",
                            "square_helix_reflection",
                            "square_layer_reflection") else
                         QColor(255, 214, 64))
                painter.setPen(QPen(color, max(1.0, 1.2 / sx)))
                rx = float(spot.get("rx", spot.get("radius_x", 4)))
                ry = float(spot.get("ry", spot.get("radius_y", 4)))
                painter.save()
                painter.translate(float(spot["x"]), float(spot["y"]))
                painter.rotate(float(spot.get("angle") or 0.0))
                painter.drawEllipse(QPointF(0, 0), rx, ry)
                painter.restore()
        else:
            quadrilaterals = assets.get("first_order_quadrilaterals", [])
            if result.get("analysis_kind") == "bilayer":
                lattice_fft = result.get("lattice_fft") or {}
                bilayer_valid = (
                    assets.get("first_order_bilayer_valid")
                    if "first_order_bilayer_valid" in assets else
                    lattice_fft.get("valid"))
                if not bilayer_valid:
                    # The numeric FFT validator has rejected two physical
                    # layers.  A second radial shell can be the same Square
                    # lattice's 45-degree diagonal reflection; never color it
                    # as another layer merely for display.
                    quadrilaterals = quadrilaterals[:1]
            mixed_families = assets.get("mixed_lattice_families") or []
            if mixed_families:
                center_value = assets.get("center") or [
                    fft.width()/2.0, fft.height()/2.0]
                quadrilaterals = []
                for family in mixed_families:
                    members = [dict(point) for point in
                               (family.get("peaks") or [])]
                    members.sort(key=lambda point: math.atan2(
                        float(point["y"])-float(center_value[1]),
                        float(point["x"])-float(center_value[0])))
                    if len(members) >= 4:
                        quadrilaterals.append(members)
            colors = tuple(QColor(value) for value in (
                "#00aee8", "#f57c00", "#63b74f", "#b24ca6",
                "#d6ad00", "#6d71d9"))
            center = assets.get("center") or [fft.width() / 2,
                                               fft.height() / 2]
            display_limit = (2 if result.get("analysis_kind") == "bilayer"
                             else 6)
            for color, polygon in zip(colors,
                                      quadrilaterals[:display_limit]):
                # A first-order family is a closed reciprocal-lattice frame.
                # Restore its cyclic order around the FFT origin, then close
                # one continuous painter path.  The pen remains dashed, while
                # NoBrush keeps the area inside the frame fully transparent.
                ordered = sorted(
                    polygon,
                    key=lambda point: math.atan2(
                        float(point["y"] if isinstance(point, dict)
                              else point[1]) - float(center[1]),
                        float(point["x"] if isinstance(point, dict)
                              else point[0]) - float(center[0])))
                painter.setPen(QPen(color, max(1.0, 3.0 / sx),
                                    Qt.PenStyle.DashLine))
                painter.setBrush(Qt.BrushStyle.NoBrush)
                points = [
                    QPointF(float(point["x"] if isinstance(point, dict)
                                  else point[0]),
                            float(point["y"] if isinstance(point, dict)
                                  else point[1]))
                    for point in ordered]
                if points:
                    outline = QPainterPath()
                    outline.moveTo(points[0])
                    for point in points[1:]:
                        outline.lineTo(point)
                    outline.closeSubpath()
                    painter.drawPath(outline)
            pore = result.get("pore_lattice") or {}
            pore_points = (pore.get("refined_peaks") or
                           pore.get("peaks") or [])
            if pore.get("valid") and len(pore_points) >= 4:
                pore_color = QColor("#70f08c")
                painter.setPen(QPen(pore_color, max(1.0, 2.5 / sx),
                                    Qt.PenStyle.DashLine))
                painter.drawPolygon(QPolygonF([
                    QPointF(float(point["x"]), float(point["y"]))
                    for point in pore_points[:4]]))
                for point in pore_points[:4]:
                    painter.drawEllipse(QPointF(float(point["x"]),
                                                float(point["y"])),
                                        5.0 / sx, 5.0 / sy)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            for color, polygon in zip(colors,
                                      quadrilaterals[:display_limit]):
                if not polygon:
                    continue
                point = polygon[0]
                px = float(point.get("x") if isinstance(point, dict) else point[0])
                py = float(point.get("y") if isinstance(point, dict) else point[1])
                painter.setPen(QPen(color, max(1.0, 2.5 / sx),
                                    Qt.PenStyle.DashLine))
                painter.drawLine(QPointF(float(center[0]), float(center[1])),
                                 QPointF(px, py))
                painter.drawEllipse(QPointF(px, py), 6.0 / sx, 6.0 / sy)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor("white"))
            painter.drawEllipse(QPointF(float(center[0]), float(center[1])),
                                5.0 / sx, 5.0 / sy)
        painter.restore()
        if not selected:
            painter.setPen(QColor("white"))
            painter.setFont(self._annotation_font(13))
            text_y = box.y() + 313 * unit
            if bilayer and metrics.get("fft_twist_reliable"):
                painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                    "FFT-derived twist: %s · FFT-derived moiré period: %s" % (
                        self._fmt(metrics.get("fft_twist_deg"), "°"),
                        self._fmt(metrics.get("fft_period_nm"), " nm")))
                text_y += 24 * unit
                if metrics.get("tem_period_reliable"):
                    painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                        "TEM-derived twist: %s · TEM-derived moiré period: %s" % (
                            self._fmt(metrics.get("tem_twist_deg"), "°"),
                            self._fmt(metrics.get("tem_period_nm"), " nm")))
                    text_y += 24 * unit
            elif bilayer:
                painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                    "FFT-derived twist unavailable · two peak sets too close")
                text_y += 24 * unit
                painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                    "TEM-derived twist unavailable · <2 reliable Moiré units")
                text_y += 24 * unit
            entries = metrics.get("fft_lattices") or []
            lattice_text = "; ".join(
                "%s a = %s" % (item.get("symmetry") or "Lattice",
                                 self._fmt(item.get("a_nm"), " nm"))
                for item in entries) or "a = —"
            if (not bilayer and metrics.get("pore_a_nm") is not None):
                lattice_text = "helix a = %s · pore a = %s" % (
                    self._fmt(metrics.get("fft_a_nm"), " nm"),
                    self._fmt(metrics.get("pore_a_nm"), " nm"))
            if (bilayer and metrics.get("fft_twist_reliable") and
                    metrics.get("tem_period_reliable")):
                tem_twist = metrics.get("tem_twist_deg")
                fft_twist = metrics.get("fft_twist_deg")
                mean_twist = ((tem_twist + fft_twist) / 2
                              if tem_twist is not None and fft_twist is not None
                              else None)
                painter.setPen(QColor("#cbd5dc"))
                painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                    "Mean a = %s · mean twist = %s" % (
                        self._fmt(metrics.get("mean_a_nm"), " nm"),
                        self._fmt(mean_twist, "°")))
                text_y += 24 * unit
                painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                    "Δtwist = %s · Δperiod = %s" % (
                        self._deviation_text(fft_twist, tem_twist, "°"),
                        self._deviation_text(metrics.get("fft_period_nm"),
                                             metrics.get("tem_period_nm"), " nm")))
            else:
                if not bilayer and spatial_domains:
                    painter.setFont(self._annotation_font(11))
                    column_width = (box.width()-32*unit)/2.0
                    rows_per_column = int(math.ceil(
                        len(spatial_domains)/2.0))
                    pixel_nm = metrics.get("pixel_size_nm")
                    for index, domain in enumerate(spatial_domains):
                        column = index//rows_per_column
                        row = index % rows_per_column
                        x = box.x()+16*unit+column*column_width
                        y = text_y+row*38*unit
                        lattice_px = domain.get("lattice_constant_px")
                        lattice_value = (
                            self._fmt(float(lattice_px)*pixel_nm, " nm")
                            if lattice_px is not None and pixel_nm else
                            self._fmt(lattice_px, " px"))
                        headline = "D%d %s · a %s · %.1f° · %.1f%%" % (
                            int(domain.get("domain_id", index+1)),
                            str(domain.get("symmetry") or "Square"),
                            lattice_value,
                            float(domain.get("orientation_deg") or 0.0),
                            100.0*float(domain.get("area_fraction", 0.0)))
                        gaps = self._domain_inter_axis_angles(domain)
                        axis_name = ("actual two-axis angles "
                                     if str(domain.get("symmetry") or
                                            "Square") == "Square" else
                                     "actual tri-axis angles ")
                        details = (axis_name + "/".join(
                            "%.1f°" % float(value) for value in gaps)
                                   if gaps else axis_name + "—")
                        painter.setPen(QColor("white"))
                        painter.drawText(QPointF(x, y), headline)
                        painter.setPen(QColor("#aebcc5"))
                        painter.drawText(QPointF(x, y+16*unit), details)
                else:
                    painter.drawText(QPointF(box.x() + 16 * unit, text_y),
                                     lattice_text)

    def _draw_scale_bar(self, painter, record):
        metrics = self._render_metrics(record)
        pixel_nm = metrics.get("pixel_size_nm")
        if not pixel_nm:
            return
        source = QImage(record["source"])
        unit, text_unit = self._uploaded_annotation_units(source)
        length_nm = 100.0
        length_px = length_nm / pixel_nm
        # Keep the bar visible in small fields while preserving its stated
        # physical length whenever the field is large enough.
        if length_px > source.width() * 0.36:
            length_nm = 50.0
            length_px = length_nm / pixel_nm
        x = 48 * unit
        y = source.height() - 85 * unit
        painter.setPen(QPen(QColor("white"), max(3.0, 7 * unit)))
        painter.drawLine(QPointF(x, y), QPointF(x + length_px, y))
        painter.setFont(self._annotation_font(14, QFont.Weight.Bold))
        painter.drawText(QRectF(x, y - 38 * unit, length_px, 30 * unit),
                         Qt.AlignmentFlag.AlignCenter,
                         "%g nm" % length_nm)

    def _draw_moire_cells(self, painter, record, draw_cells=True):
        if not self._render_metrics(record).get("tem_period_reliable"):
            return
        real = record["result"].get("moire_real_space") or {}
        if not real.get("valid"):
            return
        basis = real.get("basis_vectors_px") or []
        if draw_cells and len(basis) == 2:
            first, second = basis
            source = QImage(record["source"])
            unit, unused_text_unit = self._uploaded_annotation_units(source)
            painter.setPen(QPen(QColor(77, 255, 126, 225),
                                max(2.0, 3.0 * unit),
                                Qt.PenStyle.SolidLine))
            painter.setBrush(Qt.BrushStyle.NoBrush)
            for center in real.get("centers", []):
                cx, cy = float(center["x"]), float(center["y"])
                painter.drawPolygon(QPolygonF([
                    QPointF(cx-first[0]/2-second[0]/2,
                            cy-first[1]/2-second[1]/2),
                    QPointF(cx+first[0]/2-second[0]/2,
                            cy+first[1]/2-second[1]/2),
                    QPointF(cx+first[0]/2+second[0]/2,
                            cy+first[1]/2+second[1]/2),
                    QPointF(cx-first[0]/2+second[0]/2,
                            cy-first[1]/2+second[1]/2)]))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(112, 240, 140, 225))
        for center in real.get("centers", []):
            painter.drawEllipse(QPointF(float(center["x"]),
                                        float(center["y"])), 3.5, 3.5)
        pair = real.get("representative_pair") or {}
        if pair:
            first, second = pair.get("first") or {}, pair.get("second") or {}
            p1 = QPointF(float(first.get("x", 0)), float(first.get("y", 0)))
            p2 = QPointF(float(second.get("x", 0)), float(second.get("y", 0)))
            source = QImage(record["source"])
            unit, text_unit = self._uploaded_annotation_units(source)
            painter.setPen(QPen(QColor("#ffd234"), max(3.0, 5.0 * unit),
                                Qt.PenStyle.SolidLine))
            painter.setBrush(QColor("#ffd234"))
            painter.drawLine(p1, p2)
            painter.drawEllipse(p1, 5 * unit, 5 * unit)
            painter.drawEllipse(p2, 5 * unit, 5 * unit)
            metrics = self._render_metrics(record)
            if metrics.get("tem_period_nm"):
                text = "TEM-derived moiré period = %.3f nm" % metrics[
                    "tem_period_nm"]
                midpoint = QPointF((p1.x() + p2.x()) / 2,
                                   (p1.y() + p2.y()) / 2)
                box = QRectF(midpoint.x() - 142.5 * unit,
                             midpoint.y() + 36 * unit,
                             285 * unit, 45 * unit)
                painter.fillRect(box, QColor(16, 22, 26, 225))
                painter.setPen(QColor("white"))
                painter.setFont(self._annotation_font(
                    14, QFont.Weight.Bold))
                painter.drawText(box, Qt.AlignmentFlag.AlignCenter, text)

    def _draw_reconstruction_annotations(self, painter, record):
        source = QImage(record["source"])
        metrics = self._render_metrics(record)
        unit, text_unit = self._uploaded_annotation_units(source)
        entries = metrics.get("fft_lattices") or []
        pore = record["result"].get("pore_lattice") or {}
        twin = (record["result"].get("lattice_fft") or {}).get("twin") or {}
        has_twin = bool(twin.get("valid"))
        has_supercell = bool(pore.get("valid") and
                             metrics.get("pore_a_nm") is not None)
        bilayer = record["result"].get("analysis_kind") == "bilayer"
        entry_symmetries = {
            str(item.get("symmetry") or "Lattice") for item in entries}
        shared_bilayer_a = bool(
            bilayer and len(entries) > 1 and len(entry_symmetries) == 1)
        if has_twin:
            symmetry = (entries[0].get("symmetry") if entries else "Square")
            lattice_name = "%s twin domains" % symmetry
        elif has_supercell:
            symmetry = (entries[0].get("symmetry") if entries else "Square")
            lattice_name = "%s helix + supercell" % symmetry
        elif len(entries) > 1:
            lattice_name = ("–".join(
                item.get("symmetry") or "Lattice" for item in entries[:4]) +
                (" bilayer" if bilayer else " spatial domains"))
        else:
            symmetry = (entries[0].get("symmetry") if entries else "Square")
            lattice_name = ("%s–%s bilayer" % (symmetry, symmetry)
                            if bilayer
                            else "%s lattice" % symmetry)
        values = [(lattice_name, True)]
        if has_twin:
            orientations = twin.get("orientations_deg") or []
            values += [
                ("a = %s" % self._fmt(metrics.get("fft_a_nm"), " nm"), False),
                ("orientation 1 = %s" % self._fmt(
                    orientations[0] if len(orientations) > 0 else None, "°"),
                 False),
                ("orientation 2 = %s" % self._fmt(
                    orientations[1] if len(orientations) > 1 else None, "°"),
                 False),
                ("twin Δ = %s" % self._fmt(
                    twin.get("relative_orientation_deg"), "°"), False)]
        elif has_supercell:
            values += [
                ("helix a = %s" % self._fmt(
                    metrics.get("fft_a_nm"), " nm"), False),
                ("supercell a = %s" % self._fmt(
                    metrics.get("pore_a_nm"), " nm"), False),
                ("relation = %s× helix lattice" % (
                    pore.get("repeat_multiple") or "—"), False)]
        elif len(entries) > 1 and not shared_bilayer_a:
            values.extend(("a%d = %s" % (
                index + 1, self._fmt(item.get("a_nm"), " nm")), False)
                          for index, item in enumerate(entries[:2]))
        else:
            values.append(("a = %s" % self._fmt(
                metrics.get("mean_a_nm"), " nm"), False))
        if bilayer and metrics.get("fft_twist_reliable"):
            values += [
                ("relative twist = %s" % self._fmt(
                    metrics.get("reconstruction_twist_deg"), "°"), False),
                ("Moiré period = %s" % self._fmt(
                    metrics.get("reconstruction_period_nm"), " nm"), False)]
        box_height = max(158.0, 22.0 + 34.0 * len(values)) * unit
        box = QRectF(972 * source.width() / 1350.0,
                     28 * source.height() / 1040.0,
                     350 * unit, box_height)
        painter.setPen(QPen(QColor("#dce4e9"), max(1.0, 2 * unit)))
        painter.setBrush(QColor(9, 12, 15, 230))
        painter.drawRoundedRect(box, 10 * unit, 10 * unit)
        painter.setPen(QColor("white"))
        for index, (text, bold) in enumerate(values):
            painter.setFont(self._annotation_font(
                16 if bold else 14,
                QFont.Weight.Bold if bold else QFont.Weight.Normal))
            painter.drawText(QPointF(box.x() + 18 * unit,
                                     box.y() + (34 + index * 34) * unit), text)
        if bilayer and metrics.get("fft_twist_reliable"):
            cx = 675 * source.width() / 1350.0
            cy = 175 * source.height() / 1040.0
            half = 60 * unit
            angles = self._reconstruction_layer_angles(record, metrics)
            colors = (QColor("#ffbd3f"), QColor("#2bd4f0"))
            for angle, color in zip(angles, colors):
                # Match the TEM/FFT absolute orientation in image coordinates.
                radians = math.radians(-angle)
                ca, sa = math.cos(radians), math.sin(radians)
                painter.setPen(QPen(color, max(1.5, 3 * unit),
                                    Qt.PenStyle.DashLine))
                for x1, y1, x2, y2 in ((-half, 0, half, 0),
                                       (0, -half, 0, half)):
                    painter.drawLine(
                        QPointF(cx + ca*x1 - sa*y1, cy + sa*x1 + ca*y1),
                        QPointF(cx + ca*x2 - sa*y2, cy + sa*x2 + ca*y2))
            period = float(record["result"].get("moire_period_px") or
                           record["result"].get(
                               "fft_predicted_moire_period_px") or 0.0)
            if period > 0:
                half = period / 2.0
                painter.setPen(QPen(QColor("#70f08c"), max(1.5, 3 * unit)))
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawPolygon(QPolygonF([
                    QPointF(source.width()/2-half, source.height()/2-half),
                    QPointF(source.width()/2+half, source.height()/2-half),
                    QPointF(source.width()/2+half, source.height()/2+half),
                    QPointF(source.width()/2-half, source.height()/2+half)]))
        self._draw_scale_bar(painter, record)

    def _draw_tem_twist_marker(self, painter, record):
        metrics = self._render_metrics(record)
        if not metrics.get("tem_period_reliable"):
            return
        source = QImage(record["source"])
        unit, text_unit = self._uploaded_annotation_units(source)
        cx, cy = 310 * source.width() / 1350.0, 180 * source.height() / 1040.0
        half = 60 * unit
        angles = [0.0, float(metrics.get("fft_twist_deg") or
                             metrics.get("tem_twist_deg") or 0.0)]
        colors = (QColor("#ffbd3f"), QColor("#2bd4f0"))
        for angle, color in zip(angles, colors):
            radians = math.radians(angle)
            painter.setPen(QPen(color, max(1.5, 3 * unit),
                                Qt.PenStyle.DashLine))
            ca, sa = math.cos(radians), math.sin(radians)
            # Rotate an exact 120-by-120 orthogonal cross around its center.
            for x1, y1, x2, y2 in ((-half, 0, half, 0),
                                   (0, -half, 0, half)):
                painter.drawLine(
                    QPointF(cx + ca * x1 - sa * y1,
                            cy + sa * x1 + ca * y1),
                    QPointF(cx + ca * x2 - sa * y2,
                            cy + sa * x2 + ca * y2))
        text = "TEM-derived twist = %s" % self._fmt(
            metrics.get("tem_twist_deg"), "°")
        painter.setFont(self._annotation_font(14, QFont.Weight.Bold))
        box = QRectF(190 * source.width() / 1350.0,
                     250 * source.height() / 1040.0,
                     300 * unit, 48 * unit)
        painter.fillRect(box, QColor(9, 12, 15, 220))
        painter.setPen(QColor("white"))
        painter.drawText(box, Qt.AlignmentFlag.AlignCenter, text)

    @staticmethod
    def _boundary_vector_paths(segments):
        """Rejoin adjacent measured segments into editable vector polylines."""
        adjacency, unused = {}, set()
        for segment in segments or []:
            if len(segment) != 4:
                continue
            first = (round(float(segment[0]), 3),
                     round(float(segment[1]), 3))
            second = (round(float(segment[2]), 3),
                      round(float(segment[3]), 3))
            if first == second:
                continue
            edge = tuple(sorted((first, second)))
            unused.add(edge)
            adjacency.setdefault(first, set()).add(second)
            adjacency.setdefault(second, set()).add(first)
        paths = []
        while unused:
            degree = {}
            for edge in unused:
                for point in edge:
                    degree[point] = degree.get(point, 0)+1
            endpoints = [point for point, value in degree.items()
                         if value == 1]
            start = endpoints[0] if endpoints else next(iter(unused))[0]
            points, current, previous = [start], start, None
            while True:
                options = [neighbor for neighbor in adjacency.get(current, ())
                           if tuple(sorted((current, neighbor))) in unused]
                if not options:
                    break
                if previous is not None and len(options) > 1:
                    old_angle = math.atan2(current[1]-previous[1],
                                           current[0]-previous[0])
                    options.sort(key=lambda point: abs(math.atan2(
                        math.sin(math.atan2(point[1]-current[1],
                                            point[0]-current[0])-old_angle),
                        math.cos(math.atan2(point[1]-current[1],
                                            point[0]-current[0])-old_angle))))
                following = options[0]
                unused.remove(tuple(sorted((current, following))))
                points.append(following)
                previous, current = current, following
                if current == start:
                    break
            if len(points) < 2:
                continue
            path = QPainterPath(QPointF(*points[0]))
            for point in points[1:]:
                path.lineTo(QPointF(*point))
            paths.append(path)
        return paths

    def _draw_sample_boundary(self, painter, record):
        """Draw the physical specimen/background edge in either mode."""
        lattice = record["result"].get("lattice_fft") or {}
        spatial = (record["result"].get("orientation_domains") or
                   lattice.get("orientation_domains") or {})
        sample_boundary = spatial.get("sample_boundary") or {}
        sample_segments = sample_boundary.get("segments") or []
        if not sample_segments:
            return
        source = QImage(record["source"])
        unit, text_unit = self._uploaded_annotation_units(source)
        sample_color = QColor("#20b486")
        # A specimen boundary is a closed QPainterPath.  QPainter fills a
        # closed path with the currently inherited brush unless it is reset
        # explicitly; after an information label that brush is translucent
        # white and would therefore cover the entire specimen interior.
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for path in self._boundary_vector_paths(sample_segments):
            painter.setPen(QPen(QColor(255, 255, 255, 220),
                                max(2.0, 5.5 * unit)))
            painter.drawPath(path)
            painter.setPen(QPen(sample_color, max(1.5, 3.0 * unit),
                                Qt.PenStyle.DashLine))
            painter.drawPath(path)
        middle = sample_segments[len(sample_segments)//2]
        label_x = sample_boundary.get("label_x")
        label_y = sample_boundary.get("label_y")
        x = max(8.0, min(source.width()-230*unit,
            (float(label_x) if label_x is not None else
             (float(middle[0])+float(middle[2]))/2) + 10*unit))
        y = max(62.0, min(source.height()-38*unit,
            (float(label_y) if label_y is not None else
             (float(middle[1])+float(middle[3]))/2) - 34*unit))
        box = QRectF(x, y, 220*unit, 34*unit)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(255, 255, 255, 218))
        painter.drawRoundedRect(box, 5*unit, 5*unit)
        painter.setPen(sample_color)
        painter.setFont(self._annotation_font(14, QFont.Weight.Bold))
        painter.drawText(box, Qt.AlignmentFlag.AlignCenter,
                         "sample boundary")

    def _sample_analysis_path(self, record):
        """Return the closed vector path defining the analyzed specimen."""
        lattice = record["result"].get("lattice_fft") or {}
        spatial = (record["result"].get("orientation_domains") or
                   lattice.get("orientation_domains") or {})
        sample = spatial.get("sample_boundary") or {}
        segments = sample.get("segments") or []
        if not segments:
            return QPainterPath()
        source = QImage(record["source"])
        width = max(1.0, float(source.width()-1))
        height = max(1.0, float(source.height()-1))
        perimeter = 2.0*(width+height)
        domain_markers = [QPointF(float(domain.get("marker_x", width/2.0)),
                                  float(domain.get("marker_y", height/2.0)))
                          for domain in (spatial.get("domains") or [])]

        def snap_to_frame(point):
            x, y = float(point.x()), float(point.y())
            candidates = [
                (abs(y), QPointF(min(max(x, 0.0), width), 0.0),
                 min(max(x, 0.0), width)),
                (abs(width-x), QPointF(width, min(max(y, 0.0), height)),
                 width+min(max(y, 0.0), height)),
                (abs(height-y), QPointF(min(max(x, 0.0), width), height),
                 width+height+(width-min(max(x, 0.0), width))),
                (abs(x), QPointF(0.0, min(max(y, 0.0), height)),
                 2.0*width+height+(height-min(max(y, 0.0), height))),
            ]
            unused_distance, snapped, coordinate = min(
                candidates, key=lambda item: item[0])
            return snapped, coordinate % perimeter

        def clockwise_route(first, first_t, second, second_t):
            target = second_t
            if target <= first_t+1e-9:
                target += perimeter
            corners = [
                (0.0, QPointF(0.0, 0.0)),
                (width, QPointF(width, 0.0)),
                (width+height, QPointF(width, height)),
                (2.0*width+height, QPointF(0.0, height)),
                (perimeter, QPointF(0.0, 0.0)),
            ]
            route = [first]
            for offset in (0.0, perimeter):
                for value, corner in corners:
                    absolute = value+offset
                    if first_t+1e-9 < absolute < target-1e-9:
                        route.append(corner)
            route.append(second)
            return route

        def polygon_path(points):
            output = QPainterPath(points[0])
            for point in points[1:]:
                output.lineTo(point)
            output.closeSubpath()
            return output

        combined = QPainterPath()
        combined.setFillRule(Qt.FillRule.OddEvenFill)
        for path in self._boundary_vector_paths(segments):
            points = [QPointF(path.elementAt(index).x,
                              path.elementAt(index).y)
                      for index in range(path.elementCount())]
            if len(points) < 2:
                continue
            closed = math.hypot(points[0].x()-points[-1].x(),
                                points[0].y()-points[-1].y()) <= 3.0
            if closed:
                path.closeSubpath()
                combined.addPath(path)
                continue
            start, start_t = snap_to_frame(points[0])
            end, end_t = snap_to_frame(points[-1])
            chain = [start]+points[1:-1]+[end]
            first = polygon_path(
                chain+clockwise_route(end, end_t, start, start_t)[1:])
            reverse_route = list(reversed(
                clockwise_route(start, start_t, end, end_t)))
            second = polygon_path(chain+reverse_route[1:])
            first_hits = sum(1 for marker in domain_markers
                             if first.contains(marker))
            second_hits = sum(1 for marker in domain_markers
                              if second.contains(marker))
            if first_hits != second_hits:
                combined.addPath(first if first_hits > second_hits else second)
            else:
                expected = float(spatial.get("sample_cell_fraction", .5))
                first_area = abs(first.boundingRect().width() *
                                 first.boundingRect().height())
                second_area = abs(second.boundingRect().width() *
                                  second.boundingRect().height())
                choose_first = ((first_area >= second_area) ==
                                (expected >= .5))
                combined.addPath(first if choose_first else second)
        return combined

    def _draw_domain_area_overlay(self, painter, record):
        """Tint measured real-space domains for single-layer analysis."""
        lattice = record["result"].get("lattice_fft") or {}
        spatial = (record["result"].get("orientation_domains") or
                   lattice.get("orientation_domains") or {})
        domains = spatial.get("domains") or []
        source = QImage(record["source"])
        sample_path = self._sample_analysis_path(record)
        has_sample = not sample_path.isEmpty()
        if not domains and not has_sample:
            return
        colors = [QColor(value) for value in (
            "#00aee8", "#f57c00", "#63b74f", "#b24ca6",
            "#d6ad00", "#6d71d9")]

        if has_sample:
            full_path = QPainterPath()
            full_path.addRect(QRectF(source.rect()))
            outside_path = full_path.subtracted(sample_path)
            unit, unused_text_unit = self._uploaded_annotation_units(source)
            spacing = max(16.0, 25.0*unit)
            painter.save()
            painter.setClipPath(outside_path)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(QColor(205, 220, 230, 175),
                                max(1.2, 1.45*unit)))
            x = 0.0
            while x <= source.width():
                painter.drawLine(QPointF(x, 0.0),
                                 QPointF(x, float(source.height())))
                x += spacing
            y = 0.0
            while y <= source.height():
                painter.drawLine(QPointF(0.0, y),
                                 QPointF(float(source.width()), y))
                y += spacing
            painter.restore()

        painter.save()
        painter.setPen(Qt.PenStyle.NoPen)
        if has_sample:
            painter.setClipPath(sample_path)
        if has_sample and len(domains) <= 1:
            color = QColor(colors[0])
            color.setAlpha(38)
            painter.setBrush(color)
            painter.drawPath(sample_path)

        ordered_domains = sorted(
            enumerate(domains),
            key=lambda item: bool(item[1].get("manual_boundary_exact")))
        for domain_index, domain in ordered_domains:
            if len(domains) <= 1:
                continue
            polygons = domain.get("area_polygons") or []
            runs = domain.get("area_runs") or []
            if not polygons and not runs:
                continue
            color = QColor(colors[domain_index % len(colors)])
            color.setAlpha(38)
            painter.setBrush(color)
            if polygons:
                for polygon in polygons:
                    points = [QPointF(float(point[0]), float(point[1]))
                              for point in polygon if len(point) >= 2]
                    if len(points) >= 3:
                        painter.drawPolygon(QPolygonF(points))
            else:
                for x0, y0, x1, y1 in runs:
                    painter.drawRect(QRectF(
                        float(x0), float(y0),
                        max(0.0, float(x1)-float(x0)),
                        max(0.0, float(y1)-float(y0))))
        painter.restore()

    def _draw_single_orientation_domains(self, painter, record):
        """Draw each single-layer orientation in its measured real-space domain."""
        lattice = record["result"].get("lattice_fft") or {}
        spatial = (record["result"].get("orientation_domains") or
                   lattice.get("orientation_domains") or {})
        domains = spatial.get("domains") or []
        boundaries = spatial.get("boundaries") or []
        sample_boundary = spatial.get("sample_boundary") or {}
        source = QImage(record["source"])
        if not domains:
            # Backward-compatible fallback for older cached analyses.
            layers = lattice.get("layers") or []
            if not layers:
                return
            domains = [{
                "domain_id": index + 1,
                "orientation_deg": layer.get("orientation_deg", 0.0),
                "symmetry": layer.get("symmetry", "Square"),
                "lattice_constant_px": layer.get("lattice_constant_px"),
                "reciprocal_axis_angles_deg": layer.get(
                    "reciprocal_axis_angles_deg") or [],
                "marker_x": source.width() * (
                    .42 if len(layers) == 1 else (.28 + .32*index)),
                "marker_y": source.height() * .34,
            } for index, layer in enumerate(layers)]
        unit, text_unit = self._uploaded_annotation_units(source)
        colors = [QColor(value) for value in (
            "#00aee8", "#f57c00", "#63b74f", "#b24ca6",
            "#d6ad00", "#6d71d9")]

        boundary_color = QColor("#f53390")
        placed_boundary_labels = []
        # Boundaries are vector strokes only.  Keep both closed specimen
        # contours and any closed crystallographic-domain contours unfilled.
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for boundary_index, boundary in enumerate(boundaries, 1):
            segments = boundary.get("segments") or []
            painter.setBrush(Qt.BrushStyle.NoBrush)
            for path in self._boundary_vector_paths(segments):
                painter.setPen(QPen(QColor(255, 255, 255, 220),
                                    max(2.0, 5.5 * unit)))
                painter.drawPath(path)
                painter.setPen(QPen(boundary_color, max(1.5, 3.0 * unit),
                                    Qt.PenStyle.DashLine))
                painter.drawPath(path)
            if segments:
                middle = segments[len(segments)//2]
                label_x = boundary.get("label_x")
                label_y = boundary.get("label_y")
                x = max(8.0, min(source.width()-250*unit,
                    (float(label_x) if label_x is not None else
                     (float(middle[0])+float(middle[2]))/2) + 10*unit))
                y = max(62.0, min(source.height()-38*unit,
                    (float(label_y) if label_y is not None else
                     (float(middle[1])+float(middle[3]))/2) - 34*unit))
                text = ("domain boundary" if len(boundaries) == 1 else
                        "domain boundary %d" % boundary_index)
                box = QRectF(x, y, 240*unit, 34*unit)
                attempts = 0
                while (any(box.intersects(old)
                           for old in placed_boundary_labels) and
                       attempts < 8):
                    next_y = box.y()+42*unit
                    if next_y+box.height() > source.height()-8*unit:
                        next_y = max(62.0, box.y()-84*unit)
                    box.moveTop(next_y)
                    attempts += 1
                placed_boundary_labels.append(QRectF(box))
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(QColor(255, 255, 255, 218))
                painter.drawRoundedRect(box, 5*unit, 5*unit)
                painter.setPen(boundary_color)
                painter.setFont(self._annotation_font(
                    14, QFont.Weight.Bold))
                painter.drawText(box, Qt.AlignmentFlag.AlignCenter, text)

        # The specimen edge is not a crystallographic domain boundary.  Draw
        # it with its own colour and label while keeping the same precise,
        # editable vector-line treatment in the SVG export.
        sample_segments = sample_boundary.get("segments") or []
        sample_color = QColor("#20b486")
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for path in self._boundary_vector_paths(sample_segments):
            painter.setPen(QPen(QColor(255, 255, 255, 220),
                                max(2.0, 5.5 * unit)))
            painter.drawPath(path)
            painter.setPen(QPen(sample_color, max(1.5, 3.0 * unit),
                                Qt.PenStyle.DashLine))
            painter.drawPath(path)
        if sample_segments:
            middle = sample_segments[len(sample_segments)//2]
            label_x = sample_boundary.get("label_x")
            label_y = sample_boundary.get("label_y")
            x = max(8.0, min(source.width()-230*unit,
                (float(label_x) if label_x is not None else
                 (float(middle[0])+float(middle[2]))/2) + 10*unit))
            y = max(62.0, min(source.height()-38*unit,
                (float(label_y) if label_y is not None else
                 (float(middle[1])+float(middle[3]))/2) - 34*unit))
            box = QRectF(x, y, 220*unit, 34*unit)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(255, 255, 255, 218))
            painter.drawRoundedRect(box, 5*unit, 5*unit)
            painter.setPen(sample_color)
            painter.setFont(self._annotation_font(14, QFont.Weight.Bold))
            painter.drawText(box, Qt.AlignmentFlag.AlignCenter,
                             "sample boundary")

        half = 52 * unit
        multiple = len(domains) > 1
        pixel_size_nm = self._render_metrics(record).get("pixel_size_nm")
        for domain_index, domain in enumerate(domains):
            angle = float(domain.get("orientation_deg") or 0.0)
            symmetry = str(domain.get("symmetry") or "Square")
            color = colors[domain_index % len(colors)]
            cx = float(domain.get("marker_x", source.width()/2.0))
            cy = float(domain.get("marker_y", source.height()/2.0))
            painter.setPen(QPen(color, max(1.5, 3*unit),
                                Qt.PenStyle.DashLine))
            reciprocal_angles = domain.get(
                "reciprocal_axis_angles_deg") or []
            if symmetry == "Square":
                marker_angles = [angle, angle+90.0]
            elif len(reciprocal_angles) >= 3:
                # Preserve the measured non-60-degree Kagome distortion in
                # the marker instead of drawing an idealized hexagon.
                marker_angles = [float(value)-30.0
                                 for value in reciprocal_angles[:3]]
            else:
                marker_angles = [angle, angle+60.0, angle+120.0]
            for marker_angle in marker_angles:
                marker_radians = math.radians(-marker_angle)
                marker_ca, marker_sa = (math.cos(marker_radians),
                                        math.sin(marker_radians))
                painter.drawLine(
                    QPointF(cx-marker_ca*half, cy-marker_sa*half),
                    QPointF(cx+marker_ca*half, cy+marker_sa*half))
            domain_id = int(domain.get("domain_id", domain_index+1))
            lattice_px = domain.get("lattice_constant_px")
            lattice_nm = (float(lattice_px)*pixel_size_nm
                          if lattice_px and pixel_size_nm else None)
            a_value = (("%.2f nm" % lattice_nm) if lattice_nm else
                       (("%.1f px" % float(lattice_px))
                        if lattice_px else "—"))
            area = 100.0*float(domain.get("area_fraction", 0.0))
            gaps = self._domain_inter_axis_angles(domain)
            axis_name = ("actual two-axis angles"
                         if symmetry == "Square" else
                         "actual tri-axis angles")
            angle_values = (" / ".join("%.1f°" % float(value)
                                        for value in gaps) if gaps else "—")
            lines = [
                "Domain %d · %s · area %.1f%%" %
                (domain_id, symmetry, area),
                "a %s · orientation %.1f°" % (a_value, angle),
                "%s %s" % (axis_name, angle_values),
            ]
            painter.setFont(self._annotation_font(14, QFont.Weight.Bold))
            metrics = painter.fontMetrics()
            text_width = max(metrics.horizontalAdvance(line)
                             for line in lines)+20*unit
            line_height = max(15.0*unit, float(metrics.height()))
            text_height = line_height*len(lines)+12*unit
            x = cx+half+10*unit
            y = cy-text_height/2.0
            # Keep the information block close to its measured interior
            # marker and, whenever the measured bbox is large enough, clamp
            # it inside that domain rather than at a global image corner.
            bbox = domain.get("bbox") or []
            if len(bbox) == 4:
                left = max(8.0, float(bbox[0])-half)
                right = min(float(source.width())-8.0,
                            float(bbox[2])+half)
                top = max(8.0, float(bbox[1])-half)
                bottom = min(float(source.height())-8.0,
                             float(bbox[3])+half)
                if right-left >= text_width:
                    x = max(left, min(right-text_width, x))
                if bottom-top >= text_height:
                    y = max(top, min(bottom-text_height, y))
            x = max(8.0, min(source.width()-text_width-8.0, x))
            y = max(8.0, min(source.height()-text_height-8.0, y))
            box = QRectF(x, y, text_width, text_height)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(255, 255, 255, 218))
            painter.drawRoundedRect(box, 5*unit, 5*unit)
            painter.setPen(color)
            for line_index, line in enumerate(lines):
                line_box = QRectF(box.x()+8*unit,
                                  box.y()+6*unit+line_index*line_height,
                                  box.width()-16*unit, line_height)
                painter.drawText(line_box,
                                 Qt.AlignmentFlag.AlignVCenter |
                                 Qt.AlignmentFlag.AlignLeft,
                                 line)

    def _draw_single_twin_marker(self, painter, record):
        """Compatibility alias for older extensions."""
        self._draw_single_orientation_domains(painter, record)

    def _paint_annotations(self, painter, record, kind):
        source = QImage(record["source"])
        metrics = self._render_metrics(record)
        scale = max(1.0, min(source.width(), source.height()) / 900.0)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setFont(self._annotation_font(13, QFont.Weight.Bold))
        bilayer = record["result"].get("analysis_kind") == "bilayer"
        if kind == "original" and not bilayer:
            self._draw_domain_area_overlay(painter, record)
        if kind == "original" and bilayer:
            # Real-space unit cells and the representative TEM period belong
            # on the uploaded TEM.  Draw them before the FFT inset so no unit
            # annotation can spill across the inset panel.
            self._draw_moire_cells(painter, record, draw_cells=True)
        if kind in ("original", "ifft"):
            self._draw_fft_inset(painter, record, selected=(kind == "ifft"))
        lines = []
        if kind == "original":
            if bilayer:
                self._draw_tem_twist_marker(painter, record)
            else:
                self._draw_single_orientation_domains(painter, record)
            if bilayer and not metrics.get("fft_twist_reliable"):
                lines = ["TEM-derived twist unavailable · <2 reliable Moiré units",
                         "FFT-derived twist unavailable · peak sets too close",
                         "Only lattice constant a is reported"]
            elif bilayer and not metrics.get("tem_period_reliable"):
                lines = ["TEM field contains <2 reliable Moiré units",
                         "TEM-derived twist and moiré period unavailable · FFT-only"]
        elif kind == "reconstructed":
            self._draw_reconstruction_annotations(painter, record)
        elif kind == "ifft":
            self._draw_scale_bar(painter, record)
        if lines:
            metrics_font = painter.fontMetrics()
            margin = int(11 * scale)
            line_height = metrics_font.height() + int(4 * scale)
            box_width = max(metrics_font.horizontalAdvance(line)
                            for line in lines) + margin * 2
            box_height = line_height * len(lines) + margin
            x = (source.width() - box_width - int(28 * scale)
                 if kind == "reconstructed" else
                 (source.width() - box_width) / 2)
            y = int(28 * scale) if kind == "reconstructed" else int(250 * scale)
            painter.fillRect(QRectF(x, y, box_width, box_height),
                             QColor(0, 0, 0, 175))
            painter.setPen(QColor("white"))
            for index, line in enumerate(lines):
                painter.drawText(QPointF(x + margin,
                    y + margin + (index + 1) * line_height - 4), line)

    @staticmethod
    def _fmt(value, suffix=""):
        return "—" if value is None else "%.3f%s" % (float(value), suffix)

    def _base_image(self, record, kind):
        if kind == "reconstructed":
            return self._reconstructed_image(record)
        if kind == "ifft":
            return self._ifft_image(record)
        return QImage(record["source"]).convertToFormat(
            QImage.Format.Format_RGB32)

    def _render_analysis(self, record, kind, png_path=None, svg_path=None):
        base = self._base_image(record, kind)
        if base.isNull():
            raise ValueError("无法渲染图像：%s" % record["source"])
        raster = base.copy()
        self._normalize_annotation_dpi(raster)
        painter = QPainter(raster)
        self._paint_annotations(painter, record, kind)
        painter.end()
        if png_path:
            raster.save(str(png_path), "PNG")
        if svg_path:
            generator = QSvgGenerator()
            generator.setFileName(str(svg_path))
            generator.setSize(QSize(base.width(), base.height()))
            generator.setViewBox(QRect(0, 0, base.width(), base.height()))
            generator.setResolution(96)
            svg_painter = QPainter(generator)
            # Source and destination use the same aspect ratio and exact view
            # box. Raster TEM stays at native pixels; annotations remain vector.
            if kind == "reconstructed":
                # Paint every lattice and Moire-centre dot directly into the
                # SVG so each point remains an independently editable vector
                # object.  Preview and PNG rendering continue to use the same
                # geometry through the raster path above.
                self._reconstructed_image(record, painter=svg_painter)
            else:
                svg_painter.drawImage(
                    QRectF(0, 0, base.width(), base.height()),
                    base, QRectF(base.rect()))
            self._paint_annotations(svg_painter, record, kind)
            svg_painter.end()
            localize_svg(svg_path)
        return raster

    def _refresh_analysis_previews(self):
        record = self._current_record()
        if not record:
            return
        labels = (("original", self.tem_analysis_image),
                  ("reconstructed", self.reconstructed_analysis_image),
                  ("ifft", self.ifft_analysis_image))
        for kind, label in labels:
            image = self._render_analysis(record, kind)
            label.set_image(image)

    @staticmethod
    def _safe_stem(path):
        return "".join(character if character.isalnum() or character in "-_"
                       else "_" for character in Path(path).stem)

    def _export_analysis_record(self, record, folder):
        folder.mkdir(parents=True, exist_ok=True)
        exports = {}
        for kind, filename in (
                ("original", "original_TEM_annotated"),
                ("reconstructed", "reconstructed_bilayer"),
                ("ifft", "selected_spot_IFFT")):
            png = folder / (filename + ".png")
            svg = folder / (filename + ".svg")
            self._render_analysis(record, kind, png, svg)
            exports[kind] = {"png": str(png), "svg": str(svg)}
        metrics = self._render_metrics(record)
        report = {
            "sample": Path(record["source"]).name,
            "source": record["source"],
            "analysis_kind": record["result"].get(
                "analysis_kind", self.image_analysis_mode.currentData()),
            "overlap_mode": self._overlay_mode(),
            "overlap_percent": self._overlay_strength(),
            "metrics": metrics, "exports": exports,
            "automatic_analysis": record["result"],
        }
        (folder / "analysis.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2),
            encoding="utf-8")
        return report

    @staticmethod
    def _summary_row(report):
        """Flatten a report into the stable per-sample Bulk CSV schema."""
        metrics = report["metrics"]
        lattices = metrics.get("fft_lattices") or []
        twin = ((report.get("automatic_analysis", {}).get("lattice_fft") or
                 {}).get("twin") or {})
        twin_angles = twin.get("orientations_deg") or []
        all_single_angles = (report.get("automatic_analysis", {}).get(
            "single_layer_orientations_deg") or twin_angles)
        spatial = (report.get("automatic_analysis", {}).get(
            "orientation_domains") or {})
        row = {
            "sample": report["sample"],
            "analysis_kind": report["analysis_kind"],
            "symmetry": "; ".join(
                str(item.get("symmetry") or "") for item in lattices),
            "source_width_px": metrics.get("source_width_px"),
            "source_height_px": metrics.get("source_height_px"),
            "scale_bar_px": metrics.get("scale_bar_px"),
            "scale_bar_nm": metrics.get("scale_bar_nm"),
            "pixel_size_nm": metrics.get("pixel_size_nm"),
            "scale_source": metrics.get("scale_source"),
            "scale_warning": metrics.get("scale_warning"),
            # a is derived from the measured FFT; the TEM column is retained
            # in the file format but deliberately left blank.
            "TEM_a_nm": metrics.get("tem_a_nm"),
            "FFT_a_nm": metrics.get("fft_a_nm"),
            "Pore_a_nm": metrics.get("pore_a_nm"),
            "Twin_detected": bool(twin.get("valid")),
            "Twin_orientation_1_deg": (twin_angles[0]
                                       if len(twin_angles) > 0 else None),
            "Twin_orientation_2_deg": (twin_angles[1]
                                       if len(twin_angles) > 1 else None),
            "Twin_relative_orientation_deg": twin.get(
                "relative_orientation_deg"),
            "Single_layer_orientations_deg": ";".join(
                "%.6f" % float(value) for value in all_single_angles),
            "Spatial_domain_count": spatial.get("domain_count"),
            "Domain_boundary_count": spatial.get("boundary_count"),
            "TEM_twist_deg": metrics.get("tem_twist_deg"),
            "FFT_twist_deg": metrics.get("fft_twist_deg"),
            "TEM_period_nm": ("%.3f" % metrics["tem_period_nm"]
                              if metrics.get("tem_period_nm") is not None
                              else None),
            "FFT_period_nm": metrics.get("fft_period_nm"),
            "TEM_period_reliable": metrics.get("tem_period_reliable"),
            "overlap_mode": report["overlap_mode"],
            "overlap_percent": report["overlap_percent"],
        }
        row.update(AnalysisBulkMixin._spatial_domain_fields(
            report.get("automatic_analysis", {}),
            metrics.get("pixel_size_nm")))
        return row

    def save_image_analysis(self):
        if not self._analysis_records:
            QMessageBox.information(
                self, "没有结果",
                "请先完成批量图像分析。" if self._bulk_enabled() else
                "请先完成图像分析。")
            return
        root = (self._project_output_dir("analysis/moire_twist")
                if hasattr(self, "_project_output_dir") else None)
        if root is None:
            root = QFileDialog.getExistingDirectory(
                self, ("选择批量分析导出位置" if self._bulk_enabled() else
                       "选择当前分析结果的导出位置"),
                str(Path.home() / "Desktop"))
        if not root:
            return
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if not self._bulk_enabled():
            record = self._current_record() or self._analysis_records[0]
            sample = self._safe_stem(record["source"])
            output = Path(root) / (sample + "_TEM_analysis_" + stamp)
            counter = 1
            while output.exists():
                output = Path(root) / (
                    sample + "_TEM_analysis_%s_%d" % (stamp, counter))
                counter += 1
            self._export_analysis_record(record, output)
            QMessageBox.information(
                self, "当前分析结果已导出",
                "已按当前Overlay比例手动导出PNG、SVG和JSON。\n\n%s\n\n"
                "分析阶段没有自动写出任何SVG。" % output)
            return

        output = Path(root) / ("TEM_bulk_analysis_" + stamp)
        counter = 1
        while output.exists():
            output = Path(root) / ("TEM_bulk_analysis_%s_%d" % (stamp, counter))
            counter += 1
        output.mkdir(parents=True)
        rows, reports = [], []
        current = self.analysis_file_selector.currentIndex()
        for index, record in enumerate(self._analysis_records):
            self.analysis_file_selector.setCurrentIndex(index)
            QApplication.processEvents()
            sample = self._safe_stem(record["source"])
            folder = output / sample
            suffix = 1
            while folder.exists():
                folder = output / (sample + "_%d" % suffix); suffix += 1
            folder.mkdir()
            report = self._export_analysis_record(record, folder)
            metrics = report["metrics"]
            twin = ((record["result"].get("lattice_fft") or {}).get(
                "twin") or {})
            twin_angles = twin.get("orientations_deg") or []
            all_single_angles = (record["result"].get(
                "single_layer_orientations_deg") or twin_angles)
            spatial = record["result"].get("orientation_domains") or {}
            reports.append(report)
            row = {
                "sample": report["sample"],
                "analysis_kind": report["analysis_kind"],
                "source_width_px": metrics["source_width_px"],
                "source_height_px": metrics["source_height_px"],
                "scale_bar_px": metrics["scale_bar_px"],
                "scale_bar_nm": metrics["scale_bar_nm"],
                "TEM_a_nm": metrics["tem_a_nm"],
                "FFT_a_nm": metrics["fft_a_nm"],
                "Pore_a_nm": metrics.get("pore_a_nm"),
                "Twin_detected": bool(twin.get("valid")),
                "Twin_orientation_1_deg": (twin_angles[0]
                                           if len(twin_angles) > 0 else None),
                "Twin_orientation_2_deg": (twin_angles[1]
                                           if len(twin_angles) > 1 else None),
                "Twin_relative_orientation_deg": twin.get(
                    "relative_orientation_deg"),
                "Single_layer_orientations_deg": ";".join(
                    "%.6f" % float(value) for value in all_single_angles),
                "Spatial_domain_count": spatial.get("domain_count"),
                "Domain_boundary_count": spatial.get("boundary_count"),
                "TEM_twist_deg": metrics["tem_twist_deg"],
                "FFT_twist_deg": metrics["fft_twist_deg"],
                "TEM_period_nm": ("%.3f" % metrics["tem_period_nm"]
                                  if metrics.get("tem_period_nm") is not None
                                  else None),
                "FFT_period_nm": metrics["fft_period_nm"],
                "final_twist_deg": metrics["final_twist_deg"],
                "overlap_mode": report["overlap_mode"],
                "overlap_percent": report["overlap_percent"],
            }
            row.update(self._spatial_domain_fields(
                record["result"], metrics.get("pixel_size_nm")))
            rows.append(row)
        if current >= 0:
            self.analysis_file_selector.setCurrentIndex(current)
        fields = list(rows[0].keys())
        with (output / "bulk_analysis_summary.csv").open(
                "w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
        localize_csv(output / "bulk_analysis_summary.csv")
        (output / "bulk_analysis_summary.json").write_text(
            json.dumps(reports, ensure_ascii=False, indent=2), encoding="utf-8")
        QMessageBox.information(
            self, "批量分析已导出",
            "已导出%d个样品。\n\n%s\n\nCSV包含每个样品TEM与FFT各自的"
            "a、twist和period；所有PNG保持原始像素，SVG保持相同比例并保留"
            "可编辑矢量标注。" % (len(rows), output))

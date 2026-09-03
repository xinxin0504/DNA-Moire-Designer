"""Curved Design dialog backed by the bundled DNAxiS routing engine."""

import math
import os

from ..model.io.curved import (CURVED_SCAFFOLD_MAX_BASES, PRESETS,
                               build_rings, curved_indel_plan, curved_output_name,
                               estimated_scaffold_bases, safe_name)
import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), [
    'QDir', 'QRectF', 'QSettings', 'Qt'])
util.qtWrapImport('QtGui', globals(), [
    'QBrush', 'QColor', 'QPainter', 'QPen'])
util.qtWrapImport('QtWidgets', globals(), [
    'QComboBox', 'QDialog', 'QDialogButtonBox', 'QDoubleSpinBox',
    'QFileDialog', 'QFormLayout', 'QHBoxLayout', 'QLabel', 'QLineEdit',
    'QMessageBox', 'QPushButton', 'QScrollArea', 'QSpinBox', 'QVBoxLayout',
    'QWidget'])


class CurvedPreview(QWidget):
    def __init__(self, parent=None):
        super(CurvedPreview, self).__init__(parent)
        self._rings = []
        self.setMinimumSize(520, 300)

    def setRings(self, rings):
        self._rings = list(rings)
        self.update()

    def paintEvent(self, unused_event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor(248, 249, 251))
        if not self._rings:
            return
        max_radius = max(float(row["radius_nm"]) for row in self._rings)
        max_height = max(float(row["height_nm"]) for row in self._rings) or 1.0
        bounds = self.rect().adjusted(35, 25, -35, -25)
        scale = min(bounds.width() / (2.0 * max_radius),
                    bounds.height() / max_height)
        center_x = bounds.center().x()
        base_y = bounds.bottom()
        for ring in self._rings:
            radius = float(ring["radius_nm"]) * scale
            y_value = base_y - float(ring["height_nm"]) * scale
            ellipse_height = max(3.0, radius * 0.32)
            color = (QColor(49, 100, 169) if int(ring["layer"]) == 0
                     else QColor(123, 151, 190))
            painter.setPen(QPen(color, 1.2))
            painter.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            painter.drawEllipse(QRectF(
                center_x - radius, y_value - ellipse_height / 2.0,
                radius * 2.0, ellipse_height))


class CurvedDesignDialog(QDialog):
    def __init__(self, parent=None):
        super(CurvedDesignDialog, self).__init__(parent)
        self.setWindowTitle("Curved Design")
        self.resize(760, 760)
        self._projectEditedByUser = False
        self._lastRings = []

        outer = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        page = QWidget(scroll)
        layout = QVBoxLayout(page)
        scroll.setWidget(page)
        outer.addWidget(scroll, 1)

        intro = QLabel(
            "使用DNAxiS生成轴对称曲面DNA折纸。自动生成后会打开可编辑的"
            "cadnano设计；尺寸表示DNA双螺旋的外轮廓，基础层沿子午线"
            "按约2.8 nm取样，多层采用2.8 nm径向叠层。", self)
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.shapeBox = QComboBox(self)
        for key, label, description in PRESETS:
            self.shapeBox.addItem("%s — %s" % (label, description), key)
        form.addRow("预设形状：", self.shapeBox)

        self.latticeBox = QComboBox(self)
        self.latticeBox.addItem("Square（连续双链优先至少16 bp）",
                                "square")
        self.latticeBox.addItem("Honeycomb（连续双链优先至少14 bp）",
                                "honeycomb")
        self.latticeBox.setToolTip(
            "crossover使用所选点阵的合法相位、密度和间距规则。"
            "Honeycomb采用A/镜像A交错结构并支持1–3层。")
        form.addRow("点阵通道：", self.latticeBox)

        densityRow = QHBoxLayout()
        self.scaffoldDensityBox = QComboBox(self)
        self.scaffoldDensityCustomSpin = QSpinBox(self)
        self.scaffoldDensityCustomSpin.setRange(1, 32)
        self.scaffoldDensityCustomSpin.setValue(3)
        self.scaffoldDensityCustomSpin.setPrefix("× ")
        self.scaffoldDensityCustomSpin.setToolTip(
            "自定义scaffold crossover间距为所选点阵原生周期的整数倍。")
        self.scaffoldDensityValueLabel = QLabel(self)
        self.scaffoldDensityValueLabel.setWordWrap(True)
        densityRow.addWidget(self.scaffoldDensityBox, 1)
        densityRow.addWidget(self.scaffoldDensityCustomSpin)
        densityRow.addWidget(self.scaffoldDensityValueLabel)
        form.addRow("Scaffold crossover密度：", densityRow)

        self.heightSpin = QDoubleSpinBox(self)
        self.heightSpin.setRange(5.2, 200.0)
        self.heightSpin.setDecimals(1)
        self.heightSpin.setSingleStep(2.8)
        self.heightSpin.setValue(31.2)
        self.heightSpin.setSuffix(" nm")
        form.addRow("目标外轮廓总高度：", self.heightSpin)

        self.maximumDiameterSpin = QDoubleSpinBox(self)
        self.maximumDiameterSpin.setRange(10.0, 120.0)
        self.maximumDiameterSpin.setDecimals(1)
        self.maximumDiameterSpin.setValue(24.0)
        self.maximumDiameterSpin.setSuffix(" nm")
        form.addRow("目标外轮廓最大直径：", self.maximumDiameterSpin)

        self.minimumDiameterSpin = QDoubleSpinBox(self)
        self.minimumDiameterSpin.setRange(10.0, 120.0)
        self.minimumDiameterSpin.setDecimals(1)
        self.minimumDiameterSpin.setValue(12.2)
        self.minimumDiameterSpin.setSuffix(" nm")
        form.addRow("目标外轮廓最小直径：", self.minimumDiameterSpin)

        self.layersSpin = QSpinBox(self)
        self.layersSpin.setRange(1, 3)
        self.layersSpin.setValue(1)
        form.addRow("壁层数：", self.layersSpin)

        self.nameEdit = QLineEdit("curved", self)
        form.addRow("设计名称：", self.nameEdit)

        projectRow = QHBoxLayout()
        self.projectEdit = QLineEdit(self)
        self.projectButton = QPushButton("选择…", self)
        projectRow.addWidget(self.projectEdit, 1)
        projectRow.addWidget(self.projectButton)
        form.addRow("项目文件夹：", projectRow)
        layout.addLayout(form)

        self.capacityLabel = QLabel(self)
        self.capacityLabel.setWordWrap(True)
        self.capacityLabel.setMinimumHeight(68)
        self.capacityLabel.setMargin(8)
        layout.addWidget(self.capacityLabel)

        previewTitle = QLabel("可旋转体环形结构预览：", self)
        layout.addWidget(previewTitle)
        self.preview = CurvedPreview(self)
        layout.addWidget(self.preview)

        output = QLabel(
            "自动保存：主文件夹中的可编辑JSON；input子文件夹中的"
            "STL、design settings、modules CSV和preview PNG。所有文件名"
            "均包含形状、点阵、层数和尺寸信息。"
            "初始阶段不导出序列、PDB/mmCIF或TOP/DAT。", self)
        output.setWordWrap(True)
        layout.addWidget(output)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel, self)
        self.buttons.button(
            QDialogButtonBox.StandardButton.Ok).setText("运行Curved Design")
        outer.addWidget(self.buttons)

        for widget in (self.shapeBox, self.latticeBox, self.heightSpin,
                       self.maximumDiameterSpin, self.minimumDiameterSpin,
                       self.layersSpin):
            signal = (widget.currentIndexChanged if
                      isinstance(widget, QComboBox) else widget.valueChanged)
            signal.connect(self._refresh)
        self.nameEdit.textChanged.connect(self._nameChanged)
        self.projectEdit.textEdited.connect(self._projectWasEdited)
        self.projectButton.clicked.connect(self._chooseProject)
        self.buttons.accepted.connect(self._acceptIfValid)
        self.buttons.rejected.connect(self.reject)
        self.latticeBox.currentIndexChanged.connect(
            self._updateDiameterLimits)
        self.latticeBox.currentIndexChanged.connect(
            self._updateScaffoldDensityChoices)
        self.scaffoldDensityBox.currentIndexChanged.connect(
            self._updateScaffoldDensityState)
        self.scaffoldDensityCustomSpin.valueChanged.connect(
            self._updateScaffoldDensityState)
        self._updateScaffoldDensityChoices(refresh=False)
        self._updateDiameterLimits()
        self._refresh()

    def _defaultProject(self):
        shape = str(self.shapeBox.currentData())
        output = curved_output_name(
            self.nameEdit.text(), shape,
            self.latticeBox.currentData(), self.layersSpin.value(),
            self.heightSpin.value(), self.maximumDiameterSpin.value(),
            self.minimumDiameterSpin.value())
        desktop = os.path.join(QDir.homePath(), "Desktop")
        return os.path.join(desktop, output)

    def _nameChanged(self, unused_text):
        if not self._projectEditedByUser:
            self.projectEdit.setText(self._defaultProject())

    def _projectWasEdited(self, unused_text):
        self._projectEditedByUser = True

    def _chooseProject(self):
        selected = QFileDialog.getExistingDirectory(
            self, "选择Curved Design项目的父文件夹",
            os.path.dirname(self.projectEdit.text()) or QDir.homePath())
        if isinstance(selected, (tuple, list)):
            selected = selected[0]
        if selected:
            self._projectEditedByUser = True
            self.projectEdit.setText(os.path.join(
                str(selected), curved_output_name(
                    self.nameEdit.text(), self.shapeBox.currentData(),
                    self.latticeBox.currentData(), self.layersSpin.value(),
                    self.heightSpin.value(),
                    self.maximumDiameterSpin.value(),
                    self.minimumDiameterSpin.value())))

    def _updateDiameterLimits(self, unused_index=None):
        minimum = (10.9 if str(self.latticeBox.currentData()) ==
                   "honeycomb" else 12.2)
        for spin in (self.maximumDiameterSpin, self.minimumDiameterSpin):
            blocked = spin.blockSignals(True)
            spin.setMinimum(minimum)
            if spin.value() < minimum:
                spin.setValue(minimum)
            spin.blockSignals(blocked)
        if unused_index is not None:
            self._refresh()

    def _updateScaffoldDensityChoices(self, unused_index=None,
                                      refresh=True):
        selected = (self.scaffoldDensityBox.currentData()
                    if self.scaffoldDensityBox.count() else 1)
        native = (21 if str(self.latticeBox.currentData()) == "honeycomb"
                  else 32)
        blocked = self.scaffoldDensityBox.blockSignals(True)
        self.scaffoldDensityBox.clear()
        self.scaffoldDensityBox.addItem(
            "1/%d bp（最大合法密度）" % native, 1)
        self.scaffoldDensityBox.addItem(
            "1/%d bp（2个点阵周期）" % (2 * native), 2)
        self.scaffoldDensityBox.addItem("自定义周期倍数", -1)
        self.scaffoldDensityBox.addItem(
            "最低路由密度（仅保留必要连接）", 0)
        target = self.scaffoldDensityBox.findData(selected)
        self.scaffoldDensityBox.setCurrentIndex(max(0, target))
        self.scaffoldDensityBox.blockSignals(blocked)
        self._updateScaffoldDensityState(refresh=refresh)

    def _scaffoldDensitySettings(self):
        value = int(self.scaffoldDensityBox.currentData())
        if value == 0:
            return "minimum", 0
        if value < 0:
            value = int(self.scaffoldDensityCustomSpin.value())
        return "periodic", max(1, value)

    def _updateScaffoldDensityState(self, unused_value=None, refresh=True):
        custom = int(self.scaffoldDensityBox.currentData()) == -1
        self.scaffoldDensityCustomSpin.setEnabled(custom)
        mode, multiple = self._scaffoldDensitySettings()
        native = (21 if str(self.latticeBox.currentData()) == "honeycomb"
                  else 32)
        if mode == "minimum":
            text = "单一scaffold所需最低密度"
        else:
            text = "目标1/%d bp" % (native * multiple)
        self.scaffoldDensityValueLabel.setText(text)
        self.scaffoldDensityBox.setToolTip(
            "所有选项仍强制合法点阵相位、单一scaffold和唯一nick。"
            "若所选低密度无法形成合法整体路由，设计会明确报错而不会"
            "静默违反密度或拓扑规则。")
        if refresh:
            self._refresh()

    def _refresh(self, unused_value=None):
        try:
            rings = build_rings(
                self.shapeBox.currentData(), self.heightSpin.value(),
                self.maximumDiameterSpin.value(),
                self.minimumDiameterSpin.value(), self.layersSpin.value(),
                lattice=str(self.latticeBox.currentData()))
            lattice = str(self.latticeBox.currentData())
            rings, indelPlan = curved_indel_plan(rings, lattice)
            required = estimated_scaffold_bases(rings)
            message = ("设计长度约%d nt，共%d个DNA圆环。" %
                       (required, len(rings)))
            actualHeight = (max(float(ring["height_nm"]) for ring in rings) -
                            min(float(ring["height_nm"]) for ring in rings) +
                            2.0)
            actualDiameters = [2.0 * (float(ring["radius_nm"]) + 1.0)
                               for ring in rings]
            message += (
                " 目标外高%.2f nm、目标外径%.2f–%.2f nm；"
                "生成后实际外高约%.2f nm、实际全部helix外径%.2f–%.2f "
                "nm。" %
                (float(self.heightSpin.value()),
                 float(self.minimumDiameterSpin.value()),
                 float(self.maximumDiameterSpin.value()),
                 actualHeight, min(actualDiameters),
                 max(actualDiameters)))
            message += (" 子午线圆环间距约%.2f nm；采用%d-bp domain，"
                        "预计最大insertion/domain=%d、"
                        "deletion/domain=%d（硬上限±3）。" %
                        (float(rings[0].get("meridian_spacing_nm", 2.8)),
                         int(indelPlan["domain_size_bp"]),
                         int(indelPlan["maximum_insertion_per_domain"]),
                         int(indelPlan["maximum_deletion_per_domain"])))
            if not indelPlan.get("domain_limit_feasible", False):
                message += (" 当前几何无法满足±3/domain硬上限；请减小曲率、"
                            "增加轴向跨度或更改截面。")
            if self.maximumDiameterSpin.value() > 60.0:
                message += (" 最大直径超过60 nm；允许生成，但属于大型柔性"
                            "结构，请重点检查scaffold容量和crossover统计。")
            message += (" 所有helix先采用相同的%d bp整周期母结构；"
                        "公共parent参考长度为%d bp。该parent始终按所有"
                        "helix的最大绝对indel负载最小化选择（例如"
                        "−4/0/+4优于0/+4/+8），并非仅在触及±3/domain"
                        "上限时才调整。" %
                        (int(indelPlan["common_nominal_bp"]),
                         int(indelPlan["optimized_parent_reference_bp"])))
            densityMode, densityMultiple = self._scaffoldDensitySettings()
            if densityMode == "minimum":
                message += (" Scaffold crossover采用最低路由密度，"
                            "只保留形成单一scaffold和唯一nick所需连接。")
            else:
                nativePeriod = 21 if lattice == "honeycomb" else 32
                message += (" Scaffold crossover目标密度为1/%d bp。" %
                            (nativePeriod * densityMultiple))
            diameterShift = float(
                indelPlan.get("diameter_adjustment_nm", 0.0))
            if abs(diameterShift) > 0.005:
                message += (" 整体外径统一调整%+.2f nm。" % diameterShift)
            if indelPlan.get("cylinder_period_snap"):
                message += " 单层圆柱没有相对半径变化，因此不使用indel。"
            if required <= CURVED_SCAFFOLD_MAX_BASES:
                message += (" scaffold容量：%d/%d nt；不再受内置序列长度"
                            "限制，设计完成后可按实际需要导入或分配序列。" %
                            (required, CURVED_SCAFFOLD_MAX_BASES))
                style = ("QLabel { color:#1f6130; background:#eff8f1; "
                         "border:1px solid #9ac8a5; border-radius:5px; }")
            else:
                message += (" scaffold容量：%d/%d nt；当前设计超过Curved "
                            "Design上限，请减小尺寸或层数。" %
                            (required, CURVED_SCAFFOLD_MAX_BASES))
                style = ("QLabel { color:#8f1d1d; background:#fff1f0; "
                         "border:1px solid #e5a3a0; border-radius:5px; }")
            self.capacityLabel.setText(message)
            self.capacityLabel.setStyleSheet(style)
            self.capacityLabel.setProperty("curvedRequiredBases", required)
            self._lastRings = rings
            self.preview.setRings(rings)
        except Exception as error:
            self.capacityLabel.setText("参数无效：%s" % error)
            self.capacityLabel.setProperty("curvedRequiredBases", 0)
            self._lastRings = []
            self.preview.setRings([])
        if not self._projectEditedByUser:
            self.projectEdit.setText(self._defaultProject())

    def _acceptIfValid(self):
        if not self._lastRings:
            QMessageBox.warning(self, "Curved Design", "请先修正形状参数。")
            return
        required = int(self.capacityLabel.property(
            "curvedRequiredBases") or 0)
        if required > CURVED_SCAFFOLD_MAX_BASES:
            QMessageBox.warning(
                self, "Curved Design",
                "当前设计需要%d nt scaffold，超过%d nt上限。\n"
                "请减小尺寸或层数。" %
                (required, CURVED_SCAFFOLD_MAX_BASES))
            return
        if not safe_name(self.nameEdit.text()):
            QMessageBox.warning(self, "Curved Design", "请输入设计名称。")
            return
        if not self.projectEdit.text().strip():
            QMessageBox.warning(self, "Curved Design", "请选择项目文件夹。")
            return
        self.accept()

    def spec(self):
        densityMode, densityMultiple = self._scaffoldDensitySettings()
        return {
            "shape": str(self.shapeBox.currentData()),
            "lattice": str(self.latticeBox.currentData()),
            "height_nm": float(self.heightSpin.value()),
            "maximum_diameter_nm": float(
                self.maximumDiameterSpin.value()),
            "minimum_diameter_nm": float(
                self.minimumDiameterSpin.value()),
            "layers": int(self.layersSpin.value()),
            "scaffold_crossover_density_mode": densityMode,
            "scaffold_crossover_density_multiple": densityMultiple,
            "name": safe_name(self.nameEdit.text()),
            "project_root": os.path.abspath(self.projectEdit.text().strip())}

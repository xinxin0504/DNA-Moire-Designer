"""Modeless beginner workflow for image-to-origami design."""

from PyQt6 import QtCore, QtGui, QtWidgets

from cadnano2.model.guideddesign import (base_count_for_width,
                                         boolean_runs,
                                         cross_section_coords,
                                         estimate_scaffold_length,
                                         profile_count_for_height)


class PatternPreview(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(210)
        self._image = QtGui.QImage()
        self._spec = None

    def setImage(self, image):
        self._image = image
        self.update()

    def setSpec(self, spec):
        self._spec = spec
        self.update()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QtGui.QColor('#fafafa'))
        left = QtCore.QRectF(12, 12, self.width() * .54, self.height() - 24)
        painter.setPen(QtGui.QPen(QtGui.QColor('#b7b7b7'), 1))
        painter.drawRect(left)
        if not self._image.isNull():
            target = QtCore.QRectF(left)
            size = self._image.size()
            size.scale(int(left.width()), int(left.height()),
                       QtCore.Qt.AspectRatioMode.KeepAspectRatio)
            target.setSize(QtCore.QSizeF(size))
            target.moveCenter(left.center())
            painter.drawImage(target, self._image)
        if not self._spec:
            return
        coords = self._spec.get('coords', [])
        if not coords:
            return
        x0 = self.width() * .62
        area_w = self.width() * .35
        area_h = self.height() - 24
        max_row = max(row for row, col, profile, layer in coords) + 1
        max_col = max(col for row, col, profile, layer in coords) + 1
        diameter = max(5.0, min(20.0, area_w / max_col * .7,
                                area_h / max_row * .7))
        x_step = area_w / max(1, max_col)
        y_step = area_h / max(1, max_row)
        for row, col, profile, layer in coords:
            x = x0 + (col + .5) * x_step
            y = 12 + (row + .5) * y_step
            if self._spec['lattice'] == 'honeycomb' and ((row ^ col) & 1):
                y += min(y_step * .28, diameter * .35)
            color = QtGui.QColor('#dcecff' if layer % 2 == 0 else '#ffe1ee')
            painter.setBrush(color)
            painter.setPen(QtGui.QPen(QtGui.QColor('#56616b'), 1))
            painter.drawEllipse(QtCore.QPointF(x, y), diameter / 2,
                                diameter / 2)
        painter.setPen(QtGui.QColor('#333333'))
        painter.drawText(QtCore.QRectF(x0, 2, area_w, 20),
                         QtCore.Qt.AlignmentFlag.AlignCenter,
                         '截面层排列')


class GuidedDesignWizard(QtWidgets.QWizard):
    def __init__(self, controller):
        super().__init__(controller.window())
        self.controller = controller
        self.sourceImage = self._exampleImage('rounded')
        self.profileRuns = []
        self._lastSummary = {}
        self.setWindowTitle('智能引导设计')
        self.setWindowModality(QtCore.Qt.WindowModality.NonModal)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.resize(760, 650)
        self._addPatternPage()
        self._addGeometryPage()
        self._addPreviewPage()
        self._addBuildPage()
        self._addSequencePage()
        self.currentIdChanged.connect(self._pageChanged)

    def _page(self, title, subtitle):
        page = QtWidgets.QWizardPage()
        page.setTitle(title)
        page.setSubTitle(subtitle)
        page.setLayout(QtWidgets.QVBoxLayout())
        self.addPage(page)
        return page

    def _addPatternPage(self):
        page = self._page('1. 选择目标图案',
                          '默认情况下，螺旋方向与图案水平方向平行。')
        row = QtWidgets.QHBoxLayout()
        self.exampleCombo = QtWidgets.QComboBox()
        self.exampleCombo.addItem('圆角矩形', 'rounded')
        self.exampleCombo.addItem('圆环', 'ring')
        self.exampleCombo.addItem('三角形', 'triangle')
        self.exampleCombo.currentIndexChanged.connect(self._chooseExample)
        upload = QtWidgets.QPushButton('上传图案图片…')
        upload.clicked.connect(self._uploadImage)
        row.addWidget(QtWidgets.QLabel('简单示例：'))
        row.addWidget(self.exampleCombo, 1)
        row.addWidget(upload)
        page.layout().addLayout(row)
        self.patternPreview = PatternPreview()
        self.patternPreview.setImage(self.sourceImage)
        page.layout().addWidget(self.patternPreview, 1)
        note = QtWidgets.QLabel(
            '图片中较暗或不透明的区域会生成 DNA；白色或透明区域保持为空。'
            '增加层数时只改变截面中的螺旋排列，不会镜像目标图案本身。')
        note.setWordWrap(True)
        page.layout().addWidget(note)

    def _addGeometryPage(self):
        page = self._page('2. 设置尺寸和点阵',
                          '尺寸换算采用每个碱基 0.34 nm、螺旋中心间距 2.8 nm。')
        form = QtWidgets.QFormLayout()
        self.widthSpin = QtWidgets.QDoubleSpinBox()
        self.widthSpin.setRange(10, 1000)
        self.widthSpin.setValue(80)
        self.widthSpin.setSuffix(' nm')
        self.heightSpin = QtWidgets.QDoubleSpinBox()
        self.heightSpin.setRange(2.8, 1000)
        self.heightSpin.setValue(28)
        self.heightSpin.setSuffix(' nm')
        self.orientationCombo = QtWidgets.QComboBox()
        self.orientationCombo.addItems(['水平方向螺旋', '竖直方向螺旋'])
        self.latticeCombo = QtWidgets.QComboBox()
        self.latticeCombo.addItem('蜂窝点阵', 'honeycomb')
        self.latticeCombo.addItem('方形点阵', 'square')
        self.layersSpin = QtWidgets.QSpinBox()
        self.layersSpin.setRange(1, 40)
        self.layersSpin.setValue(1)
        self.directionCombo = QtWidgets.QComboBox()
        self.directionCombo.addItem('Z 方向——A／镜像A交替', 'z')
        self.directionCombo.addItem('右上点阵方向', 'up-right')
        self.directionCombo.addItem('右下点阵方向', 'down-right')
        self.latticeCombo.currentIndexChanged.connect(self._latticeChanged)
        form.addRow('沿螺旋方向的图案最大尺寸：', self.widthSpin)
        form.addRow('垂直于螺旋方向的图案最大尺寸：', self.heightSpin)
        form.addRow('螺旋方向：', self.orientationCombo)
        form.addRow('点阵类型：', self.latticeCombo)
        form.addRow('层数：', self.layersSpin)
        form.addRow('蜂窝层方向：', self.directionCombo)
        page.layout().addLayout(form)
        explanation = QtWidgets.QLabel(
            '方形点阵会沿截面方向复制同一个平面图案。蜂窝点阵的 Z 方向层从上'
            '到下按照 A、镜像A、A、镜像A 交替；只有螺旋位置交替，图案不翻转。')
        explanation.setWordWrap(True)
        page.layout().addWidget(explanation)
        page.layout().addStretch(1)

    def _addPreviewPage(self):
        page = self._page('3. 预览螺旋布局',
                          '生成前请检查采样后的图案轮廓和截面层排列。')
        self.layoutPreview = PatternPreview()
        page.layout().addWidget(self.layoutPreview, 1)
        self.previewText = QtWidgets.QLabel()
        self.previewText.setWordWrap(True)
        page.layout().addWidget(self.previewText)

    def _addBuildPage(self):
        page = self._page('4. 生成并检查 DNA 设计',
                          '引导窗口会保持打开，各步骤之间可以在 cadnano 主界面手动修改。')
        self.generateButton = QtWidgets.QPushButton('生成螺旋布局')
        self.scaffoldButton = QtWidgets.QPushButton(
            '生成闭环骨架链并保留一个切口（首个 crossover 位于 base 50 之后）')
        self.stapleButton = QtWidgets.QPushButton(
            '添加短链、短链 crossover 并自动优化断点')
        self.scaffoldButton.setEnabled(False)
        self.stapleButton.setEnabled(False)
        self.generateButton.clicked.connect(self._generateLayout)
        self.scaffoldButton.clicked.connect(self._buildScaffold)
        self.stapleButton.clicked.connect(self._buildStaples)
        page.layout().addWidget(self.generateButton)
        page.layout().addWidget(self.scaffoldButton)
        page.layout().addWidget(self.stapleButton)
        self.buildStatus = QtWidgets.QLabel(
            '请从空白文档开始。每一步生成的结果都可以在普通二维视图中继续编辑。')
        self.buildStatus.setWordWrap(True)
        self.buildStatus.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        page.layout().addWidget(self.buildStatus)
        page.layout().addStretch(1)

    def _addSequencePage(self):
        page = self._page('5. 添加序列、保存和导出',
                          '只显示长度不短于当前设计的骨架链序列。')
        row = QtWidgets.QHBoxLayout()
        self.sequenceCombo = QtWidgets.QComboBox()
        refresh = QtWidgets.QPushButton('刷新可用序列')
        apply_sequence = QtWidgets.QPushButton('应用所选骨架链序列')
        refresh.clicked.connect(self._refreshSequences)
        apply_sequence.clicked.connect(self._applySequence)
        row.addWidget(self.sequenceCombo, 1)
        row.addWidget(refresh)
        page.layout().addLayout(row)
        page.layout().addWidget(apply_sequence)
        self.sequenceStatus = QtWidgets.QLabel(
            '如果所选序列长于设计，只使用设计所需的部分；多出的碱基会被忽略，'
            '不会在切口处增加 insertion 或单链环。')
        self.sequenceStatus.setWordWrap(True)
        page.layout().addWidget(self.sequenceStatus)
        exports = QtWidgets.QGridLayout()
        actions = [
            ('保存（不保存序列）', self.controller.actionSaveSlot),
            ('另存为（包含序列）',
             self.controller.actionSaveWithSequencesSlot),
            ('导出 Illustrator 矢量图', self.controller.actionIllustratorSlot),
            ('导出序列 XLSX', self.controller.actionExportStaplesSlot),
            ('导出结构（oxDNA 可选）',
             self.controller.actionGuidedExportPdbSlot),
        ]
        for index, (text, slot) in enumerate(actions):
            button = QtWidgets.QPushButton(text)
            button.clicked.connect(slot)
            exports.addWidget(button, index // 2, index % 2)
        page.layout().addLayout(exports)
        page.layout().addStretch(1)

    def spec(self):
        lattice = self.latticeCombo.currentData()
        image = self.sourceImage.transformed(QtGui.QTransform().rotate(90)) \
            if self.orientationCombo.currentIndex() else self.sourceImage
        bases = base_count_for_width(self.widthSpin.value())
        profiles = profile_count_for_height(self.heightSpin.value())
        self.profileRuns = self._sampleImage(image, profiles, bases)
        direction = self.directionCombo.currentData()
        coords = cross_section_coords(profiles, self.layersSpin.value(),
                                      lattice, direction)
        return {'lattice': lattice, 'base_count': bases,
                'profile_count': profiles, 'layers': self.layersSpin.value(),
                'direction': direction, 'coords': coords,
                'profile_runs': self.profileRuns,
                'estimated_length': estimate_scaffold_length(
                    self.profileRuns, self.layersSpin.value())}

    def _pageChanged(self, page_id):
        if page_id == 2:
            spec = self.spec()
            self.layoutPreview.setImage(self.sourceImage)
            self.layoutPreview.setSpec(spec)
            occupied = sum(1 for runs in spec['profile_runs'] if runs)
            self.previewText.setText(
                '{lattice_name} · {layers} 层 · 每层采样 '
                '{profile_count} 条轮廓螺旋（其中 {occupied} 条被图案占用）· '
                '共 {helix_count} 条点阵螺旋 · 预计需要约 '
                '{estimated_length} nt 骨架链。首个自动骨架链 crossover '
                '将位于 cadnano base 51 或之后。'
                .format(occupied=occupied,
                        helix_count=len(spec['coords']),
                        lattice_name=('蜂窝点阵' if spec['lattice'] ==
                                      'honeycomb' else '方形点阵'), **spec))
        elif page_id == 4:
            self._refreshSequences()

    def _sampleImage(self, image, profiles, bases):
        if image.isNull():
            return [[] for unused in range(profiles)]
        rows = []
        for profile in range(profiles):
            y = min(image.height() - 1,
                    int((profile + .5) * image.height() / profiles))
            values = []
            for x in range(image.width()):
                pixel = image.pixelColor(x, y)
                darkness = (pixel.red() + pixel.green() + pixel.blue()) / 3
                values.append(pixel.alpha() > 20 and darkness < 235)
            rows.append(boolean_runs(values, bases))
        return rows

    def _exampleImage(self, name):
        image = QtGui.QImage(420, 220, QtGui.QImage.Format.Format_ARGB32)
        image.fill(QtGui.QColor('white'))
        painter = QtGui.QPainter(image)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor('#2f3439'))
        if name == 'ring':
            painter.drawEllipse(45, 25, 330, 170)
            painter.setBrush(QtGui.QColor('white'))
            painter.drawEllipse(115, 75, 190, 70)
        elif name == 'triangle':
            painter.drawPolygon(QtGui.QPolygonF([
                QtCore.QPointF(210, 18), QtCore.QPointF(390, 200),
                QtCore.QPointF(30, 200)]))
        else:
            painter.drawRoundedRect(35, 28, 350, 164, 42, 42)
        painter.end()
        return image

    def _chooseExample(self, unused_index):
        self.sourceImage = self._exampleImage(self.exampleCombo.currentData())
        self.patternPreview.setImage(self.sourceImage)

    def _uploadImage(self):
        filename, unused_filter = QtWidgets.QFileDialog.getOpenFileName(
            self, '选择目标图案', '',
            '图案图片 (*.png *.jpg *.jpeg *.bmp *.tif *.tiff)')
        if filename:
            image = QtGui.QImage(filename)
            if image.isNull():
                QtWidgets.QMessageBox.warning(self, '智能引导设计',
                                              '无法读取所选图片。')
                return
            self.sourceImage = image
            self.patternPreview.setImage(image)

    def _latticeChanged(self, unused_index):
        self.directionCombo.setEnabled(
                                self.latticeCombo.currentData() == 'honeycomb')

    def _generateLayout(self):
        result = self.controller.guidedGenerateLayout(self.spec())
        self._showResult(result)
        if result.get('ok'):
            self.scaffoldButton.setEnabled(True)

    def _buildScaffold(self):
        result = self.controller.guidedBuildScaffold()
        self._showResult(result)
        if result.get('ok'):
            self.stapleButton.setEnabled(True)
            self._refreshSequences()

    def _buildStaples(self):
        self._showResult(self.controller.guidedBuildStaples())

    def _showResult(self, result):
        self._lastSummary = result
        self.buildStatus.setText(result.get('message', ''))
        if not result.get('ok'):
            QtWidgets.QMessageBox.warning(self, '智能引导设计',
                                          result.get('message', '无法继续。'))

    def _refreshSequences(self):
        records = self.controller.guidedCompatibleSequences()
        self.sequenceCombo.clear()
        for name, length, excess in records:
            self.sequenceCombo.addItem('%s — %d nt（多出并忽略：%d nt）' %
                                       (name, length, excess), name)
        if not records:
            self.sequenceCombo.addItem('请先生成并确认骨架链', None)

    def _applySequence(self):
        name = self.sequenceCombo.currentData()
        if not name:
            return
        result = self.controller.guidedApplySequence(name)
        self.sequenceStatus.setText(result.get('message', ''))
        if not result.get('ok'):
            QtWidgets.QMessageBox.warning(self, '智能引导设计',
                                          result.get('message', '无法应用序列。'))

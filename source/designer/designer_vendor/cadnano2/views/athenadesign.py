"""In-cadnano ATHENA wireframe design dialog."""

import math
import os

from cadnano2.data.dnasequences import sequences
from cadnano2.model.io.athena import (inspect_ply_dimension, preset_shapes,
                                      read_ply_preview_mesh,
                                      estimate_scaffold_minimum,
                                      recommended_engine, safe_name)
import cadnano2.util as util

util.qtWrapImport(
    'QtCore', globals(), ['QDir', 'QPointF', 'QRectF', 'Qt'])
util.qtWrapImport(
    'QtGui', globals(), ['QBrush', 'QColor', 'QPainter', 'QPen'])
util.qtWrapImport(
    'QtWidgets', globals(), [
        'QComboBox', 'QDialog', 'QDialogButtonBox', 'QFileDialog',
        'QFormLayout', 'QHBoxLayout', 'QLabel', 'QLineEdit', 'QPushButton',
        'QMessageBox', 'QScrollArea', 'QSpinBox', 'QVBoxLayout', 'QWidget'])


class WireframeCylinderPreview(QWidget):
    """Lightweight, mouse-rotatable target-geometry cylinder preview."""
    def __init__(self, parent=None):
        super(WireframeCylinderPreview, self).__init__(parent)
        self._vertices = []
        self._edges = []
        self._edgeType = "DX"
        self._dimension = "3D"
        self._title = ""
        self._error = ""
        self._yaw = -0.65
        self._pitch = 0.42
        self._zoom = 1.0
        self._lastPosition = None
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        self.setToolTip(
            "Drag to rotate, use the mouse wheel to zoom, and "
            "double-click to reset the view.")

    @staticmethod
    def _subtract(first, second):
        return tuple(first[index] - second[index] for index in range(3))

    @staticmethod
    def _add(first, second):
        return tuple(first[index] + second[index] for index in range(3))

    @staticmethod
    def _scale(vector, value):
        return tuple(component * value for component in vector)

    @staticmethod
    def _cross(first, second):
        return (first[1] * second[2] - first[2] * second[1],
                first[2] * second[0] - first[0] * second[2],
                first[0] * second[1] - first[1] * second[0])

    @staticmethod
    def _normalise(vector):
        length = math.sqrt(sum(component * component
                               for component in vector))
        if length <= 1e-12:
            return (0.0, 0.0, 0.0)
        return tuple(component / length for component in vector)

    def _resetView(self):
        if self._dimension == "2D":
            self._yaw = 0.0
            self._pitch = 0.0
        else:
            self._yaw = -0.65
            self._pitch = 0.42
        self._zoom = 1.0

    def setMesh(self, vertices, faces, dimension, edge_type, title):
        mins = [min(vertex[axis] for vertex in vertices)
                for axis in range(3)]
        maxs = [max(vertex[axis] for vertex in vertices)
                for axis in range(3)]
        centre = tuple((mins[axis] + maxs[axis]) * 0.5
                       for axis in range(3))
        centred = [self._subtract(vertex, centre) for vertex in vertices]
        radius = max(math.sqrt(sum(component * component
                                   for component in vertex))
                     for vertex in centred)
        radius = max(radius, 1e-9)
        normalised = [self._scale(vertex, 1.0 / radius)
                      for vertex in centred]
        edges = set()
        for face in faces:
            for first, second in zip(face, face[1:] + face[:1]):
                edges.add(tuple(sorted((first, second))))
        geometry_changed = (normalised != self._vertices or
                            sorted(edges) != self._edges)
        self._vertices = normalised
        self._edges = sorted(edges)
        self._dimension = str(dimension)
        self._edgeType = str(edge_type).upper()
        self._title = str(title)
        self._error = ""
        if geometry_changed:
            self._resetView()
        self.update()

    def setError(self, message):
        self._error = str(message)
        self._vertices = []
        self._edges = []
        self.update()

    def _rotate(self, point):
        cos_yaw, sin_yaw = math.cos(self._yaw), math.sin(self._yaw)
        x_value = cos_yaw * point[0] + sin_yaw * point[2]
        z_value = -sin_yaw * point[0] + cos_yaw * point[2]
        cos_pitch = math.cos(self._pitch)
        sin_pitch = math.sin(self._pitch)
        y_value = cos_pitch * point[1] - sin_pitch * z_value
        depth = sin_pitch * point[1] + cos_pitch * z_value
        return x_value, y_value, depth

    def _project(self, point, scale, centre_x, centre_y):
        rotated = self._rotate(point)
        return (QPointF(centre_x + rotated[0] * scale,
                        centre_y - rotated[1] * scale), rotated[2])

    def _cylinders(self):
        """Return one visual cylinder per PLY edge, independent of DX/6HB."""
        result = []
        for first_index, second_index in self._edges:
            result.append((self._vertices[first_index],
                           self._vertices[second_index], 0))
        return result

    def paintEvent(self, unused_event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor("#fafbfc"))
        painter.setPen(QPen(QColor("#d2d8df"), 1.0))
        painter.setBrush(QBrush(QColor("#ffffff")))
        painter.drawRoundedRect(
            QRectF(1, 1, self.width() - 2, self.height() - 2), 8, 8)
        if self._error:
            painter.setPen(QColor("#a13030"))
            painter.drawText(
                QRectF(24, 24, self.width() - 48, self.height() - 48),
                Qt.AlignmentFlag.AlignCenter,
                "Preview unavailable\n%s" % self._error)
            painter.end()
            return
        if not self._vertices:
            painter.end()
            return

        painter.setPen(QColor("#637889"))
        painter.drawText(
            QRectF(10, 7, self.width() - 20, 18),
            Qt.AlignmentFlag.AlignCenter,
            "Drag: rotate   •   Wheel: zoom   •   Double-click: reset")
        scale = min(self.width() - 54, self.height() - 84) * 0.47 * self._zoom
        centre_x = self.width() * 0.5
        centre_y = (self.height() - 32) * 0.5 + 7
        cylinders = []
        for first, second, cylinder_index in self._cylinders():
            point_a, depth_a = self._project(
                first, scale, centre_x, centre_y)
            point_b, depth_b = self._project(
                second, scale, centre_x, centre_y)
            cylinders.append(((depth_a + depth_b) * 0.5,
                              point_a, point_b, cylinder_index))
        cylinders.sort(key=lambda item: item[0])
        colours = ["#4f91bd"]
        cylinder_width = 6.0
        for unused_depth, point_a, point_b, cylinder_index in cylinders:
            outline = QPen(QColor("#204a68"), cylinder_width + 2.0)
            outline.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(outline)
            painter.drawLine(point_a, point_b)
            body = QPen(QColor(colours[cylinder_index % len(colours)]),
                        cylinder_width)
            body.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(body)
            painter.drawLine(point_a, point_b)

        projected_vertices = [self._project(
            vertex, scale, centre_x, centre_y) for vertex in self._vertices]
        projected_vertices.sort(key=lambda item: item[1])
        node_radius = 4.2
        painter.setPen(QPen(QColor("#204a68"), 1.2))
        painter.setBrush(QBrush(QColor("#d9ebf7")))
        for point, unused_depth in projected_vertices:
            painter.drawEllipse(point, node_radius, node_radius)

        painter.setPen(QColor("#263746"))
        painter.drawText(
            QRectF(12, self.height() - 43, self.width() - 24, 36),
            Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap,
            self._title)
        painter.end()

    def mousePressEvent(self, event):
        self._lastPosition = (event.position()
                              if hasattr(event, "position") else event.pos())
        self.setCursor(Qt.CursorShape.ClosedHandCursor)
        event.accept()

    def mouseMoveEvent(self, event):
        if self._lastPosition is None:
            return
        position = (event.position()
                    if hasattr(event, "position") else event.pos())
        delta = position - self._lastPosition
        self._lastPosition = position
        self._yaw += float(delta.x()) * 0.012
        self._pitch += float(delta.y()) * 0.012
        self._pitch = max(-1.5, min(1.5, self._pitch))
        self.update()
        event.accept()

    def mouseReleaseEvent(self, event):
        self._lastPosition = None
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        event.accept()

    def wheelEvent(self, event):
        steps = float(event.angleDelta().y()) / 120.0
        self._zoom *= math.pow(1.12, steps)
        self._zoom = max(0.35, min(4.0, self._zoom))
        self.update()
        event.accept()

    def mouseDoubleClickEvent(self, event):
        self._resetView()
        self.update()
        event.accept()


class AthenaDesignDialog(QDialog):
    def __init__(self, parent=None):
        super(AthenaDesignDialog, self).__init__(parent)
        self.setWindowTitle("Wireframe Design")
        self.setMinimumSize(680, 520)
        self.resize(760, 760)
        self._presets = preset_shapes()
        self._customPath = None
        self._projectPathWasEdited = False

        outerLayout = QVBoxLayout(self)
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        contentWidget = QWidget(self.scrollArea)
        layout = QVBoxLayout(contentWidget)
        self.scrollArea.setWidget(contentWidget)
        outerLayout.addWidget(self.scrollArea, 1)

        intro = QLabel(
            "Choose a bundled paper/example shape or import an ASCII PLY "
            "wireframe. cadnano detects 2D/3D, recommends the compatible "
            "ATHENA engine, and opens the generated JSON automatically.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.shapeBox = QComboBox(self)
        for preset in self._presets:
            self.shapeBox.addItem(
                "%s — %s" % (preset["dimension"], preset["label"]),
                preset["path"])
        form.addRow("Shape:", self.shapeBox)

        self.previewLabel = WireframeCylinderPreview(self)
        self.previewLabel.setFixedSize(500, 280)
        self.previewLabel.setToolTip(
            "Drag to rotate, use the mouse wheel to zoom, and "
            "double-click to reset the cylinder preview.")

        sourceRow = QHBoxLayout()
        self.pathEdit = QLineEdit(self)
        self.pathEdit.setReadOnly(True)
        self.browseShapeButton = QPushButton("Import PLY…", self)
        sourceRow.addWidget(self.pathEdit, 1)
        sourceRow.addWidget(self.browseShapeButton)
        form.addRow("Input mesh:", sourceRow)

        self.dimensionLabel = QLabel(self)
        form.addRow("Detected shape:", self.dimensionLabel)

        self.edgeTypeBox = QComboBox(self)
        self.edgeTypeBox.addItem("DX / 2HB — flexible, lightweight", "DX")
        self.edgeTypeBox.addItem("6HB — rigid", "6HB")
        form.addRow("Edge structure:", self.edgeTypeBox)

        self.engineLabel = QLabel(self)
        form.addRow("Recommended program:", self.engineLabel)

        self.edgeLengthSpin = QSpinBox(self)
        self.edgeLengthSpin.setRange(21, 420)
        self.edgeLengthSpin.setValue(42)
        self.edgeLengthSpin.setSuffix(" bp")
        form.addRow("Shortest edge length:", self.edgeLengthSpin)

        self.scaffoldBox = QComboBox(self)
        for name in sorted(sequences):
            sequence = "".join(char for char in sequences[name].upper()
                               if char in "ACGT")
            self.scaffoldBox.addItem(
                "%s — %d nt" % (name, len(sequence)), name)
        form.addRow("Scaffold template:", self.scaffoldBox)
        scaffoldHelp = QLabel(
            "The selected sequence supplies the length required by ATHENA. "
            "The generated cadnano JSON is opened without applied sequence; "
            "you can assign the final scaffold later with Seq.")
        scaffoldHelp.setWordWrap(True)
        form.addRow("", scaffoldHelp)

        self.scaffoldCapacityLabel = QLabel(self)
        self.scaffoldCapacityLabel.setWordWrap(True)
        self.scaffoldCapacityLabel.setMinimumHeight(78)
        self.scaffoldCapacityLabel.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.scaffoldCapacityLabel.setMargin(10)

        self.nameEdit = QLineEdit(self)
        form.addRow("Design name:", self.nameEdit)

        projectRow = QHBoxLayout()
        self.projectEdit = QLineEdit(self)
        self.browseProjectButton = QPushButton("Choose…", self)
        projectRow.addWidget(self.projectEdit, 1)
        projectRow.addWidget(self.browseProjectButton)
        form.addRow("Project folder:", projectRow)
        layout.addLayout(form)

        capacityHeading = QLabel("Capacity check:", self)
        layout.addWidget(capacityHeading)
        layout.addWidget(self.scaffoldCapacityLabel)

        previewHeading = QLabel("Shape preview:", self)
        previewRow = QHBoxLayout()
        previewRow.addStretch(1)
        previewRow.addWidget(self.previewLabel)
        previewRow.addStretch(1)
        layout.addWidget(previewHeading)
        layout.addLayout(previewRow)

        outputHelp = QLabel(
            "Automatic output: JSON, input PLY, parameters, and cylindrical/"
            "routing/pseudo-atomic BILD models. Initial PDB/mmCIF, sequence "
            "CSV and scaffold_sequence files are not saved; use Wireframe "
            "3D Export after assigning the final sequence.")
        outputHelp.setWordWrap(True)
        layout.addWidget(outputHelp)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel, parent=self)
        self.buttons.button(
            QDialogButtonBox.StandardButton.Ok).setText(
                "Run Wireframe Design")
        outerLayout.addWidget(self.buttons)

        self.shapeBox.currentIndexChanged.connect(self._presetChanged)
        self.browseShapeButton.clicked.connect(self._browseShape)
        self.edgeTypeBox.currentIndexChanged.connect(self._edgeTypeChanged)
        self.edgeLengthSpin.valueChanged.connect(self._refreshPreview)
        self.edgeLengthSpin.valueChanged.connect(
            self._refreshScaffoldCapacity)
        self.scaffoldBox.currentIndexChanged.connect(
            self._refreshScaffoldCapacity)
        self.nameEdit.textChanged.connect(self._nameChanged)
        self.projectEdit.textEdited.connect(self._projectEdited)
        self.browseProjectButton.clicked.connect(self._browseProject)
        self.buttons.accepted.connect(self._acceptIfValid)
        self.buttons.rejected.connect(self.reject)
        self._presetChanged(0)

    def _currentPath(self):
        return self._customPath or str(self.shapeBox.currentData())

    def _presetChanged(self, unused_index):
        self._customPath = None
        path = str(self.shapeBox.currentData())
        self.pathEdit.setText(path)
        self.nameEdit.setText(
            safe_name(os.path.splitext(os.path.basename(path))[0]))
        self._refreshDimension()

    def _browseShape(self):
        selected = QFileDialog.getOpenFileName(
            self, "Import wireframe", QDir.homePath(),
            "PLY wireframe (*.ply)")
        if isinstance(selected, (tuple, list)):
            selected = selected[0]
        if not selected:
            return
        self._customPath = os.path.abspath(str(selected))
        self.pathEdit.setText(self._customPath)
        self.nameEdit.setText(safe_name(
            os.path.splitext(os.path.basename(self._customPath))[0]))
        self._refreshDimension()

    def _refreshDimension(self):
        try:
            dimension = inspect_ply_dimension(self._currentPath())
            self.dimensionLabel.setText(
                "%s wireframe" % ("2D planar" if dimension == "2D"
                                  else "3D spatial"))
            self.dimensionLabel.setProperty("athenaDimension", dimension)
            self._refreshPreview()
            self._refreshScaffoldCapacity()
        except Exception as error:
            self.dimensionLabel.setText("Invalid mesh: %s" % error)
            self.dimensionLabel.setProperty("athenaDimension", "")
            self._showPreviewError(str(error))
            self.scaffoldCapacityLabel.setText("Unavailable for invalid mesh.")
        self._updateEngine()

    def _updateEngine(self, unused_index=None):
        dimension = self.dimensionLabel.property("athenaDimension")
        edge_type = self.edgeTypeBox.currentData()
        if dimension:
            engine = recommended_engine(str(dimension), str(edge_type))
            shown = "DAEDALUS" if engine == "DAEDALUS2" else engine
            self.engineLabel.setText(shown)
            self.engineLabel.setProperty("athenaEngine", engine)
        else:
            self.engineLabel.setText("Unavailable")
            self.engineLabel.setProperty("athenaEngine", "")

    def _edgeTypeChanged(self, unused_index=None):
        self._updateEngine()
        self._refreshSuggestedProjectPath()
        self._refreshPreview()
        self._refreshScaffoldCapacity()

    def _refreshScaffoldCapacity(self, unused_value=None):
        try:
            minimum = estimate_scaffold_minimum(
                self._currentPath(), self.edgeTypeBox.currentData(),
                self.edgeLengthSpin.value())
            scaffold_name = str(self.scaffoldBox.currentData())
            sequence = "".join(
                char for char in sequences[scaffold_name].upper()
                if char in "ACGT")
            available = len(sequence)
        except Exception as error:
            self.scaffoldCapacityLabel.setText(
                "Capacity could not be estimated: %s" % error)
            self.scaffoldCapacityLabel.setStyleSheet(
                "QLabel { color: #8f1d1d; background: #fff1f0; "
                "border: 1px solid #e5a3a0; border-radius: 5px; }")
            self.scaffoldCapacityLabel.setProperty(
                "wireframeScaffoldMinimum", 0)
            return
        self.scaffoldCapacityLabel.setProperty(
            "wireframeScaffoldMinimum", minimum)
        if available < minimum:
            self.scaffoldCapacityLabel.setText(
                "Too short: selected scaffold has %d nt; this geometry "
                "requires at least about %d nt, plus possible unpaired "
                "vertex bases. Choose a longer scaffold or reduce the "
                "shortest edge length." % (available, minimum))
            self.scaffoldCapacityLabel.setStyleSheet(
                "QLabel { color: #8f1d1d; background: #fff1f0; "
                "border: 1px solid #e5a3a0; border-radius: 5px; "
                "font-weight: 600; }")
        else:
            self.scaffoldCapacityLabel.setText(
                "Selected scaffold: %d nt. Geometric minimum: about %d nt; "
                "the final routing may include additional unpaired vertex "
                "bases." % (available, minimum))
            self.scaffoldCapacityLabel.setStyleSheet(
                "QLabel { color: #1f6130; background: #eff8f1; "
                "border: 1px solid #9ac8a5; border-radius: 5px; }")

    def _previewTitle(self):
        if self._customPath:
            name = os.path.splitext(os.path.basename(self._customPath))[0]
        else:
            name = str(self.shapeBox.currentText()).split("—", 1)[-1].strip()
        dimension = str(
            self.dimensionLabel.property("athenaDimension") or "")
        edge_type = str(self.edgeTypeBox.currentData() or "DX").upper()
        edge_label = "DX / 2HB" if edge_type == "DX" else "6HB"
        return "%s   |   %s   |   %s   |   shortest edge %d bp" % (
            name, dimension, edge_label, self.edgeLengthSpin.value())

    def _showPreviewError(self, message):
        self.previewLabel.setError(message)

    def _refreshPreview(self, unused_value=None):
        """Render the selected PLY as an interactive cylinder preview."""
        try:
            vertices, faces = read_ply_preview_mesh(self._currentPath())
            dimension = str(
                self.dimensionLabel.property("athenaDimension") or
                inspect_ply_dimension(self._currentPath()))
        except Exception as error:
            self._showPreviewError(str(error))
            return

        self.previewLabel.setMesh(
            vertices, faces, dimension, self.edgeTypeBox.currentData(),
            self._previewTitle())

    def _outputName(self):
        return "%s_%s" % (
            safe_name(self.nameEdit.text()),
            str(self.edgeTypeBox.currentData()).upper())

    def _refreshSuggestedProjectPath(self):
        if not self._projectPathWasEdited:
            self.projectEdit.setText(os.path.join(
                QDir.homePath(), "Desktop", self._outputName()))

    def _nameChanged(self, value):
        self._refreshSuggestedProjectPath()

    def _projectEdited(self, unused_value):
        self._projectPathWasEdited = True

    def _browseProject(self):
        parent = QFileDialog.getExistingDirectory(
            self, "Choose wireframe project parent folder",
            os.path.dirname(str(self.projectEdit.text())) or QDir.homePath())
        if not parent:
            return
        self._projectPathWasEdited = True
        self.projectEdit.setText(os.path.join(
            str(parent), self._outputName()))

    def _acceptIfValid(self):
        if not self._currentPath() or not os.path.isfile(self._currentPath()):
            self.pathEdit.setFocus()
            return
        if not self.dimensionLabel.property("athenaDimension"):
            return
        if not str(self.projectEdit.text()).strip():
            self.projectEdit.setFocus()
            return
        minimum = int(self.scaffoldCapacityLabel.property(
            "wireframeScaffoldMinimum") or 0)
        scaffold_name = str(self.scaffoldBox.currentData())
        available = len("".join(
            char for char in sequences[scaffold_name].upper()
            if char in "ACGT"))
        if minimum and available < minimum:
            QMessageBox.warning(
                self, "Scaffold too short",
                "The selected scaffold has %d nt, but this geometry needs "
                "at least about %d nt. Choose a longer scaffold or reduce "
                "the shortest edge length." % (available, minimum))
            self.scaffoldBox.setFocus()
            return
        self.accept()

    def spec(self):
        scaffold_name = str(self.scaffoldBox.currentData())
        return {
            "ply_path": self._currentPath(),
            "dimension": str(
                self.dimensionLabel.property("athenaDimension")),
            "edge_type": str(self.edgeTypeBox.currentData()),
            "engine": str(self.engineLabel.property("athenaEngine")),
            "edge_length": int(self.edgeLengthSpin.value()),
            "scaffold_name": scaffold_name,
            "scaffold_sequence": sequences[scaffold_name],
            "name": safe_name(str(self.nameEdit.text())),
            "project_root": os.path.abspath(
                str(self.projectEdit.text()).strip())}

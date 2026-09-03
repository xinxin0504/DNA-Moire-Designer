from cadnano2.views import styles
import cadnano2.util as util
util.qtWrapImport('QtCore', globals(), ['QRectF', 'Qt'])
util.qtWrapImport('QtGui', globals(), ['QBrush', 'QFont'])
util.qtWrapImport('QtWidgets', globals(),  ['QColorDialog',
                                            'QGraphicsItem',
                                            'QGraphicsSimpleTextItem'])

_font = QFont(styles.thefont, 12, QFont.Weight.Bold)


class ColorPanel(QGraphicsItem):
    _scafColors = styles.scafColors
    _stapColors = styles.stapColors
    _pen = Qt.PenStyle.NoPen

    def __init__(self, parent=None, master=None):
        super(ColorPanel, self).__init__(parent)
        self._master = master if master is not None else self
        if master is None:
            self._peers = [self]
        else:
            master._peers.append(self)
        self.rect = QRectF(0, 0, 30, 30)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIgnoresTransformations)
        if master is None:
            self._scafColorIndex = -1  # painttool cycles to 0
            self._stapColorIndex = -1
            self._scafColor = self._scafColors[self._scafColorIndex]
            self._stapColor = self._stapColors[self._stapColorIndex]
            self._scafBrush = QBrush(self._scafColor)
            self._stapBrush = QBrush(self._stapColor)
        self._initLabel()
        self.hide()

    def _initLabel(self):
        self._label = label = QGraphicsSimpleTextItem("scaf\nstap", parent=self)
        label.setPos(32, 0)
        label.setFont(_font)
        # label.setBrush(_labelbrush)
        # label.hide()

    def boundingRect(self):
        return self.rect

    def paint(self, painter, option, widget=None):
        master = self._master
        painter.setPen(self._pen)
        painter.setBrush(master._scafBrush)
        painter.drawRect(0, 0, 30, 15)
        painter.setBrush(master._stapBrush)
        painter.drawRect(0, 15, 30, 15)

    def nextColor(self):
        master = self._master
        master._stapColorIndex += 1
        if master._stapColorIndex == len(self._stapColors):
            master._stapColorIndex = 0
        master._stapColor = self._stapColors[master._stapColorIndex]
        master._stapBrush.setColor(master._stapColor)
        master._updatePeers()

    def prevColor(self):
        self._master._stapColorIndex -= 1

    def color(self):
        return self._master._stapColor

    def scafColorName(self):
        return self._master._scafColor.name()

    def stapColorName(self):
        return self._master._stapColor.name()

    def changeScafColor(self):
        self._master._updatePeers()

    def changeStapColor(self):
        self._master._updatePeers()

    def _dialogParent(self):
        views = self.scene().views() if self.scene() is not None else []
        return views[0].window() if views else None

    def _updatePeers(self):
        for panel in self._peers:
            panel.update()

    def mousePressEvent(self, event):
        event.accept()
        self.setFocus()
        master = self._master
        if event.pos().y() < 15:
            newColor = QColorDialog.getColor(
                master._scafColor, self._dialogParent(),
                "Select scaffold color")
            if newColor.isValid() and \
                    newColor.name() != master._scafColor.name():
                master._scafColor = newColor
                master._scafBrush = QBrush(newColor)
                if not newColor in self._scafColors:
                    self._scafColors.insert(
                        master._scafColorIndex, newColor)
                master._updatePeers()
        else:
            newColor = QColorDialog.getColor(
                master._stapColor, self._dialogParent(),
                "Select staple color")
            if newColor.isValid() and \
                    newColor.name() != master._stapColor.name():
                master._stapColor = newColor
                master._stapBrush = QBrush(newColor)
                if not newColor in self._stapColors:
                    self._stapColors.insert(
                        master._stapColorIndex, newColor)
                master._updatePeers()

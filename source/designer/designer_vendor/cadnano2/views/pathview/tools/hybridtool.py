from .abstractpathtool import AbstractPathTool


class HybridTool(AbstractPathTool):
    """Two-click connector for free endpoints in different lattice parts."""
    def __init__(self, controller):
        super(HybridTool, self).__init__(controller)
        self._pending = None

    def __repr__(self):
        return "hybridTool"

    def setActive(self, willBeActive, oldTool=None):
        super(HybridTool, self).setActive(willBeActive, oldTool)
        if not willBeActive:
            self._pending = None

    def endpointClicked(self, strand, end):
        document = strand.document()
        if not document.isHybrid():
            self._window.statusBar().showMessage(
                "Cross-lattice connections are available only in Hybrid mode.",
                6000)
            return
        connection = (strand.connection3p() if end == '3p'
                      else strand.connection5p())
        if connection is not None:
            if connection.part() is not strand.part():
                self._pending = None
                if document.removeHybridConnection(strand, end):
                    self._window.statusBar().showMessage(
                        "Cross-lattice connection removed.", 6000)
                else:
                    self._window.statusBar().showMessage(
                        "Cross-lattice connection could not be removed.",
                        6000)
                return
            self._window.statusBar().showMessage(
                "Choose a free 5' or 3' endpoint.", 6000)
            return
        if self._pending is None:
            self._pending = (strand, end)
            lattice = 'Honeycomb' if strand.part()._step == 21 else 'Square'
            self._window.statusBar().showMessage(
                "%s %s endpoint selected; choose the complementary endpoint "
                "in the other lattice." % (lattice, end), 10000)
            return
        firstStrand, firstEnd = self._pending
        self._pending = None
        created, error = document.createHybridConnection(
            firstStrand, firstEnd, strand, end)
        if created:
            self._window.statusBar().showMessage(
                "Cross-lattice connection created.", 6000)
        else:
            self._window.statusBar().showMessage(error, 10000)

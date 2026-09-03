"""Regression check: Guided Design ignores scaffold sequence excess."""

from PyQt6.QtWidgets import QApplication

qtApp = QApplication.instance() or QApplication([])

from cadnano2.controllers import documentcontroller as module


class FakeOligo:
    def __init__(self, length):
        self._length = length
        self.applied = None

    def length(self):
        return self._length

    def applySequence(self, sequence):
        self.applied = sequence


class FakeController:
    def __init__(self, oligo):
        self._oligo = oligo

    def _guidedScaffoldOligos(self):
        return [self._oligo]


oligo = FakeOligo(6)
controller = FakeController(oligo)
module.guidedSequences['长度测试'] = 'AACCGGTTTT'
result = module.DocumentController.guidedApplySequence(controller, '长度测试')
assert result['ok']
assert oligo.applied == 'AACCGG'
assert '4 nt 已忽略' in result['message']
assert not hasattr(controller, '_guidedLoopInsertion')
print('guided sequence length smoke: OK')

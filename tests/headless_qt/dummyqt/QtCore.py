class QObject(object):
    """Minimal QObject stand-in used only by headless file-I/O tests."""
    def __init__(self, *args, **kwargs):
        self._dummy_parent = kwargs.get('parent')

    def setParent(self, parent):
        self._dummy_parent = parent

    def deleteLater(self):
        pass

class Qt(object):
    pass

class pyqtSignal(object):
    def __init__(self, *args):
        """ We don't actually do anything with argtypes because
        the real Qt will perform checks in the Gui version of
        cadnano which should suffice. """
        self.signalEmitters = {}
        self.argtypes = args
    def __get__(self, emitter, owner=None):
        if emitter is None:
            return self
        signal = self.signalEmitters.get(emitter)
        if signal is None:
            signal = self.signalEmitters[emitter] = pyqtBoundSignal()
        return signal

class pyqtBoundSignal(object):
    def __init__(self):
        self.targets = set()
    def connect(self, target):
        self.targets.add(target)
    def disconnect(self, target):
        self.targets.remove(target)
    def emit(self, *args):
        for t in self.targets:
            t(*args)

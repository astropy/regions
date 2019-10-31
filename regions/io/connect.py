from astropy.io.registry import UnifiedReadWrite, read, write


class ShapeListRead(UnifiedReadWrite):

    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'read')

    def __call__(self, *args, **kwargs):
        return read(self._cls, *args, **kwargs)


class ShapeListWrite(UnifiedReadWrite):

    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'write')

    def __call__(self, *args, **kwargs):
        write(self._instance, *args, **kwargs)

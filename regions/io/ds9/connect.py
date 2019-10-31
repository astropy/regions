from astropy.io import registry

from .read import read_ds9
from .write import write_ds9
from ..core import ShapeList


DS9_SIGNATURE = '# Region file format: DS9'


def is_ds9(origin, path, fileobj, *args, **kwargs):
    if fileobj is not None:
        pos = fileobj.tell()
        sig = fileobj.read(len(DS9_SIGNATURE))
        fileobj.seek(pos)
        return sig == DS9_SIGNATURE or sig == DS9_SIGNATURE.encode()
    else:
        return path is not None and path.lower().endswith((
            '.ds9', '.reg', '.ds9.gz', '.reg.gz'))


registry.register_reader('ds9', ShapeList, read_ds9)
registry.register_writer('ds9', ShapeList, write_ds9)
registry.register_identifier('ds9', ShapeList, is_ds9)

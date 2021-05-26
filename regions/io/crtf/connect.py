# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io import registry

from ..core import ShapeList
from .read import read_crtf
from .write import write_crtf


CRTF_SIGNATURE = '#CRTF'


def is_crtf(origin, path, fileobj, *args, **kwargs):  # noqa pylint: disable=unused-argument
    """
    Identify a CRTF region file.

    Parameters
    ----------
    origin : {'read', 'write'}
        A string identifying whether the file is to be opened for
        reading or writing.

    path : str
        The path to the file.

    fileobj : file-like or `None`
        An open file object to read the file's contents, or `None` if
        the file could not be opened.

    *args : tuple
        Positional arguments for the ``read`` or ``write`` function.

    **kwargs : dict
        Keyword arguments for the ``read`` or ``write`` function.

    Returns
    -------
    result : bool
        Returns `True` if the given file is a CRTF region file.
    """
    if fileobj is not None:
        pos = fileobj.tell()
        sig = fileobj.read(len(CRTF_SIGNATURE))
        fileobj.seek(pos)
        return sig == CRTF_SIGNATURE or sig == CRTF_SIGNATURE.encode()
    else:
        return (path is not None
                and path.lower().endswith(('.crtf', '.crtf.gz')))


registry.register_reader('crtf', ShapeList, read_crtf)
registry.register_writer('crtf', ShapeList, write_crtf)
registry.register_identifier('crtf', ShapeList, is_crtf)

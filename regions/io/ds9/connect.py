# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = []

DS9_SIGNATURE = '# Region file format: DS9'


def is_ds9(origin, path, fileobj, *args, **kwargs):  # noqa pylint: disable=unused-argument
    """
    Identify a DS9 region file.

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
        Returns `True` if the given file is a DS9 region file.
    """
    if fileobj is not None:
        pos = fileobj.tell()
        sig = fileobj.read(len(DS9_SIGNATURE))
        fileobj.seek(pos)
        return sig == DS9_SIGNATURE or sig == DS9_SIGNATURE.encode()
    else:
        return (path is not None
                and path.lower().endswith(('.ds9', '.reg', '.ds9.gz',
                                           '.reg.gz')))

# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = []


def is_stcs(methodname, filepath):
    """
    Identify an STC-S region file.

    Note: STC-S files cannot be reliably auto-detected from content alone,
    as they lack unique file signatures and use keywords that appear in
    other region formats (DS9, CRTF, etc.). Auto-detection is only based
    on file extensions.

    Parameters
    ----------
    methodname : {'read', 'write'}
        The method name called that needs auto-identification.

    filepath : str
        The path to the file.

    Returns
    -------
    result : bool
        Returns `True` if the given file has an STC-S file extension.
    """
    all_exten = ('.stcs', '.stc', '.stcs.txt', '.stc.txt')
    exten = {'read': all_exten, 'write': all_exten[0:2]}

    if methodname in ['read', 'write']:
        if isinstance(filepath, str):
            return filepath.lower().endswith(exten[methodname])

    return False

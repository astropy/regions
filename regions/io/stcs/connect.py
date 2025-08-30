# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.data import get_readable_fileobj

from regions.core import Region, Regions
from regions.core.registry import RegionsRegistry

__all__ = []


@RegionsRegistry.register(Region, 'identify', 'stcs')
@RegionsRegistry.register(Regions, 'identify', 'stcs')
def is_stcs(methodname, filepath):
    """
    Identify an STC-S region file.

    Parameters
    ----------
    methodname : {'read', 'write'}
        The method name called that needs auto-identification.

    filepath : str
        The path to the file.

    Returns
    -------
    result : bool
        Returns `True` if the given file is an STC-S region file.
    """
    all_exten = ('.stcs', '.stc', '.stcs.txt', '.stc.txt')
    exten = {'read': all_exten, 'write': all_exten[0:2]}

    if methodname == 'write':
        return filepath.lower().endswith(exten[methodname])

    elif methodname == 'read':
        if (isinstance(filepath, str)
                and filepath.lower().endswith(exten[methodname])):
            return True
        else:
            # Check file content for STC-S keywords
            try:
                with get_readable_fileobj(filepath, encoding='utf-8') as fileobj:
                    # Read first few lines to check for STC-S content
                    content = fileobj.read(512)  # Read first 512 characters
                    if content:
                        # Look for common STC-S keywords
                        stcs_keywords = ['Circle', 'Ellipse', 'Box', 'Polygon', 'Position',
                                       'ICRS', 'FK5', 'FK4', 'GALACTIC', 'ECLIPTIC']
                        content_upper = content.upper()
                        return any(keyword.upper() in content_upper for keyword in stcs_keywords)
                    return False
            except (UnicodeDecodeError, OSError):
                return False

    else:
        return False

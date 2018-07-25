
import distutils.version as vers
import pytest

from numpy.testing import assert_allclose
import astropy.version as astrov
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table


from ..read import FITSRegionParser
from ..core import FITSRegionParserError

_ASTROPY_MINVERSION = vers.LooseVersion('1.1')
_ASTROPY_VERSION = vers.LooseVersion(astrov.version)

implemented_region_types = ('ellipse', 'circle', 'box', 'polygon', 'point',
                            'annulus', 'elliptannulus')


@pytest.mark.parametrize('filename', ['data/region.fits'])
def test_file_fits(filename):

    filename = get_pkg_data_filename(filename)
    table = Table.read(filename)

    shapes = FITSRegionParser(table, 'warn').shapes

    assert shapes[0].region_type == 'circle'
    assert shapes[1].region_type == 'rectangle'
    assert shapes[2].region_type == 'ellipse'

    assert_allclose(shapes[0].coord[:2], list(table['X', 'Y'][0]))
    assert_allclose(shapes[1].coord[:2], list(table['X', 'Y'][1]))
    assert_allclose(shapes[2].coord[:2], list(table['X', 'Y'][2]))

    assert_allclose(shapes[0].coord[2:], [table['R'][0][0]])
    assert_allclose(shapes[1].coord[2:], list(table['R'][1]) + [table['ROTANG'][1]])
    assert_allclose(shapes[2].coord[2:], list(table['R'][2]) + [table['ROTANG'][2]])

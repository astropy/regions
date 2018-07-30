
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


@pytest.mark.parametrize('filename', ['data/fits_region_sample.fits'])
def test_file_fits(filename):

    filename = get_pkg_data_filename(filename)
    table = Table.read(filename)
    table['X'] = table['X'].astype(list)
    table['Y'] = table['Y'].astype(list)
    table.add_row([[1, 2, 3, 4], [5, 6, 7, 8], 'polygon', 0, 0, 8])

    shapes = FITSRegionParser(table, 'warn').shapes

    assert shapes[0].region_type == 'circle'
    assert shapes[1].region_type == 'rectangle'
    assert shapes[2].region_type == 'ellipse'
    assert shapes[3].region_type == 'rectangle'
    assert shapes[4].region_type == 'circleannulus'
    assert shapes[5].region_type == 'point'
    assert shapes[6].region_type == 'point'
    assert shapes[7].region_type == 'polygon'

    for x in range(7):
        assert_allclose(shapes[x].coord[:2], list(table['X', 'Y'][x]))

    assert_allclose(shapes[7].coord, [1, 5, 2, 6, 3, 7, 4, 8])

    assert_allclose(shapes[0].coord[2:], [table['R'][0][0]])
    for x in range(1, 5):
        assert_allclose(shapes[1].coord[2:], list(table['R'][1]) + [table['ROTANG'][1]])

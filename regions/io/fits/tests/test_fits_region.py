
import distutils.version as vers
import pytest

from numpy.testing import assert_allclose
import astropy.version as astrov
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from ..read import FITSRegionParser
from ..write import fits_region_objects_to_table
from ..core import FITSRegionParserError

from ....shapes.circle import CircleSkyRegion

_ASTROPY_MINVERSION = vers.LooseVersion('1.1')
_ASTROPY_VERSION = vers.LooseVersion(astrov.version)

implemented_region_types = ('ellipse', 'circle', 'box', 'polygon', 'point',
                            'annulus', 'elliptannulus')


@pytest.mark.parametrize('filename', ['data/region.fits'])
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
        assert_allclose([x.value for x in shapes[1].coord[2:]], list(table['R'][1]) + [table['ROTANG'][1]])

    regs = shapes.to_regions()
    table_ouput = fits_region_objects_to_table(regs)

    shape_ouput = FITSRegionParser(table_ouput).shapes

    for i in range(len(shapes)):
        assert shapes[i].region_type == shape_ouput[i].region_type
        assert shapes[i].coord == shape_ouput[i].coord
        assert shapes[i].meta == shape_ouput[i].meta


def test_only_pixel_regions():

    reg_sky = CircleSkyRegion(SkyCoord(1, 2, unit='deg'), 5 * u.deg)

    with pytest.raises(TypeError) as err:
        fits_region_objects_to_table([reg_sky])

    print(str(err))
    assert 'Every region must be a pixel region' in str(err)


def test_valid_columns():

    t = Table([[1, 2, 3]], names=('a'))

    with pytest.raises(ValueError) as err:
        FITSRegionParser(t)

    assert "This table has an invalid column name: 'a'" in str(err)


def test_valid_row():

    x = [1]
    y = [2]
    shapes = ['CIRCLE']
    t = Table([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    t['X'].unit = 'pix'
    t['Y'].unit = 'pix'

    with pytest.raises(ValueError) as err:
        FITSRegionParser(t)

    assert "The column: 'R' is missing in the table" in str(err)

    t[0]['SHAPE'] = 'PONT'

    with pytest.raises(ValueError) as err:
        FITSRegionParser(t)

    assert "'PONT' is not a valid FITS Region type" in str(err)

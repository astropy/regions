import pytest

from numpy.testing import assert_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

from ..read import FITSRegionParser, read_fits_region
from ..write import fits_region_objects_to_table
from ..core import FITSRegionParserError
from ...core import to_shape_list

from ....shapes.circle import CircleSkyRegion

implemented_region_types = ('ellipse', 'circle', 'box', 'polygon', 'point',
                            'annulus', 'elliptannulus')


@pytest.mark.parametrize('filename', ['data/fits_region.fits'])
def test_file_fits(filename):

    filename = get_pkg_data_filename(filename)
    table = Table.read(filename)

    shapes = FITSRegionParser(table, 'warn').shapes

    assert shapes[0].region_type == 'circle'
    assert shapes[1].region_type == 'rectangle'
    assert shapes[2].region_type == 'ellipse'
    assert shapes[3].region_type == 'rectangle'
    assert shapes[4].region_type == 'circleannulus'
    assert shapes[5].region_type == 'point'
    assert shapes[6].region_type == 'point'
    assert shapes[7].region_type == 'polygon'
    assert shapes[8].region_type == 'rectangle'

    for x in range(8):
        assert_allclose(shapes[x].coord[:2],
                        [table['X'][x][0], table['Y'][x][0]])

    assert_allclose(shapes[7].coord, [1, 5, 2, 6, 3, 7, 4, 8])

    assert_allclose(shapes[0].coord[2:], [table['R'][0][0]])

    for x in range(1, 4):
        assert_allclose([val.value for val in shapes[x].coord[2:]],
                        list(table['R'][x][:2]) + [table['ROTANG'][x]])

    regs = shapes.to_regions()
    table_ouput = fits_region_objects_to_table(regs)

    shape_ouput = FITSRegionParser(table_ouput).shapes

    for i in range(len(shapes)):
        assert shapes[i].region_type == shape_ouput[i].region_type
        assert shapes[i].coord == shape_ouput[i].coord
        assert shapes[i].meta == shape_ouput[i].meta

    # Reading the regions directly from file and converting to sky regions.
    regs_sky = read_fits_region(filename)
    with fits.open(filename) as pf:
        header = pf[1].header
        wcs = WCS(header, keysel=['image', 'binary', 'pixel'])
        regs_pix = [reg.to_pixel(wcs) for reg in regs_sky]
        shapes_roundtrip = to_shape_list(regs_pix, 'image')

    for i in range(len(shapes)):
        assert shapes[i].region_type == shapes_roundtrip[i].region_type
        assert_allclose(shapes[i].coord[:-1], shapes_roundtrip[i].coord[:-1])


def test_only_pixel_regions():

    reg_sky = CircleSkyRegion(SkyCoord(1, 2, unit='deg'), 5 * u.deg)

    with pytest.raises(TypeError) as err:
        fits_region_objects_to_table([reg_sky])

    assert 'Every region must be a pixel region' in str(err)


def test_valid_columns():

    t = Table([[1, 2, 3]], names=('a'))

    with pytest.raises(FITSRegionParserError) as err:
        FITSRegionParser(t)

    assert "This table has an invalid column name: 'a'" in str(err)


def test_valid_row():

    x = [1]
    y = [2]
    shapes = ['CIRCLE']
    t = Table([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    t['X'].unit = 'pix'
    t['Y'].unit = 'pix'

    with pytest.raises(FITSRegionParserError) as err:
        FITSRegionParser(t)

    assert "The column: 'R' is missing in the table" in str(err)

    t[0]['SHAPE'] = 'PONT'

    with pytest.raises(FITSRegionParserError) as err:
        FITSRegionParser(t)

    assert "'PONT' is not a valid FITS Region type" in str(err)

    t['ROTANG'] = [[20, 30]]
    t['ROTANG'].unit = 'deg'
    t[0]['SHAPE'] = 'PIE'

    with pytest.raises(FITSRegionParserError) as err:
        FITSRegionParser(t)

    assert "'PIE' is currently not supported in regions" in str(err)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.wcs import WCS
from numpy.testing import assert_allclose
import pytest

from ....shapes import CircleSkyRegion
from ...core import to_shape_list
from ..core import FITSRegionParserError
from ..read import FITSRegionParser, read_fits_region
from ..write import fits_region_objects_to_table

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
    table_output = fits_region_objects_to_table(regs)
    shape_output = FITSRegionParser(table_output).shapes

    for i, shape in enumerate(shapes):
        assert shape.region_type == shape_output[i].region_type
        assert shape.coord == shape_output[i].coord
        assert shape.meta == shape_output[i].meta

    # Reading the regions directly from file and converting to sky
    # regions.
    regs_sky = read_fits_region(filename)
    with fits.open(filename) as hdulist:
        header = hdulist[1].header
        wcs = WCS(header, keysel=['image', 'binary', 'pixel'])
        regs_pix = [reg.to_pixel(wcs) for reg in regs_sky]
        shapes_roundtrip = to_shape_list(regs_pix, 'image')

    for i, shape in enumerate(shapes):
        assert shape.region_type == shapes_roundtrip[i].region_type
        assert_allclose(shape.coord[:-1], shapes_roundtrip[i].coord[:-1])


def test_only_pixel_regions():
    reg_sky = CircleSkyRegion(SkyCoord(1, 2, unit='deg'), 5 * u.deg)

    with pytest.raises(TypeError) as excinfo:
        fits_region_objects_to_table([reg_sky])

    estr = 'Every region must be a pixel region'
    assert estr in str(excinfo.value)


def test_valid_columns():
    t = Table([[1, 2, 3]], names=('a'))

    with pytest.raises(FITSRegionParserError) as excinfo:
        FITSRegionParser(t)

    estr = 'This table has an invalid column name: "a"'
    assert estr in str(excinfo.value)


def test_valid_row():
    x = [1]
    y = [2]
    shapes = ['CIRCLE']
    t = Table([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    t['X'].unit = 'pix'
    t['Y'].unit = 'pix'

    with pytest.raises(FITSRegionParserError) as excinfo:
        FITSRegionParser(t)

    estr = 'The column "R" is missing in the table'
    assert estr in str(excinfo.value)

    t[0]['SHAPE'] = 'PONT'

    with pytest.raises(FITSRegionParserError) as excinfo:
        FITSRegionParser(t)

    estr = '"PONT" is not a valid FITS Region type'
    assert estr in str(excinfo.value)

    t['ROTANG'] = [[20, 30]]
    t['ROTANG'].unit = 'deg'
    t[0]['SHAPE'] = 'PIE'

    with pytest.raises(FITSRegionParserError) as excinfo:
        FITSRegionParser(t)

    estr = '"PIE" is currently not supported'
    assert estr in str(excinfo.value)

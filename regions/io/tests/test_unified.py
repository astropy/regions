# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.data import get_pkg_data_filename, get_pkg_data_fileobj
import pytest

from ...core import Regions


parametrize_filenames = pytest.mark.parametrize(
    'filename', ['../crtf/tests/data/CRTFgeneral.crtf',
                 '../ds9/tests/data/ds9.color.reg',
                 '../fits/tests/data/fits_region.fits'])


@parametrize_filenames
def test_read_filename(filename):
    filename = get_pkg_data_filename(filename)
    result = Regions.read(filename)
    assert isinstance(result, Regions)


# this does not work yet
@parametrize_filenames
@pytest.mark.parametrize('encoding',
                         ['binary',
                          pytest.param(None, marks=pytest.mark.xfail)])
def test_read_fileobj(filename, encoding):
    with get_pkg_data_fileobj(filename, encoding=encoding) as fileobj:
        result = Regions.read(fileobj)
    assert isinstance(result, Regions)


@pytest.mark.parametrize('file_type', ['crtf', 'ds9'])
def test_write_format(file_type, tmpdir):
    infilename = get_pkg_data_filename('../crtf/tests/data/CRTFgeneral.crtf')
    outfilename = str(tmpdir / 'region')
    shapelist = Regions.read(infilename)
    shapelist.write(outfilename, format=file_type)

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import pytest
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS

from regions.core import RegionMeta, RegionVisual
from regions.shapes.tests.test_common import BaseTestSphericalSkyRegion
from regions.shapes.whole_sky import WholeSphericalSkyRegion


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestWholeSphericalSkyRegion(BaseTestSphericalSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = WholeSphericalSkyRegion(meta=meta, visual=visual)
    inside = [(3 * u.deg, 4 * u.deg), (1 * u.deg, 3 * u.deg)]
    outside = []
    expected_repr = ('<WholeSphericalSkyRegion()>')
    expected_str = ('Region: WholeSphericalSkyRegion')

    def test_del(self):
        # test that Region shape attributes (descriptors) cannot be
        # deleted)
        if len(self.reg._params) > 0:
            with pytest.raises(AttributeError):
                delattr(self.reg, self.reg._params[0])

    def test_copy(self):
        reg = self.reg.copy()
        assert isinstance(reg, WholeSphericalSkyRegion)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):
        try:
            self.reg.to_sky(wcs)
        except NotImplementedError:
            pytest.xfail()

        try:
            self.reg.to_pixel(wcs)
        except NotImplementedError:
            pytest.xfail()

    def test_frame_transformation(self):
        reg2 = self.reg.transform_to('galactic')
        assert isinstance(reg2, WholeSphericalSkyRegion)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg

    def test_bounding_circle(self):
        bounding_circle = self.reg.bounding_circle
        assert bounding_circle is None

    def test_bounding_lonlat(self):
        bounding_lonlat = self.reg.bounding_lonlat
        assert bounding_lonlat is None

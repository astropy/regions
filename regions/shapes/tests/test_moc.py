# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

import pytest
import tempfile
import numpy as np
import copy
import os

from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import ICRS
from astropy.io import fits
from astropy_healpix import HEALPix

from ...core import PixCoord, RegionMask
from ..._utils.examples import make_example_dataset
from ...tests.helpers import make_simple_wcs

from ..moc import MOCSkyRegion

from .utils import HAS_MATPLOTLIB


@pytest.fixture()
def get_random_skycoords():
    def f(size):
        return SkyCoord(np.random.uniform(0, 360, size),
                        np.random.uniform(-90, 90, size),
                        unit='deg')

    return f


@pytest.fixture()
def get_random_lonlat():
    def f(size):
        return np.random.uniform(0, 360, size) * u.deg, \
               np.random.uniform(-90, 90, size) * u.deg

    return f


@pytest.fixture()
def mocs():
    moc1 = {'1': [0]}
    moc1_increased = {'1': [0, 1, 2, 17, 22]}
    moc2 = {'1': [30]}
    moc2_increased = {'1': [8, 30, 31, 28, 43]}

    return dict(moc1=MOCSkyRegion.from_json(moc1),
                moc1_increased=MOCSkyRegion.from_json(moc1_increased),
                moc2=MOCSkyRegion.from_json(moc2),
                moc2_increased=MOCSkyRegion.from_json(moc2_increased))


@pytest.fixture()
def mocs_op():
    moc1 = MOCSkyRegion.from_json({
        '0': [0, 2, 3, 4, 5]
    })
    moc2 = MOCSkyRegion.from_json({
        '0': [0, 1, 7, 4, 3]
    })
    moc3 = MOCSkyRegion.from_json({
        '0': [6, 3]
    })
    return dict(first=moc1, second=moc2, third=moc3)


@pytest.fixture()
def wcs():
    config = dict(crpix=(0, 0), crval=(0, 0), cdelt=(-5, 5), shape=(18, 36))
    dataset = make_example_dataset(config=config)
    return dataset.wcs

class TestMOC(object):
    """Test involving the manipulation of MOC objects"""

    def setup_class(self):
        self.galex = MOCSkyRegion.from_fits(get_pkg_data_filename('data/P-GALEXGR6-AIS-FUV.fits'))


    @pytest.mark.parametrize("size", [
        1000,
        10000,
        50000
    ])
    def test_from_skycoord(self, get_random_skycoords, size):
        skycoord = get_random_skycoords(size)
        moc = MOCSkyRegion.from_skycoord(skycoord, max_level=7)


    @pytest.mark.parametrize("size", [
        1000,
        10000,
        50000
    ])
    def test_from_lonlat(self, get_random_lonlat, size):
        lon, lat = get_random_lonlat(size)
        moc = MOCSkyRegion.from_lonlat(lon=lon, lat=lat, max_level=7)


    def test_from_fits(self):
        assert self.galex


    def test_write_and_from_json(self):
        # A dictionary of ('order', [ipix]) key-value pairs
        data = self.galex.serialize(format='json')
        moc_from_serialization = MOCSkyRegion.from_json(data)
        assert self.galex == moc_from_serialization


    def test_write_to_fits(self):
        hdulist = self.galex.serialize(format='fits')
        assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


    def test_write_to_json(self):
        data = self.galex.serialize(format='json')
        assert isinstance(data, dict)


    def test_contains(self):
        # Create a new WCS object. The number of axes must be set
        # from the start
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)

        level = 4
        size = 20
        ipix = np.random.randint(0, 12*4**level, size)
        npix = np.arange(12*4**level)
        ipix_out = np.setdiff1d(npix, ipix)

        moc = MOCSkyRegion.from_json({str(level): list(ipix)})

        hp = HEALPix(nside=(1 << level), order='nested', frame=ICRS())
        lon, lat = hp.healpix_to_lonlat(ipix)
        lon_out, lat_out = hp.healpix_to_lonlat(ipix_out)

        test_in = moc.contains(ra=lon, dec=lat)
        assert test_in.all()
        test_out = moc.contains(ra=lon_out, dec=lat_out)
        assert not test_out.any()


    def test_add_neighbours(self, mocs):
        assert mocs['moc1'].add_neighbours() == mocs['moc1_increased']
        assert mocs['moc2'].add_neighbours() == mocs['moc2_increased']


    def test_remove_neighbours(self, mocs):
        assert mocs['moc1_increased'].remove_neighbours() == mocs['moc1']
        assert mocs['moc2_increased'].remove_neighbours() == mocs['moc2']


    def test_neighbours(self, mocs):
        moc1 = copy.deepcopy(mocs['moc1'])
        moc2 = copy.deepcopy(mocs['moc2'])
        assert moc1.add_neighbours().remove_neighbours() == mocs['moc1']
        assert moc2.add_neighbours().remove_neighbours() == mocs['moc2']


    def test_sky_fraction(self):
        moc = MOCSkyRegion.from_json({
            '0': [0, 1, 2, 3, 4, 5]
        })
        assert moc.sky_fraction == 0.5


    def test_union(self, mocs_op):
        assert mocs_op['first'].union(mocs_op['second'], mocs_op['third']) == MOCSkyRegion.from_json({
            '0': [0, 1, 2, 3, 4, 5, 6, 7]
        })


    def test_intersection(self, mocs_op):
        assert mocs_op['first'].intersection(mocs_op['second'], mocs_op['third']) == MOCSkyRegion.from_json({
            '0': [3]
        })


    def test_difference(self, mocs_op):
        assert mocs_op['first'].difference(mocs_op['second'], mocs_op['third']) == MOCSkyRegion.from_json({
            '0': [2, 5]
        })


    def test_degrade_moc(self):
        precise_moc = MOCSkyRegion.from_json({
            '1': [4, 21]
        })
        degraded_moc = MOCSkyRegion.degrade_to_order(precise_moc, 0)
        assert degraded_moc == MOCSkyRegion.from_json({'0': [1, 5]})


    def test_complement(self):
        assert self.galex.complement().complement() == self.galex


    @pytest.mark.parametrize("pos_x, pos_y, expected_result", [
        (-6, 4, True),
        (0, 0, False),
        (0.6, 3.4, False),
    ])
    def test_contains_pix_moc(self, wcs, pos_x, pos_y, expected_result):
        pixel = PixCoord(pos_x, pos_y)

        reg = self.galex.to_pixel(wcs)
        assert reg.contains(pixel) == expected_result


    def test_to_sky(self, wcs):
        reg_px = self.galex.to_pixel(wcs)
        reg_cast_back_sky = reg_px.to_sky(wcs)

        assert reg_cast_back_sky == self.galex


    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_artist_moc(self, wcs):
        reg = self.galex.to_pixel(wcs)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': wcs})

        reg.plot(ax=ax, alpha=0.5, fill=True, color='r')

        plt.axis('equal')
        plt.xlabel('lon')
        plt.ylabel('lat')
        plt.xlim([0, 20])
        plt.ylim([0, 20])
        plt.grid(color="black", linestyle="dotted")


    def test_to_mask(self, wcs):
        reg_px = self.galex.to_pixel(wcs)
        mask = reg_px.to_mask()

        assert mask.shape == (33, 65)

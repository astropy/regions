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

from ..moc import MOC

def get_random_skycoords(size):
    return SkyCoord(ra=np.random.uniform(0, 360, size),
                    dec=np.random.uniform(-90, 90, size),
                    unit="deg")


@pytest.fixture()
def skycoords_gen_f():
    def gen_f(size):
        return SkyCoord(np.random.uniform(0, 360, size), np.random.uniform(-90, 90, size), unit='deg')

    return gen_f


@pytest.fixture()
def lonlat_gen_f():
    def gen_f(size):
        return np.random.uniform(0, 360, size) * u.deg, np.random.uniform(-90, 90, size) * u.deg

    return gen_f


@pytest.fixture()
def mocs():
    moc1 = {'1': [0]}
    moc1_increased = {'0': [0], '1': [17, 19, 22, 23, 35]}
    moc2 = {'1': [30]}
    moc2_increased = {'0': [7], '1': [8, 9, 25, 43, 41]}

    return dict(moc1=MOC.from_json(moc1),
                moc1_increased=MOC.from_json(moc1_increased),
                moc2=MOC.from_json(moc2),
                moc2_increased=MOC.from_json(moc2_increased))


@pytest.fixture()
def mocs_op():
    moc1 = MOC.from_json({
        '0': [0, 2, 3, 4, 5]
    })
    moc2 = MOC.from_json({
        '0': [0, 1, 7, 4, 3]
    })
    moc3 = MOC.from_json({
        '0': [6, 3]
    })
    return dict(first=moc1, second=moc2, third=moc3)


class TestMOC(object):
    """Test involving the manipulation of MOC objects"""

    def setup_class(self):
        self.moc_from_fits = MOC.from_fits(get_pkg_data_filename('data/P-GALEXGR6-AIS-FUV.fits'));


    @pytest.mark.parametrize("size", [
        1000,
        10000,
        50000
    ])
    def test_from_skycoords(self, skycoords_gen_f, size):
        skycoords = skycoords_gen_f(size)
        moc = MOC.from_skycoords(skycoords, max_norder=7)


    @pytest.mark.parametrize("size", [
        1000,
        10000,
        50000
    ])
    def test_from_lonlat(self, lonlat_gen_f, size):
        lon, lat = lonlat_gen_f(size)

        moc = MOC.from_lonlat(lon=lon, lat=lat, max_norder=7)


    def test_from_fits(self):
        assert self.moc_from_fits


    def test_write_and_from_json(self):
        tmp_file = tempfile.NamedTemporaryFile()
        self.moc_from_fits.write(tmp_file.name, format='json', write_to_file=True)

        with open(tmp_file.name, 'r') as moc_file:
            import json
            moc_d = json.load(moc_file)
            moc2 = MOC.from_json(json_moc=moc_d)
            assert self.moc_from_fits == moc2


    def test_write_to_fits(self):
        hdulist = self.moc_from_fits.write(format='fits')
        assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


    def test_write_to_json(self):
        moc_json = self.moc_from_fits.write(format='json')
        assert isinstance(moc_json, dict)


    def test_contains(self):
        order = 4
        size = 20
        healpix_arr = np.random.randint(0, 12*4**order, size)
        all_healpix_arr = np.arange(12*4**order)
        healpix_outside_arr = np.setdiff1d(all_healpix_arr, healpix_arr)

        moc = MOC.from_json(json_moc={str(order): list(healpix_arr)})

        hp = HEALPix(nside=(1 << order), order='nested', frame=ICRS())
        lon, lat = hp.healpix_to_lonlat(healpix_arr)
        lon_out, lat_out = hp.healpix_to_lonlat(healpix_outside_arr)

        should_be_inside_arr = moc.contains(ra=lon, dec=lat)
        assert should_be_inside_arr.all()
        should_be_outside_arr = moc.contains(ra=lon_out, dec=lat_out)
        assert not should_be_outside_arr.any()

        # test keep_inside field
        should_be_outside_arr = moc.contains(ra=lon, dec=lat, keep_inside=False)
        assert not should_be_outside_arr.any()
        should_be_inside_arr = moc.contains(ra=lon_out, dec=lat_out, keep_inside=False)
        assert should_be_inside_arr.all()


    def test_add_neighbours(self, mocs):
        mocs['moc1'].add_neighbours()
        assert mocs['moc1'] == mocs['moc1_increased']

        mocs['moc2'].add_neighbours()
        assert mocs['moc2'] == mocs['moc2_increased']


    def test_remove_neighbours(self, mocs):
        mocs['moc1_increased'].remove_neighbours()
        mocs['moc2_increased'].remove_neighbours()
        assert mocs['moc1_increased'] == mocs['moc1']
        assert mocs['moc2_increased'] == mocs['moc2']


    def test_neighbours(self, mocs):
        moc1 = copy.deepcopy(mocs['moc1'])
        moc2 = copy.deepcopy(mocs['moc2'])
        moc1.add_neighbours().remove_neighbours()
        moc2.add_neighbours().remove_neighbours()
        assert moc1 == mocs['moc1']
        assert moc2 == mocs['moc2']


    def test_sky_fraction(self):
        moc = MOC.from_json({
            '0': [0, 1, 2, 3, 4, 5]
        })
        assert moc.sky_fraction == 0.5


    def test_union(self, mocs_op):
        assert mocs_op['first'].union(mocs_op['second'], mocs_op['third']) == MOC.from_json({
            '0': [0, 1, 2, 3, 4, 5, 6, 7]
        })


    def test_intersection(self, mocs_op):
        assert mocs_op['first'].intersection(mocs_op['second'], mocs_op['third']) == MOC.from_json({
            '0': [3]
        })


    def test_difference(self, mocs_op):
        assert mocs_op['first'].difference(mocs_op['second'], mocs_op['third']) == MOC.from_json({
            '0': [2, 5]
        })


    def test_degrade_moc(self):
        precise_moc = MOC.from_json({
            '1': [4, 21]
        })
        degraded_moc = precise_moc.degrade_to_order(new_order=0)
        assert degraded_moc == MOC.from_json({'0': [1, 5]})


    def test_complement(self):
        assert self.moc_from_fits.complement().complement() == self.moc_from_fits


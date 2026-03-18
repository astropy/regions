# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure. It needs to live inside the package in order for it to
# get picked up when running the tests inside an interpreter using
# packagename.test

import astropy.units as u
import pytest
from astropy.wcs import WCS

from regions._utils.optional_deps import HAS_GWCS
from regions.tests.helpers import (WCS_CDELT_ARCSEC, WCS_CENTER,
                                   make_simple_wcs, make_sip_wcs)

try:
    from pytest_astropy_header.display import (PYTEST_HEADER_MODULES,
                                               TESTED_VERSIONS)
    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


def pytest_configure(config):
    if ASTROPY_HEADER:
        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the
        # list of packages for which version numbers are displayed when
        # running the tests.
        PYTEST_HEADER_MODULES.clear()
        deps = ['NumPy', 'Matplotlib', 'Astropy', 'Shapely']
        for dep in deps:
            PYTEST_HEADER_MODULES[dep] = dep.lower()

        from regions import __version__
        TESTED_VERSIONS['regions'] = __version__


@pytest.fixture
def simple_wcs():
    """Non-distorted TAN WCS aligned with the celestial axes."""
    return make_simple_wcs(WCS_CENTER, WCS_CDELT_ARCSEC * u.arcsec, 20)


@pytest.fixture
def rotated_wcs():
    """Non-distorted TAN WCS with a 25-degree rotation (CD matrix)."""
    return make_simple_wcs(WCS_CENTER, WCS_CDELT_ARCSEC * u.arcsec, 20,
                           rotation_deg=25.0)


@pytest.fixture
def sip_wcs():
    """TAN WCS with small SIP distortion terms."""
    return make_sip_wcs()


@pytest.fixture
def nonsquare_wcs():
    """Non-distorted TAN WCS with non-square pixels (0.03 x 0.05 deg)."""
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [10.5, 10.5]
    wcs.wcs.crval = [WCS_CENTER.ra.deg, WCS_CENTER.dec.deg]
    wcs.wcs.cdelt = [-0.03, 0.05]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    return wcs


@pytest.fixture
def gwcs_obj():
    """Return a simple GWCS object with a TAN projection (requires gwcs)."""
    if not HAS_GWCS:
        pytest.skip('gwcs is required')
    from regions.tests.helpers import make_gwcs
    return make_gwcs()

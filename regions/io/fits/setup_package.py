# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    parser_test = ['data/*.fits']
    return {'regions.io.fits.tests': parser_test}

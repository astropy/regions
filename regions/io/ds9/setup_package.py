# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    parser_test = ['data/*.reg']
    return {'regions.io.ds9.tests': parser_test}

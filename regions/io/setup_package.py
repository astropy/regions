# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

def get_package_data():
    parser_test = [os.path.join('data', 'ds9.fk5.reg')]
    return {'regions.io.tests': parser_test}

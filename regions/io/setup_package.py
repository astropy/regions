# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

def get_package_data():
    parser_test = [os.path.join('data', 'ds9.fk5.reg'),
                   os.path.join('data', 'ds9.fk5.strip.reg'),
                   os.path.join('data', 'ds9.fk5.hms.reg'),
                   os.path.join('data', 'ds9.fk5.hms.strip.reg'),
                   os.path.join('data', 'fk5_reference.reg'),
                   os.path.join('data', 'ds9.galactic.reg'),
                   os.path.join('data', 'ds9.galactic.strip.reg'),
                   os.path.join('data', 'ds9.galactic.hms.reg'),
                   os.path.join('data', 'ds9.galactic.hms.strip.reg'),
                   os.path.join('data', 'galactic_reference.reg')]

    return {'regions.io.tests': parser_test}

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
import string
import re
import copy
from collections import OrderedDict
from warnings import warn

from astropy import units as u
from astropy import coordinates
from astropy import log

from .core import *
from ..core import *




# All CASA files start with '#CRTF' . It may also include the version number like '#CRTFv0' .
regex_begin = re.compile(r'^#CRTFv?[\d]?$')

# Comment Format :
regex_comment = re.compile(r'^#.*$')

# Meta attributes of a region
regex_global = re.compile(r'^global\s+(?P<parameters>.*)?')

# Coordinate Format : "[x, y]"
regex_coordinate = re.compile(r'\[(?P<x>[\w.]*?)\s*[,]\s*(?P<y>[\w.]*?)\]')

# last length Format. For Ex : radius of a circle
regex_length = re.compile(r'(?:\[.*\])+[,]\s*([\w.]*)\]')

# region format
regex_region = re.compile(r'(?P<include>[+-])?(?P<regiontype>[a-z]*?)\[.*]')

# line format
regex_line = re.compile(r'(?P<region>[+-]?[a-z]*?\[.*\])(?:\s*[,]\s*(?P<parameters>.*))?')


def read_crtf(filename, errors='strict'):
    """Read a CRTF region file and a list of region objects.

    Parameters
    ----------
    filename : str
        The file path
    errors : ``warn``, ``ignore``, ``strict``
      The error handling scheme to use for handling parsing errors.
      The default is 'strict', which will raise a ``CRTFRegionParserError``.
      ``warn`` will raise a warning, and ``ignore`` will do nothing
      (i.e., be silent).

    Returns
    -------
    regions : list
        Python list of `regions.Region` objects.
    """
    with open(filename) as fh:
        if regex_begin.search(fh.readline()):
            region_string = fh.read()
            parser = CRTFParser(region_string, errors='strict')
            return parser.shapes.to_region()
        else :
            raise CRTFRegionParserError('Not in the given format')


class CRTFParser:



    # List of valid coordinate system
    # TODO : There are still many reference systems to support
    coordinate_systems = ['J2000', 'icrs', 'galactic', 'supergal', 'image', 'ecliptic']

    coordsys_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                coordinates.frame_transform_graph.get_names()))
    coordsys_mapping['J2000'] = 'fk5'
    coordsys_mapping['supergal'] = 'supergalactic'
    coordsys_mapping['ecliptic'] = 'geocentrictrueecliptic'

    def __init__(self, region_string, errors='strict'):
        if errors not in ('strict', 'ignore', 'warn'):
            msg = "``errors`` must be one of strict, ignore, or warn; is {}"
            raise ValueError(msg.format(errors))
        self.region_string = region_string
        self.errors = errors

        # Global states
        self.global_meta = {}

        # Results
        self.shapes = ShapeList()

    def __str__(self):
        ss = self.__class__.__name__
        ss += '\nErrors: {}'.format(self.errors)
        ss += '\nCoordsys: {}'.format(self.coordsys)
        ss += '\nGlobal meta: {}'.format(self.global_meta)
        ss += '\nShapes: {}'.format(self.shapes)
        ss += '\n'
        return ss

    def set_coordsys(self, coorsdysy):
        pass

    def parse_line(self, line):
        """Parse one line"""
        log.debug('Parsing {}'.format(line))
        # Skip blanks
        if line == '':
            return

        # Skip comments
        if regex_comment.search(line):
            return

        # Special case / header: parse global parameters into metadata
        global_parameters = regex_global.search(line)
        if global_parameters:
            self.global_meta = self.parse_meta(global_parameters.group('parameters'))
            return

        # Try to parse the line
        crtf_line = regex_line.search(line)
        if crtf_line:
            helper = CRTFRegionParser(self.global_meta, line)
            helper.parse()
            self.shapes.append(helper.shape)
        else:
            self._raise_error("Not a valid CRTF line '{0}'.".format(line))
            return

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, CRTFRegionParserWarning)
        elif self.errors == 'strict':
            raise CRTFRegionParserError(msg)

    def run(self):
        """Run all steps"""
        for line in self.region_string.split('\n'):
            self.parse_line(line)
            log.debug('Global state: {}'.format(self))

    def parse_meta(self, line):
        pass


class CRTFRegionParser:

    def __init__(self, global_meta, line):
        self.global_meta = global_meta
        self.line = line

        self.coordsys = self.global_meta.get('coordsys', 'image')
        self.meta_str = None
        self.coord_str = None
        self.type = False    # Contains where the definition is an annotation
        self.region_type = None
        self.meta = None
        self.shape = None
        self.include = None

    def parse(self):
        self.split_line()
        self.convert_coordinates()
        self.convert_meta()
        self.make_shape()
        log.debug(self)

    def split_line(self):
        pass

    def convert_coordinates(self):
        pass

    def convert_meta(self):
        pass

    def make_shape(self):
        pass






class CoordinateParser(object):

    @staticmethod
    def parse_angular_length_quantity(string_rep):
        """
        Given a string that is either a number or a number and a unit, return a
        Quantity of that string.  e.g.:

            23.9 -> 23.9*u.deg
            50" -> 50*u.arcsec
        """
        unit_mapping = {
            '"': u.arcsec,
            'arcsec':u.arcsec,
            "'": u.arcmin,
            'arcmin':u.arcmin,
            'deg':u.deg,
            'rad': u.rad,
            'pix': u.dimensionless_unscaled,
        }
        has_unit = string_rep[-1] not in string.digits
        if has_unit:
            unit = unit_mapping[string_rep[-1]]
            return u.Quantity(float(string_rep[:-1]), unit=unit)
        else:
            return u.Quantity(float(string_rep), unit=u.deg)




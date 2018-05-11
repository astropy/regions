# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
import re
import copy
import itertools
from warnings import warn

from astropy import units as u
from astropy import coordinates

__all__ = [
    'CRTFParser',
    'CRTFRegionParser',
]

from .core import CRTFRegionParserError, CRTFRegionParserWarning
from ..core import Shape, ShapeList

# All CASA files start with '#CRTF' . It may also include the version number like '#CRTFv0' .
regex_begin = re.compile(r'^#CRTFv?[\d]?$')

# Comment Format :
regex_comment = re.compile(r'^#.*$')

# Meta attributes of a region
regex_global = re.compile(r'^global\s+(?P<parameters>.*)?')

# Coordinate Format : "[x, y]"
regex_coordinate = re.compile(r'\[([\w.+-:]*?)\s*[,]\s*([\w.+-:]*?)\]')

# last length Format. For Ex : radius of a circle
regex_length = re.compile(r'(?:\[[^=]*\])+[,]\s*([\w.+-]*)\]')

regex_meta = re.compile(r'((?:\w+\s*=\s*[^,\[\]]+)|(?:\w+\s*=\s*\[.*?\]))')

# region format
regex_region = re.compile(r'(?P<include>[+-])?(?P<type>ann(?=\s))?(?P<regiontype>[a-z]*?)\[[^=]*]')

# line format
regex_line = re.compile(r'(?P<region>[+-]?(?:ann(?=\s))?[a-z]*?\[[^=]*\])(?:\s*[,]\s*(?P<parameters>.*))?')


def read_crtf(filename, errors='strict'):
    """
    Read a CRTF region file and a list of region objects.

    Parameters
    ----------
    filename : str
        The file path
    errors : ``warn``, ``ignore``, ``strict``
      The error handling scheme to use for handling parsing errors.
      The default is 'strict', which will raise a ``CRTFRegionParserError``.
      ``warn`` will raise a ``CRTFRegionParserWarning``, and ``ignore`` will do nothing
      (i.e., be silent).

    Returns
    -------
    regions : list
        Python list of `regions.Region` objects.
    """
    with open(filename) as fh:
        if regex_begin.search(fh.readline()):
            region_string = fh.read()
            parser = CRTFParser(region_string, errors)
            return parser.shapes.to_regions()
        else:
            raise CRTFRegionParserError('Every CRTF Region must start with "#CRTF" ')


class CRTFParser:

    valid_definition = ('box', 'centerbox', 'rotbox', 'poly', 'circle', 'annulus', 'ellipse',
                        'line', 'vector', 'text', 'symbol')

    valid_global_keys = ['coord', 'frame', 'corr', 'veltype', 'restfreq', 'linewidth', 'linestyle', 'symsize',
                         'symthick', 'color', 'font', 'fontsize', 'fontstyle', 'usetex', 'labelpos','labelcolor',
                         'labeloff', 'range']

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
        self.run()

    def __str__(self):
        ss = self.__class__.__name__
        ss += '\nErrors: {}'.format(self.errors)
        ss += '\nGlobal meta: {}'.format(self.global_meta)
        ss += '\nShapes: {}'.format(self.shapes)
        ss += '\n'
        return ss

    def parse_line(self, line):
        """
        Parse one line
        """

        # Skip blanks
        if line == '':
            return

        # Skip comments
        if regex_comment.search(line):
            return

        # Special case / header: parse global parameters into metadata
        global_parameters = regex_global.search(line)
        if global_parameters:
            self.parse_global_meta(global_parameters.group('parameters'))
            return

        # Try to parse the line
        crtf_line = regex_line.search(line)

        if crtf_line:
            region = regex_region.search(crtf_line.group('region'))
            type = region.group('type') or 'reg'
            include = region.group('include') or '+'
            region_type = region.group('regiontype')

            if region_type in self.valid_definition:
                helper = CRTFRegionParser(self.global_meta, include, type, region_type,
                                          *crtf_line.group('region', 'parameters'))
                helper.parse()
                self.shapes.append(helper.shape)
            else:
                self._raise_error("Not a valid CRTF Region type '{0}'.".format(region_type))

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
            self.parse_line(line.lower())


    def parse_global_meta(self, global_meta_str):

        if global_meta_str:
            global_meta_str = regex_meta.findall(global_meta_str + ',')
        if global_meta_str:
            for par in global_meta_str:
                par = par.strip()
                par = par.rstrip(',')
                par = par.split('=')
                val1 = par[0].strip()
                val2 = par[1].strip()
                if val1 in self.valid_global_keys :
                    if val1 in ('range', 'corr', 'labeloff'):
                        val2 = val2[1:-1]
                        val2 = val2.split(",")
                        val2 = [x.strip() for x in val2 if x]
                    self.global_meta[val1] = val2
                else:
                    self._raise_error('{0} is not a valid global meta key').format(val1)


class CRTFRegionParser:

    # List of valid coordinate system
    # TODO : There are still many reference systems to support

    coordinate_systems = ['j2000', 'icrs', 'galactic', 'supergal', 'image', 'ecliptic']

    coordsys_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                coordinates.frame_transform_graph.get_names()))
    coordsys_mapping['j2000'] = 'fk5'
    coordsys_mapping['supergal'] = 'supergalactic'
    coordsys_mapping['ecliptic'] = 'geocentrictrueecliptic'

    language_spec = {'circle': ['c', 'l'],
                     'box': ['c', 'c'],
                     'centerbox': ['c', 'pl'],
                     'rotbox': ['c', 'pl', 'l'],
                     'poly': itertools.cycle('c'),
                     'annulus': ['c', 'pl'],
                     'ellipse': ['c', 'pl', 'l'],
                     'line': ['c', 'c'],
                     'vector': ['c', 'c'],
                     'symbol': ['c', 's'],
                     'text': ['c', 's']}

    def __init__(self, global_meta, include, type, region_type, reg_str, meta_str, errors='strict'):

        self.global_meta = global_meta
        self.reg_str = reg_str
        self.meta_str = meta_str
        self.errors = errors

        self.coord = None
        self.coordsys = None
        self.coord_str = None
        self.type = type or 'reg'
        self.region_type = region_type
        self.meta = copy.deepcopy(global_meta)
        self.shape = None
        self.include = include or '+'

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, CRTFRegionParserWarning)
        elif self.errors == 'strict':
            raise CRTFRegionParserError(msg)

    def parse(self):

        self.convert_meta()
        self.coordsys = self.meta.get('coord', 'image')
        self.set_coordsys()
        self.convert_coordinates()
        self.make_shape()

    def set_coordsys(self):
        """
        Mapping to astropy's coordinate system name

        # TODO: needs expert attention (Most reference systems are not mapped)
        """
        if self.coordsys in self.coordsys_mapping:
            self.coordsys = self.coordsys_mapping[self.coordsys]

    def convert_coordinates(self):

        coord_list_str = regex_coordinate.findall(self.reg_str) + regex_length.findall(self.reg_str)
        coord_list = []

        if self.region_type == 'poly':
            if len(coord_list_str) < 4 or coord_list_str[0] != coord_list_str[-1]:
                self._raise_error('{} not in proper format'.format(self.reg_str))
        else:
            if len(coord_list_str) != len(self.language_spec[self.region_type]):
                self._raise_error('{} not in proper format'.format(self.reg_str))

        for x, y in zip(self.language_spec[self.region_type], coord_list_str):

            if x == 'c':
                if len(y) == 2:
                    coord_list.append(CoordinateParser.parse_coordinate(y[0]))
                    coord_list.append(CoordinateParser.parse_coordinate(y[1]))
                else:
                    self._raise_error('{} not in proper format'.format(self.reg_str))
            if x == 'pl':
                if len(y) == 2:
                    coord_list.append(CoordinateParser.parse_angular_length_quantity(y[0]))
                    coord_list.append(CoordinateParser.parse_angular_length_quantity(y[1]))
                else:
                    self._raise_error('{} not in proper format'.format(self.reg_str))
            if x == 'l':
                if isinstance(y, str):
                    coord_list.append(CoordinateParser.parse_angular_length_quantity(y))
                else:
                    self._raise_error('{} not in proper format'.format(self.reg_str))

        self.coord = coord_list

    def convert_meta(self):

        if self.meta_str:
            self.meta_str = regex_meta.findall(self.meta_str + ',')
        if self.meta_str:
            for par in self.meta_str:
                par = par.strip()
                par = par.rstrip(',')
                par = par.split("=")
                val1 = par[0].strip()
                val2 = par[1].strip()
                if val1 in CRTFParser.valid_global_keys or val1 == 'label':
                    if val1 in ('range', 'corr', 'labeloff'):
                        val2 = val2[1:-1]
                        val2 = val2.split(',')
                        val2 = [x.strip() for x in val2]
                    self.meta[val1] = val2
                else:
                    self._raise_error('{0} is not a valid meta key'.format(val1))
        self.meta['include'] = self.include

    def make_shape(self):
        """
        Make shape object
        """
        self.shape = Shape('CRTF', coordsys=self.coordsys,
                           region_type=self.region_type,
                           coord=self.coord,
                           meta=self.meta,
                           composite=False,
                           include=self.include != '-',
                           )


class CoordinateParser(object):

    @staticmethod
    def parse_coordinate(string_rep):
        """
        Parse a single coordinate
        """
        # Any CRTF coordinate representation (sexagesimal or degrees)

        if 'pix' in string_rep:
            return u.Quantity(string_rep[:-3], u.dimensionless_unscaled)
        if 'h' in string_rep or 'rad' in string_rep:
            return coordinates.Angle(string_rep)
        else:
            return coordinates.Angle(string_rep, u.deg)

    @staticmethod
    def parse_angular_length_quantity(string_rep):
        """
        Given a string that is a number and a unit, return a
        Quantity of that string.Raise an Error If there is no unit.  e.g.:
            50" -> 50*u.arcsec
        """
        unit_mapping = {
            '"': u.arcsec,
            "'": u.arcmin,
        }
        regex_str = re.compile(r'([0-9+-,.]*)(.*)')
        str = regex_str.search(string_rep)
        unit = str.group(2)
        if unit:
            if unit in unit_mapping:
                return u.Quantity(str.group(1), unit=unit_mapping[unit])
            return u.Quantity(str.group(1))
        else:
            raise CRTFRegionParserError('Units must be specified for {0} '.format(string_rep))

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import itertools
import re
import string
from warnings import warn

from astropy.coordinates import Angle, frame_transform_graph
import astropy.units as u
from astropy.utils.data import get_readable_fileobj

from ..core import Shape, ShapeList, reg_mapping
from .core import (DS9RegionParserError, DS9RegionParserWarning,
                   valid_symbols_ds9)

__all__ = ['read_ds9', 'DS9Parser', 'DS9RegionParser', 'CoordinateParser']

# Regular expression to extract region type or coodinate system
regex_global = re.compile('^#? *(-?)([a-zA-Z0-9]+)')

# Regular expression to extract meta attributes
regex_meta = re.compile(r'([a-zA-Z]+)(=)({.*?}|\'.*?\'|\".*?\"|[0-9\s]+\s?|[^=\s]+\s?[0-9]*)\s?')  # noqa

# Regular expression to strip parenthesis
regex_paren = re.compile('[()]')

# Regular expression to split coordinate strings
regex_splitter = re.compile('[, ]')


def read_ds9(filename, errors='strict'):
    """
    Read a DS9 region file in as a list of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The file path.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.DS9RegionParserError`. 'warn' will raise a
        `~regions.DS9RegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.

    Examples
    --------
    >>> from regions import read_ds9
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> file = get_pkg_data_filename('data/physical_reference.reg', package='regions.io.ds9.tests')
    >>> regs = read_ds9(file, errors='warn')
    >>> print(regs[0])
    Region: CirclePixelRegion
    center: PixCoord(x=330.0, y=1090.0)
    radius: 40.0
    """
    with get_readable_fileobj(filename) as fh:
        region_string = fh.read()

    parser = DS9Parser(region_string, errors=errors)
    return parser.shapes.to_regions()


class CoordinateParser:
    """
    Helper class to structure coordinate parser.
    """

    @staticmethod
    def parse_coordinate(string_rep, unit):
        """
        Parse a single coordinate.
        """
        # explicit radian ('r') value
        if string_rep[-1] == 'r':
            return Angle(string_rep[:-1], unit=u.rad)

        # explicit image ('i') and physical ('p') pixels
        elif string_rep[-1] in ['i', 'p']:
            return u.Quantity(string_rep[:-1]) - 1

        # Any ds9 coordinate representation (sexagesimal or degrees)
        elif 'd' in string_rep or 'h' in string_rep:
            return Angle(string_rep)

        elif unit == 'hour_or_deg':
            if ':' in string_rep:
                spl = tuple(float(x) for x in string_rep.split(':'))
                return Angle(spl, u.hourangle)
            else:
                ang = float(string_rep)
                return Angle(ang, u.deg)

        elif unit.is_equivalent(u.deg):
            # return Angle(string_rep, unit=unit)
            if ':' in string_rep:
                ang = tuple(float(x) for x in string_rep.split(':'))
            else:
                ang = float(string_rep)
            return Angle(ang, u.deg)

        elif unit.is_equivalent(u.dimensionless_unscaled):
            return u.Quantity(float(string_rep), unit) - 1

        else:
            return u.Quantity(float(string_rep), unit)

    @staticmethod
    def parse_angular_length_quantity(string_rep, unit=u.deg):
        """
        Parse a string into a Quantity object.

        Given a string that is a number and a unit, return a Quantity of
        that string, e.g.:

        * 23.9 -> 23.9 * u.deg
        * 50" -> 50 * u.arcsec
        """
        unit_mapping = {'"': u.arcsec,
                        "'": u.arcmin,
                        'd': u.deg,
                        'r': u.rad,
                        'i': u.dimensionless_unscaled,
                        'p': u.dimensionless_unscaled}
        has_unit = string_rep[-1] not in string.digits
        if has_unit:
            unit = unit_mapping[string_rep[-1]]
            return u.Quantity(float(string_rep[:-1]), unit=unit)
        else:
            return u.Quantity(float(string_rep), unit=unit)


class DS9Parser:
    """
    Parse a DS9 string.

    This class transforms a DS9 string to a
    `~regions.io.core.ShapeList`. The result is stored as ``shapes``
    attribute.

    Each line is tested for either containing a region type or a
    coordinate system. If a coordinate system is found the global
    coordsys state of the parser is modified. If a region type is found
    the `~regions.DS9RegionParser` is invokes to transform the line into
    a `~regions.Shape` object.

    Parameters
    ----------
    region_string : str
        A DS9 region string.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.DS9RegionParserError`. 'warn' will raise a
        `~regions.DS9RegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).

    Examples
    --------
    >>> from regions import DS9Parser
    >>> reg_str = 'image\\n circle(331.00,1091.00,40.00) # dashlist=8 3 select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 source=1 text={Circle} tag={foo} tag={foo bar} This is a Comment color=pink width=3 font="times 10 normal roman"'
    >>> regs = DS9Parser(reg_str, errors='warn').shapes.to_regions()
    >>> print(regs[0])
    Region: CirclePixelRegion
    center: PixCoord(x=330.0, y=1090.0)
    radius: 40.0
    """

    # List of valid coordinate system (all lowercase)
    coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs',
                          'physical', 'image', 'ecliptic', 'j2000']
    coordinate_systems += [f'wcs{letter}' for letter in string.ascii_lowercase]

    # Map to convert coordinate system names
    coordsys_mapping = dict(zip(frame_transform_graph.get_names(),
                                frame_transform_graph.get_names()))
    coordsys_mapping['ecliptic'] = 'geocentrictrueecliptic'
    coordsys_mapping['j2000'] = 'fk5'

    def __init__(self, region_string, errors='strict'):
        if errors not in ('strict', 'ignore', 'warn'):
            raise ValueError('errors must be one of "strict", "ignore", '
                             'or "warn"')
        self.region_string = region_string
        self.errors = errors

        # Global states
        self.coordsys = None
        self.global_meta = {}

        # Results
        self.shapes = ShapeList()

        self.run()

    def __str__(self):
        ss = self.__class__.__name__
        ss += f'\nErrors: {self.errors}'
        ss += f'\nCoordsys: {self.coordsys}'
        ss += f'\nGlobal meta: {self.global_meta}'
        ss += f'\nShapes: {self.shapes}'
        ss += '\n'
        return ss

    def set_coordsys(self, coordsys):
        """
        Transform coordinate system
        """
        # TODO: needs expert attention
        if coordsys in self.coordsys_mapping:
            self.coordsys = self.coordsys_mapping[coordsys]
        else:
            self.coordsys = coordsys

    def run(self):
        """
        Run all the steps.
        """
        for line_ in self.region_string.split('\n'):
            for line in line_.split(';'):
                self.parse_line(line)

    def parse_line(self, line):
        """
        Parse a single line.
        """
        # Skip blanks
        if line == '':
            return

        # Skip comments
        if line[0] == '#':
            return

        # Special case / header: parse global parameters into metadata
        if line.lstrip()[:6].lower() == 'global':
            self.global_meta = self.parse_meta(line.lower())
            # global_meta can specify "include=1"; never seen other options
            # used but presumably 0 means false
            include = not (self.global_meta.get('include')
                           in ('0', 'False', False))
            self.global_meta['include'] = include

            return

        # Try to parse the line
        region_type_search = regex_global.search(line.lower())
        if region_type_search:
            include = region_type_search.groups()[0]
            region_type = region_type_search.groups()[1]
        else:
            self._raise_error(f'No region type found for line "{line}".')
            return

        if region_type in self.coordinate_systems:
            # Found coord system definition
            self.set_coordsys(region_type)
            return
        if region_type not in DS9RegionParser.language_spec:
            self._raise_error(f'Region type "{region_type}" was found, but '
                              'it is not one of the supported region types.')
            return
        else:
            # Found region specification,
            region_end = region_type_search.span()[1]
            self.parse_region(include, region_type, region_end, line)

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, DS9RegionParserWarning)
        elif self.errors == 'strict':
            raise DS9RegionParserError(msg)

    @staticmethod
    def parse_meta(meta_str):
        """
        Parse the metadata for a single DS9 region string.

        Parameters
        ----------
        meta_str : str
            Meta string, the metadata is everything after the closing
            parenthesis of the region coordinate specification. All
            metadata is specified as key=value pairs separated by
            whitespace, but sometimes the values can also be whitespace
            separated.

        Returns
        -------
        meta : dict
            A dictionary containing the metadata.
        """
        # keys must be lower-casae, but data can be any case
        keys_vals = [(x.lower(), y) for x, _, y
                     in regex_meta.findall(meta_str.strip())]
        extra_text = regex_meta.split(meta_str.strip())[-1]
        result = {}
        for key, val in keys_vals:
            # regex can include trailing whitespace or inverted commas
            # remove it
            val = val.strip().strip("'").strip('"')
            if key == 'text':
                val = val.lstrip('{').rstrip('}')
            if key in result:
                if key == 'tag':
                    result[key].append(val)
                else:
                    raise ValueError(f'Duplicate key {key} found')
            else:
                if key == 'tag':
                    result[key] = [val]
                else:
                    result[key] = val
        if extra_text:
            result['comment'] = extra_text

        return result

    def parse_region(self, include, region_type, region_end, line):
        """
        Extract a Shape from a region string.
        """
        if self.coordsys is None:
            raise DS9RegionParserError('No coordinate system specified and a '
                                       'region has been found.')

        helper = DS9RegionParser(coordsys=self.coordsys, include=include,
                                 region_type=region_type,
                                 region_end=region_end,
                                 global_meta=self.global_meta, line=line)
        helper.parse()
        self.shapes.append(helper.shape)


# Global definitions to improve readability
radius = CoordinateParser.parse_angular_length_quantity
width = CoordinateParser.parse_angular_length_quantity
height = CoordinateParser.parse_angular_length_quantity
angle = CoordinateParser.parse_angular_length_quantity
coordinate = CoordinateParser.parse_coordinate


class DS9RegionParser:
    """
    Parse a DS9 region string.

    This will turn a line containing a DS9 region into a
    `~regions.Shape` object.

    Parameters
    ----------
    coordsys : str
        The coordinate system.

    include : {'', '-'}
        Flag at the beginning of the line

    region_type : str
        Region type.

    region_end : int
        Coordinate of the end of the region name. This is passed in
        order to handle whitespace correctly.

    global_meta : dict
        Global metadata.

    line : str
        The line to parse.
    """

    # Coordinate unit transformations
    coordinate_units = {'fk5': ('hour_or_deg', u.deg),
                        'fk4': ('hour_or_deg', u.deg),
                        'icrs': ('hour_or_deg', u.deg),
                        'geocentrictrueecliptic': (u.deg, u.deg),
                        'galactic': (u.deg, u.deg),
                        'physical': (u.dimensionless_unscaled,
                                     u.dimensionless_unscaled),
                        'image': (u.dimensionless_unscaled,
                                  u.dimensionless_unscaled),
                        'wcs': (u.dimensionless_unscaled,
                                u.dimensionless_unscaled)}

    for letter in string.ascii_lowercase:
        coordinate_units[f'wcs{letter}'] = (u.dimensionless_unscaled,
                                            u.dimensionless_unscaled)

    # DS9 language specification. This defines how a certain region is read.
    language_spec = {'point': (coordinate, coordinate),
                     'text': (coordinate, coordinate),
                     'circle': (coordinate, coordinate, radius),
                     # This is a special case to deal with n elliptical annuli
                     'ellipse': itertools.chain((coordinate, coordinate),
                                                itertools.cycle((radius,))),
                     'box': (coordinate, coordinate, width, height, angle),
                     'polygon': itertools.cycle((coordinate,)),
                     'line': (coordinate, coordinate, coordinate, coordinate),
                     'annulus': itertools.chain((coordinate, coordinate),
                                                itertools.cycle((radius,))),
                     }

    def __init__(self, coordsys, include, region_type, region_end,
                 global_meta, line):
        self.coordsys = coordsys
        self.include = include
        self.region_type = region_type
        self.region_end = region_end
        self.global_meta = global_meta
        self.line = line

        self.meta_str = None
        self.coord_str = None
        self.composite = None
        self.coord = None
        self.meta = None
        self.shape = None

    def __str__(self):
        ss = self.__class__.__name__
        ss += f'\nLine : {self.line}'
        ss += f'\nRegion end : {self.region_end}'
        ss += f'\nMeta string : {self.meta_str}'
        ss += f'\nCoord string: {self.coord_str}'
        ss += f'\nShape: {self.shape}'
        ss += '\n'
        return ss

    def parse(self):
        """
        Parse the region string.
        """
        self.parse_composite()
        self.split_line()
        self.convert_coordinates()
        self.convert_meta()
        self.make_shape()

    def parse_composite(self):
        """
        Determine whether the region is composite.
        """
        self.composite = '||' in self.line

    def split_line(self):
        """
        Split line into coordinates and meta string.
        """
        # coordinate of the # symbol or end of the line (-1) if not found
        hash_or_end = self.line.find('#')
        temp = self.line[self.region_end:hash_or_end].strip(' |')
        # force all coordinate names (circle, etc) to be lower-case
        self.coord_str = regex_paren.sub('', temp).lower()

        # don't want any meta_str if there is no metadata found
        if hash_or_end >= 0:
            self.meta_str = self.line[hash_or_end:]
        else:
            self.meta_str = ''

    def convert_coordinates(self):
        """
        Convert coordinate string to objects.
        """
        coord_list = []
        # strip out 'null' elements, i.e., ''. It might be possible to
        # eliminate these some other way, i.e., with regex directly.
        # We need to copy in order not to burn up the iterators.
        elements = [x for x in regex_splitter.split(self.coord_str) if x]
        element_parsers = self.language_spec[self.region_type]
        for ii, (element, element_parser) in enumerate(zip(elements,
                                                           element_parsers)):
            if element_parser is coordinate:
                unit = self.coordinate_units[self.coordsys][ii % 2]
                coord_list.append(element_parser(element, unit))
            elif (self.coordinate_units[self.coordsys][0]
                  is u.dimensionless_unscaled):
                coord_list.append(
                    element_parser(element, unit=u.dimensionless_unscaled))
            else:
                coord_list.append(element_parser(element))

        if self.region_type in ['ellipse', 'box'] and len(coord_list) % 2 == 1:
            coord_list[-1] = CoordinateParser.parse_angular_length_quantity(
                elements[len(coord_list) - 1])

        # Reset iterator for ellipse and annulus
        # Note that this cannot be done with copy.deepcopy on python2
        if self.region_type in ['ellipse', 'annulus']:
            self.language_spec[self.region_type] = itertools.chain(
                (coordinate, coordinate), itertools.cycle((radius,)))

        self.coord = coord_list

    def convert_meta(self):
        """
        Convert meta string to dict.
        """
        meta_ = DS9Parser.parse_meta(self.meta_str)
        self.meta = copy.deepcopy(self.global_meta)
        self.meta.update(meta_)
        # the 'include' is not part of the metadata string;
        # it is pre-parsed as part of the shape type and should always
        # override the global one
        self.include = (self.meta.get('include', True)
                        if self.include == '' else self.include != '-')
        self.meta['include'] = self.include

    def make_shape(self):
        """
        Make shape object.
        """
        # In DS9, ellipse can also represents an elliptical annulus
        # For elliptical annulus angle is optional.
        if self.region_type == 'ellipse':
            self.coord[2:] = [x * 2 for x in self.coord[2:]]
            if len(self.coord) % 2 == 1:  # This checks if angle is present
                self.coord[-1] /= 2

        if 'point' in self.meta:
            point = self.meta['point'].split(' ')
            if len(point) > 1:
                self.meta['symsize'] = point[1]
            self.meta['point'] = valid_symbols_ds9[point[0]]

        if 'font' in self.meta:
            fonts = self.meta['font'].split(' ')
            keys = ['font', 'fontsize', 'fontstyle', 'fontweight']
            for i, val in enumerate(fonts):
                self.meta[keys[i]] = val

        self.meta.pop('coord', None)

        self.shape = Shape(coordsys=self.coordsys,
                           region_type=reg_mapping['DS9'][self.region_type],
                           coord=self.coord, meta=self.meta,
                           composite=self.composite, include=self.include)

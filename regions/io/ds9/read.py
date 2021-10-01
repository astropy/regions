# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import itertools
import re
import string
from warnings import warn

from astropy.coordinates import Angle, frame_transform_graph
import astropy.units as u
from astropy.utils.data import get_readable_fileobj

from ...core import Regions
from ...core.registry import RegionsRegistry
from ..core import _Shape, _ShapeList, reg_mapping
from .core import (DS9RegionParserError, DS9RegionParserWarning,
                   valid_symbols_ds9)

__all__ = []

# Regular expression to extract region type or coordinate system
regex_global = re.compile('^#? *(-?)([a-zA-Z0-9]+)')

# Regular expression to extract meta attributes
regex_meta = re.compile(r'([a-zA-Z]+)(\s*)(=)(\s*)({.*?}|\'.*?\'|\".*?\"|'
                        r'[0-9\s]+\s?|[^=\s]+\s?[0-9]*)\s?')

# Regular expression to strip parenthesis
regex_paren = re.compile('[()]')

# Regular expression to split coordinate strings
regex_splitter = re.compile('[, ]')


@RegionsRegistry.register(Regions, 'read', 'ds9')
def _read_ds9(filename, errors='strict', cache=False):
    """
    Read a DS9 region file in as a list of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The filename of the file to access.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.DS9RegionParserError`. 'warn' will raise a
        `~regions.DS9RegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.
    """
    with get_readable_fileobj(filename, cache=cache) as fh:
        region_string = fh.read()
        return _parse_ds9(region_string, errors=errors)


@RegionsRegistry.register(Regions, 'parse', 'ds9')
def _parse_ds9(region_string, errors='strict'):
    parser = _DS9Parser(region_string, errors=errors)
    return Regions(parser.shapes.to_regions())


class _DS9CoordinateParser:
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


class _DS9Parser:
    """
    Parse a DS9 string.

    This class transforms a DS9 string to a
    `~regions.io.core.ShapeList`. The result is stored as ``shapes``
    attribute.

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
    """

    # Each line is tested for either containing a region type or a
    # coordinate system. If a coordinate system is found the global
    # coordsys state of the parser is modified. If a region type is
    # found the ``_DS9RegionParser`` is invokes to transform the line
    # into a `~regions.Shape` object.

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
        self.shapes = _ShapeList()

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
            # split on all semicolons, except those between {} braces
            # (ds9 text strings)
            for line in _split_semicolon(line_):
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
        if region_type not in _DS9RegionParser.language_spec:
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
        # keys must be lower-case, but data can be any case
        keys_vals = [(x.lower(), y) for x, _, _, _, y
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

        helper = _DS9RegionParser(coordsys=self.coordsys, include=include,
                                  region_type=region_type,
                                  region_end=region_end,
                                  global_meta=self.global_meta, line=line)
        helper.parse()
        self.shapes.append(helper.shape)


# Global definitions to improve readability
radius = _DS9CoordinateParser.parse_angular_length_quantity
width = _DS9CoordinateParser.parse_angular_length_quantity
height = _DS9CoordinateParser.parse_angular_length_quantity
angle = _DS9CoordinateParser.parse_angular_length_quantity
coordinate = _DS9CoordinateParser.parse_coordinate


class _DS9RegionParser:
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
                                                itertools.cycle((radius,)))}

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
        # index of the # symbol or end of the line (-1) if not found
        hash_idx = self.line.find('#')
        if hash_idx == -1:
            temp = self.line[self.region_end:].strip(' |')
            self.meta_str = ''  # no metadata found
        else:
            temp = self.line[self.region_end:hash_idx].strip(' |')
            self.meta_str = self.line[hash_idx:]

        # force all coordinate names (circle, etc) to be lower-case
        self.coord_str = regex_paren.sub('', temp).lower()

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
            coord_list[-1] = _DS9CoordinateParser.parse_angular_length_quantity(
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
        meta_ = _DS9Parser.parse_meta(self.meta_str)
        self.meta = copy.deepcopy(self.global_meta)
        self.meta.update(meta_)
        # the 'include' is not part of the metadata string;
        # it is pre-parsed as part of the shape type and should always
        # override the global one
        self.include = (self.meta.get('include', True)
                        if self.include == '' else self.include != '-')
        self.meta['include'] = self.include
        self.meta['default_style'] = 'ds9'

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
            keys = ['font', 'fontsize', 'fontweight', 'fontstyle']
            for i, val in enumerate(fonts):
                self.meta[keys[i]] = val

        self.meta.pop('coord', None)

        self.shape = _Shape(coordsys=self.coordsys,
                            region_type=reg_mapping['DS9'][self.region_type],
                            coord=self.coord, meta=self.meta,
                            composite=self.composite, include=self.include)


def _find_text_delim_idx(regstr):
    """
    Find the indices of the DS9 text field delimiters ({}, '', or "") in
    a string.
    """
    pattern = re.compile(r'(text\s*=\s*[{\'"])')
    idx0 = []
    delim = []
    start_idx = []
    for match in pattern.finditer(regstr):
        idx0.append(match.start())
        end_delim = match.group()[-1]
        if end_delim == '{':
            end_delim = '}'
        delim.append(end_delim)
        start_idx.append(match.span()[1])

    idx1 = []
    for sidx, char in zip(start_idx, delim):
        idx1.append(regstr.find(char, sidx))

    return idx0, idx1


def _split_semicolon(regstr):
    r"""
    Split a DS9 region string on semicolons.  The line is not split on
    semicolons found in a text field.

    This turned out to be a very trickly problem (attempts with regex failed)
      * text fields are delimited by {}, '', or ""
      * the text delimiters do not have to be consistent within a file
      * region strings can have unpaired ' (arcmin) or " (arcsec)
      * region strings can have "#" in a color value (cannot split on #)
      * the text field can contain "text={str}", "text='str'",
        'test="str"' etc., e.g., text={text="hello"}
      * the delimiter characters (escaped or not) used by the text field
        are not allowed within the text field, e.g., "text={my
        field\{test\}}" and "text='my field \'test\''" are invalid.
        However, "text={my field, 'test'}" is valid

    This code finds the text field delimiters and then finds the indices
    of the opening and closing delimiters. Text fields that contain
    "text={str}", etc. are included (as a smaller range between the
    larger text field range), but this is fine because all we need are
    index ranges where to exclude semicolons (for splitting). Semicolons
    found at indices between the open/close delimiter indices are
    excluded from splitting.
    """
    idx0, idx1 = _find_text_delim_idx(regstr)

    semi_idx = [pos for pos, char in enumerate(regstr) if char == ';']
    fidx = []
    for i in semi_idx:
        for i0, i1 in zip(idx0, idx1):
            if i >= i0 and i <= i1:
                break
        else:
            fidx.append(i + 1)
    fidx.insert(0, 0)

    return [regstr[i:j].rstrip(';') for i, j in zip(fidx, fidx[1:] + [None])]

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
import string
import itertools
import re
import copy
from astropy import units as u
from astropy import coordinates
from astropy.coordinates import BaseCoordinateFrame
from astropy import log
from astropy.utils.exceptions import AstropyUserWarning
from warnings import warn
from ..shapes import circle, rectangle, polygon, ellipse, point
from ..core import PixCoord

__all__ = [
    'read_ds9',
    'ds9_string_to_objects',
    'ds9_string_to_region_list',
    'ds9_region_list_to_objects',
    'DS9RegionParserWarning',
    'DS9RegionParserError',
]


class DS9RegionParserWarning(AstropyUserWarning):
    """
    A generic warning class for DS9 region parsing inherited from astropy's
    warnings
    """


class DS9RegionParserError(ValueError):
    """
    A generic error class for DS9 region parsing
    """


def read_ds9(filename, errors='strict'):
    """Read a ds9 region file in as a list of region objects.

    Parameters
    ----------
    filename : str
        The file path
    errors : ``warn``, ``ignore``, ``strict``
      The error handling scheme to use for handling parsing errors.
      The default is 'strict', which will raise a ``DS9RegionParserError``.
      ``warn`` will raise a warning, and ``ignore`` will do nothing
      (i.e., be silent).

    Returns
    -------
    regions : list
        Python list of `regions.Region` objects.
    """
    with open(filename) as fh:
        region_string = fh.read()

    return ds9_string_to_objects(region_string, errors=errors)


def parse_coordinate(string_rep, unit):
    """
    Parse a single coordinate
    """
    # Any ds9 coordinate representation (sexagesimal or degrees)
    if 'd' in string_rep or 'h' in string_rep:
        return coordinates.Angle(string_rep)
    elif unit is 'hour_or_deg':
        if ':' in string_rep:
            spl = tuple([float(x) for x in string_rep.split(":")])
            return coordinates.Angle(spl, u.hourangle)
        else:
            ang = float(string_rep)
            return coordinates.Angle(ang, u.deg)
    elif unit.is_equivalent(u.deg):
        # return coordinates.Angle(string_rep, unit=unit)
        if ':' in string_rep:
            ang = tuple([float(x) for x in string_rep.split(":")])
        else:
            ang = float(string_rep)
        return coordinates.Angle(ang, u.deg)
    else:
        return u.Quantity(float(string_rep), unit)


unit_mapping = {
    '"': u.arcsec,
    "'": u.arcmin,
    'r': u.rad,
    'i': u.dimensionless_unscaled,
}


def parse_angular_length_quantity(string_rep):
    """
    Given a string that is either a number or a number and a unit, return a
    Quantity of that string.  e.g.:

        23.9 -> 23.9*u.deg
        50" -> 50*u.arcsec
    """
    has_unit = string_rep[-1] not in string.digits
    if has_unit:
        unit = unit_mapping[string_rep[-1]]
        return u.Quantity(float(string_rep[:-1]), unit=unit)
    else:
        return u.Quantity(float(string_rep), unit=u.deg)


# these are the same function, just different names
radius = parse_angular_length_quantity
width = parse_angular_length_quantity
height = parse_angular_length_quantity
angle = parse_angular_length_quantity

# For the sake of readability in describing the spec, parse_coordinate etc. are renamed here
coordinate = parse_coordinate
language_spec = {'point': (coordinate, coordinate),
                 'circle': (coordinate, coordinate, radius),
                 # This is a special case to deal with n elliptical annuli
                 'ellipse': itertools.chain((coordinate, coordinate), itertools.cycle((radius,))),
                 'box': (coordinate, coordinate, width, height, angle),
                 'polygon': itertools.cycle((coordinate,)),
                 }

coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 'physical', 'image', 'ecliptic']
coordinate_systems += ['wcs{0}'.format(letter) for letter in string.ascii_lowercase]

coordsys_name_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                 coordinates.frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic'  # needs expert attention TODO

hour_or_deg = 'hour_or_deg'
coordinate_units = {'fk5': (hour_or_deg, u.deg),
                    'fk4': (hour_or_deg, u.deg),
                    'icrs': (hour_or_deg, u.deg),
                    'geocentrictrueecliptic': (u.deg, u.deg),
                    'galactic': (u.deg, u.deg),
                    'physical': (u.dimensionless_unscaled, u.dimensionless_unscaled),
                    'image': (u.dimensionless_unscaled, u.dimensionless_unscaled),
                    'wcs': (u.dimensionless_unscaled, u.dimensionless_unscaled),
                    }
for letter in string.ascii_lowercase:
    coordinate_units['wcs{0}'.format(letter)] = (u.dimensionless_unscaled, u.dimensionless_unscaled)

region_type_or_coordsys_re = re.compile("^#? *(-?)([a-zA-Z0-9]+)")

paren = re.compile("[()]")


def strip_paren(string_rep):
    return paren.sub("", string_rep)


def ds9_region_list_to_objects(region_list):
    """Given a list of parsed region tuples, product a list of astropy objects.

    TODO: show example what a "region list" is.

    Parameters
    ----------
    region_list : list
        List of TODO???

    Returns
    -------
    regions : list
        List of `regions.Region` objects
    """
    viz_keywords = ['color', 'dashed', 'width', 'point', 'font']

    output_list = []
    for region_type, coord_list, meta in region_list:

        # TODO: refactor, possible on the basis of # of parameters + sometimes handle corner cases

        if region_type == 'circle':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = circle.CircleSkyRegion(coord_list[0], coord_list[1])
            elif isinstance(coord_list[0], PixCoord):
                reg = circle.CirclePixelRegion(coord_list[0], coord_list[1])
            else:
                raise DS9RegionParserError("No central coordinate")
        elif region_type == 'ellipse':
            # Do not read elliptical annuli for now
            if len(coord_list) > 4:
                continue
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = ellipse.EllipseSkyRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            elif isinstance(coord_list[0], PixCoord):
                reg = ellipse.EllipsePixelRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            else:
                raise DS9RegionParserError("No central coordinate")
        elif region_type == 'polygon':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = polygon.PolygonSkyRegion(coord_list[0])
            elif isinstance(coord_list[0], PixCoord):
                reg = polygon.PolygonPixelRegion(coord_list[0])
            else:
                raise DS9RegionParserError("No central coordinate")
        elif region_type in ('rectangle', 'box'):
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = rectangle.RectangleSkyRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            elif isinstance(coord_list[0], PixCoord):
                reg = rectangle.RectanglePixelRegion(coord_list[0], coord_list[1], coord_list[2], coord_list[3])
            else:
                raise DS9RegionParserError("No central coordinate")
        elif region_type == 'point':
            if isinstance(coord_list[0], BaseCoordinateFrame):
                reg = point.PointSkyRegion(coord_list[0])
            elif isinstance(coord_list[0], PixCoord):
                reg = point.PointPixelRegion(coord_list[0])
            else:
                raise DS9RegionParserError("No central coordinate")
        else:
            # Note: this should effectively never happen, because it would
            # imply that the line_parser found a region that didn't match the
            # above region types.  However, this can help with development,
            # since we could in theory implement more region types in the line
            # parser and forget to add them here.
            warn("Skipping region with coords {0} because its type '{1}'"
                 " is not recognized."
                 .format(str(coord_list), region_type),
                 DS9RegionParserWarning
                 )
            continue
        reg.visual = {key: meta[key] for key in meta.keys() if key in viz_keywords}
        reg.meta = {key: meta[key] for key in meta.keys() if key not in viz_keywords}
        output_list.append(reg)
    return output_list


def ds9_string_to_objects(region_string, errors='strict'):
    """Parse ds9 region string to region objects

    Parameters
    ----------
    region_string : str
        DS9 region string
    errors : ``warn``, ``ignore``, ``strict``
      The error handling scheme to use for handling parsing errors.
      The default is 'strict', which will raise a ``ValueError``.
      ``warn`` will raise a warning, and ``ignore`` will do nothing
      (i.e., be silent).

    Returns
    -------
    regions : list
        List of `~regions.Region` objects

    Examples
    --------
    TODO
    """
    region_list = ds9_string_to_region_list(region_string, errors=errors)
    regions = ds9_region_list_to_objects(region_list)
    return regions


def ds9_string_to_region_list(region_string, errors='strict'):
    """Parse a DS9 region string.

    Parameters
    ----------
    region_string : str
        DS9 region string
    errors : ``warn``, ``ignore``, ``strict``
      The error handling scheme to use for handling skipped entries
      in a region file that were not parseable.
      The default is 'strict', which will raise a ``ValueError``.
      ``warn`` will raise a warning, and ``ignore`` will do nothing
      (i.e., be silent).

    Returns
    -------
    list of (region type, coord_list, meta, composite, include) tuples
    region_type : str
    coord_list : list of coordinate objects
    meta : metadata dict
    composite : bool
        indicates whether region is a composite region
    include : bool
        Whether the region is included (False -> excluded)
    """
    coordsys = None
    regions = []
    composite_region = None

    global_meta = {}

    # ds9 regions can be split on \n or ;
    lines = []
    for line_ in region_string.split('\n'):
        for line in line_.split(";"):
            lines.append(line)
            parsed = line_parser(line, coordsys, errors=errors)
            if parsed in coordinate_systems:
                coordsys = parsed
            elif parsed and (parsed[0] == 'global'):
                # set some global metadata from the 'global' header
                _, global_meta = parsed
            elif parsed:
                # local meta must override global meta
                meta = copy.copy(global_meta)
                region_type, coordlist, meta_, composite, include = parsed
                meta.update(meta_)
                meta['include'] = include
                log.debug("Region type = {0}.  Composite={1}"
                          .format(region_type, composite))
                if composite and composite_region is None:
                    composite_region = [(region_type, coordlist)]
                elif composite:
                    composite_region.append((region_type, coordlist))
                elif composite_region is not None:
                    composite_region.append((region_type, coordlist))
                    regions.append(composite_region)
                    composite_region = None
                else:
                    regions.append((region_type, coordlist, meta))

    return regions


def line_parser(line, coordsys=None, errors='strict'):
    """
    Parse a single ds9 region line into a string

    Parameters
    ----------
    line : str
        A single ds9 region contained in a string
    coordsys : str
        The global coordinate system name declared at the top of the ds9 file
    errors : ``warn``, ``ignore``, ``strict``
      The error handling scheme to use for handling skipped entries
      in a region file that were not parseable.
      The default is 'strict', which will raise a ``DS9RegionParserError``.
      ``warn`` will raise a warning, and ``ignore`` will do nothing
      (i.e., be silent).

    Returns
    -------
    (region_type, parsed_return, parsed_meta, composite, include)
    region_type : str
    coord_list : list of coordinate objects
    meta : metadata dict
    composite : bool
        indicates whether region is a composite region
    include : bool
        Whether the region is included (False -> excluded)
    """
    if errors not in ('strict', 'ignore', 'warn'):
        raise ValueError("``errors`` must be one of strict, ignore, or warn")

    if '# Region file format' in line and line[0] == '#':
        # This is just a file format line, we can safely skip it
        return

    # special case / header: parse global parameters into metadata
    if line.lstrip()[:6] == 'global':
        return global_parser(line)

    region_type_search = region_type_or_coordsys_re.search(line)
    if region_type_search:
        include = region_type_search.groups()[0]
        region_type = region_type_search.groups()[1]
    else:
        # if there's no line, it's just blank, so don't warn
        if line:
            # but otherwise, this should probably always raise a warning?
            # at least until we identify common cases for it
            warn("No region type found for line '{0}'.".format(line),
                 DS9RegionParserWarning)
        return

    if region_type in coordinate_systems:
        return region_type  # outer loop has to do something with the coordinate system information
    elif region_type in language_spec:
        if coordsys is None:
            raise DS9RegionParserError("No coordinate system specified and a"
                                       " region has been found.")

        if "||" in line:
            composite = True
        else:
            composite = False

        # end_of_region_name is the coordinate of the end of the region's name, e.g.:
        # circle would be 6 because circle is 6 characters
        end_of_region_name = region_type_search.span()[1]
        # coordinate of the # symbol or end of the line (-1) if not found
        hash_or_end = line.find("#")
        coords_etc = strip_paren(line[end_of_region_name:hash_or_end].strip(" |"))
        meta_str = line[hash_or_end:]

        parsed_meta = meta_parser(meta_str)

        if coordsys in coordsys_name_mapping:
            parsed = type_parser(coords_etc, language_spec[region_type],
                                 coordsys_name_mapping[coordsys])

            # Reset iterator for ellipse annulus
            if region_type == 'ellipse':
                language_spec[region_type] = itertools.chain((coordinate, coordinate), itertools.cycle((radius,)))

            parsed_angles = [(x, y)
                             for x, y in zip(parsed[:-1:2], parsed[1::2])
                             if (isinstance(x, coordinates.Angle) and
                                 isinstance(y, coordinates.Angle))
                            ]
            frame = coordinates.frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])

            lon, lat = zip(*parsed_angles)
            if hasattr(lon, '__len__') and hasattr(lon, '__lat__') and len(lon) == 1 and len(lat==1):
                # force entries to be scalar if they are length-1
                lon, lat = u.Quantity(lon[0]), u.Quantity(lat[0])
            else:
                # otherwise, they are vector quantitites
                lon, lat = u.Quantity(lon), u.Quantity(lat)
            sphcoords = coordinates.UnitSphericalRepresentation(lon, lat)
            coords = frame(sphcoords)

            return region_type, [coords] + parsed[len(coords) * 2:], parsed_meta, composite, include
        else:
            parsed = type_parser(coords_etc, language_spec[region_type],
                                 coordsys)
            if region_type == 'polygon':
                # have to special-case polygon in the phys coord case
                # b/c can't typecheck when iterating as in sky coord case
                coord = PixCoord(parsed[0::2], parsed[1::2])
                parsed_return = [coord]
            else:
                parsed = [_.value for _ in parsed]
                coord = PixCoord(parsed[0], parsed[1])
                parsed_return = [coord] + parsed[2:]

            # Reset iterator for ellipse annulus
            if region_type == 'ellipse':
                language_spec[region_type] = itertools.chain((coordinate, coordinate), itertools.cycle((radius,)))

            return region_type, parsed_return, parsed_meta, composite, include
    else:
        # This will raise warnings even if the first line is acceptable,
        # e.g. something like:
        # # Region file format: DS9 version 4.1
        # That behavior is unfortunate, but there's not a great workaround
        # except to let the user set `errors='ignore'`
        if errors in ('warn', 'strict'):
            message = ("Region type '{0}' was identified, but it is not one of "
                       "the known region types.".format(region_type))
            if errors == 'strict':
                raise DS9RegionParserError(message)
            else:
                warn(message, DS9RegionParserWarning)


def type_parser(string_rep, specification, coordsys):
    """
    For a given region line in which the type has already been determined,
    parse the coordinate definition

    Parameters
    ----------
    string_rep : str
        The string containing the coordinates.  For example, if your region is
        `circle(1,2,3)` this string would be `1,2,3`
    specification : iterable
        An iterable of coordinate specifications.  For example, for a circle,
        this would be a list of (coordinate, coordinate, radius).  Each
        individual specification should be a function that takes a string and
        returns the appropriate astropy object.  See ``language_spec`` for the
        definition of the grammar used here.
    coordsys : str
        The string name of the global coordinate system

    Returns
    -------
    coord_list : list
        The list of astropy coordinates and/or quantities representing radius,
        width, etc. for the region
    """
    coord_list = []
    splitter = re.compile("[, ]")
    # strip out "null" elements, i.e. ''.  It might be possible to eliminate
    # these some other way, i.e. with regex directly, but I don't know how.
    elements = [x for x in splitter.split(string_rep) if x]
    for ii, (element, element_parser) in enumerate(zip(elements, specification)):
        if element_parser is coordinate:
            unit = coordinate_units[coordsys][ii % 2]
            coord_list.append(element_parser(element, unit))
        else:
            coord_list.append(element_parser(element))

    return coord_list


# match an x=y pair (where y can be any set of characters) that may or may not
# be followed by another one
meta_token = re.compile("([a-zA-Z]+)(=)([^= ]+) ?")


# meta_spec = {'color': color,
#            }
# global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
# ruler(+175:07:14.900,+50:56:21.236,+175:06:52.643,+50:56:11.190) ruler=physical physical color=white font="helvetica 12 normal roman" text={Ruler}


def meta_parser(meta_str):
    """Parse the metadata for a single ds9 region string.

    The metadata is everything after the close-paren of the region coordinate specification.
    All metadata is specified as key=value pairs separated by whitespace, but
    sometimes the values can also be whitespace separated.
    """
    meta_token_split = [x for x in meta_token.split(meta_str.strip()) if x]
    equals_inds = [i for i, x in enumerate(meta_token_split) if x == '=']
    result = {meta_token_split[ii - 1]:
              " ".join(meta_token_split[ii + 1:jj - 1 if jj is not None else None])
              for ii, jj in zip(equals_inds, equals_inds[1:] + [None])}

    return result


def global_parser(line):
    return "global", meta_parser(line)

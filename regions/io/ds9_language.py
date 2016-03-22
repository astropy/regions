import string
import re
from astropy import units as u
from astropy import coordinates

def coordinate(string_rep, unit):
    # Any ds9 coordinate representation (sexagesimal or degrees)
    return coordinates.Angle(string_rep, unit=unit)

unit_mapping = {'"': u.arcsec,
                "'": u.arcmin}

def angular_length_quantity(string_rep):
    has_unit = string_rep[-1] not in string.digits
    if has_unit:
        unit = unit_mapping[string_rep[-1]]
        return u.Quantity(float(string_rep[:-1]), unit=unit)
    else:
        return u.Quantity(float(string_rep), unit=u.deg)

# these are the same function, just different names
radius = angular_length_quantity
width = angular_length_quantity
height = angular_length_quantity
angle = angular_length_quantity

language_spec = {'point': (coordinate, coordinate),
                 'circle': (coordinate, coordinate, radius),
                 'box': (coordinate, coordinate, width, height, angle),
                }

coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 'physical', 'image']
coordinate_systems += ['wcs{0}'.format(letter) for letter in string.ascii_lowercase]

coordinate_units = {'fk5': (u.hour, u.deg),
                    'fk4': (u.hour, u.deg),
                    'icrs': (u.hour, u.deg),
                    'galactic': (u.deg, u.deg),
                    'physical': (u.dimensionless_unscaled, u.dimensionless_unscaled),
                    'image': (u.dimensionless_unscaled, u.dimensionless_unscaled),
                    'wcs': (u.dimensionless_unscaled, u.dimensionless_unscaled),
                   }
for letter in string.ascii_lowercase:
    coordinate_units['wcs{0}'.format(letter)] = (u.dimensionless_unscaled, u.dimensionless_unscaled)

# circle(1.5, 3.6, 1.2)

region_type_or_coordsys_re = re.compile("#? *([a-zA-Z0-9]+)")

paren = re.compile("[()]")

def strip_paren(string_rep):
    return paren.sub("", string_rep)

def ds9_parser(filename):
    coordsys = None
    regions = []

    with open(filename,'r') as fh:
        for line in fh:
            parsed = line_parser(line, coordsys)
            if parsed in coordinate_systems:
                coordsys = parsed
            elif parsed:
                region_type, coordlist = parsed
                regions.append((region_type, coordlist))

    return regions

def line_parser(line, coordsys=None):
    region_type_search = region_type_or_coordsys_re.search(line)
    if region_type_search:
        region_type = region_type_search.groups()[0]
    else:
        return

    if region_type in coordinate_systems:
        return region_type # outer loop has to do something with the coordinate system information
    elif region_type in language_spec:
        if coordsys is None:
            raise ValueError("No coordinate system specified and a region has been found.")
        return region_type, type_parser(strip_paren(line[region_type_search.span()[1]:]),
                           language_spec[region_type], coordsys)

def type_parser(string_rep, specification, coordsys):
    coord_list = []
    splitter = re.compile("[, ]")
    for ii, (element, element_parser) in enumerate(zip(splitter.split(string_rep), specification)):
        if element_parser is coordinate:
            unit = coordinate_units[coordsys][ii]
            coord_list.append(element_parser(element, unit))
        else:
            coord_list.append(element_parser(element))

    return coord_list

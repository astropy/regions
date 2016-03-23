import string
import itertools
import re
from astropy import units as u
from astropy import coordinates

def coordinate(string_rep, unit):
    # Any ds9 coordinate representation (sexagesimal or degrees)
    if 'd' in string_rep or 'h' in string_rep:
        return coordinates.Angle(string_rep)
    elif unit is 'hour_or_deg':
        if ':' in string_rep:
            return coordinates.Angle(string_rep, unit=u.hour)
        else:
            return coordinates.Angle(string_rep, unit=u.deg)
    elif unit.is_equivalent(u.deg):
        return coordinates.Angle(string_rep, unit=unit)
    else:
        return u.Quantity(float(string_rep), unit)

unit_mapping = {'"': u.arcsec,
                "'": u.arcmin,
                'r': u.rad,
                'i': u.dimensionless_unscaled,
               }

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
                 'polygon': itertools.cycle((coordinate, )),
                }

coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 'physical', 'image', 'ecliptic']
coordinate_systems += ['wcs{0}'.format(letter) for letter in string.ascii_lowercase]

coordsys_name_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                 coordinates.frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic' # needs expert attention TODO

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

region_type_or_coordsys_re = re.compile("#? *([a-zA-Z0-9]+)")

paren = re.compile("[()]")

def strip_paren(string_rep):
    return paren.sub("", string_rep)

def ds9_parser(filename):
    """
    Parse a complete ds9 .reg file

    Returns
    -------
    list of (region type, coord_list) tuples
    """
    coordsys = None
    regions = []
    composite_region = None

    with open(filename,'r') as fh:
        for line_ in fh:
            # ds9 regions can be split on \n or ;
            for line in line_.split(";"):
                parsed = line_parser(line, coordsys)
                if parsed in coordinate_systems:
                    coordsys = parsed
                elif parsed:
                    region_type, coordlist, meta, composite = parsed
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
            coords = coordinates.SkyCoord([(x, y)
                                           for x, y in zip(parsed[:-1:2], parsed[1::2])
                                           if isinstance(x, coordinates.Angle) and
                                           isinstance(x, coordinates.Angle)], frame=coordsys_name_mapping[coordsys])
            return region_type, [coords] + parsed[len(coords)*2:], parsed_meta, composite
        else:
            return region_type, type_parser(coords_etc, language_spec[region_type],
                                            coordsys), parsed_meta, composite

def type_parser(string_rep, specification, coordsys):
    coord_list = []
    splitter = re.compile("[, ]")
    for ii, (element, element_parser) in enumerate(zip(splitter.split(string_rep), specification)):
        if element_parser is coordinate:
            unit = coordinate_units[coordsys][ii % 2]
            coord_list.append(element_parser(element, unit))
        else:
            coord_list.append(element_parser(element))

    return coord_list


# match an x=y pair (where y can be any set of characters) that may or may not
# be followed by another one
meta_token = re.compile("([a-zA-Z]+)(=)([^= ]+) ?")

#meta_spec = {'color': color,
#            }
# global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
# ruler(+175:07:14.900,+50:56:21.236,+175:06:52.643,+50:56:11.190) ruler=physical physical color=white font="helvetica 12 normal roman" text={Ruler}
    


def meta_parser(meta_str):
    meta_token_split = [x for x in meta_token.split(meta_str.strip()) if x]
    equals_inds = [i for i, x in enumerate(meta_token_split) if x is '=']
    result = {meta_token_split[ii-1]:
              " ".join(meta_token_split[ii+1:jj-1 if jj is not None else None])
              for ii,jj in zip(equals_inds, equals_inds[1:]+[None])}

    return result

if __name__ == "__main__":
    # simple tests for now...
    import glob
    results = {}
    for fn in glob.glob('/Users/adam/Downloads/tests/regions/*.reg'):
        print(fn)
        results[fn] = ds9_parser(fn)

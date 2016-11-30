# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy import coordinates
from .. import shapes

__all__ = [
    'write_ds9',
    'ds9_objects_to_string',
]

coordsys_name_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                 coordinates.frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic'  # needs expert attention TODO


def ds9_objects_to_string(regions, coordsys='fk5', fmt='.4f', radunit='deg'):
    """Convert list of regions ato ds9 region strings.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects

    Returns
    -------
    region_string : str
        ds9 region string


    Examples
    --------
    TODO
    """

    ids = {
        shapes.circle.CircleSkyRegion: 'skycircle',
        shapes.circle.CirclePixelRegion: 'pixcircle',
        shapes.ellipse.EllipseSkyRegion: 'skyellipse',
        shapes.ellipse.EllipsePixelRegion: 'pixellipse',
        shapes.polygon.PolygonSkyRegion: 'skypolygon',
        shapes.polygon.PolygonPixelRegion: 'pixpolygon',
    }

    if radunit == 'arcsec':
        if coordsys in coordsys_name_mapping.keys():
            radunitstr = '"'
        else:
            raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
    else:
        radunitstr = ''

    ds9_strings = {
        'circle': 'circle({x:' + fmt + '},{y:' + fmt + '},{r:' + fmt + '}' + radunitstr + ')\n',
        'ellipse': 'ellipse({x:' + fmt + '},{y:' + fmt + '},{r1:' + fmt + '}' + radunitstr + ',{r2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})\n',
        'polygon': 'polygon({c})\n',
    }

    output = '# Region file format: DS9 astropy/regions\n'
    output += '{}\n'.format(coordsys)

    # convert coordsys string to coordsys object
    if coordsys in coordsys_name_mapping:
        frame = coordinates.frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
    else:
        # for pixel/image/physical frames
        frame = None

    for reg in regions:
        if isinstance(reg, shapes.circle.CircleSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            r = float(reg.radius.to(radunit).value)
            output += ds9_strings['circle'].format(**locals())
        elif isinstance(reg, shapes.circle.CirclePixelRegion):
            x = reg.center.x
            y = reg.center.y
            r = reg.radius
            output += ds9_strings['circle'].format(**locals())
        elif isinstance(reg, shapes.ellipse.EllipseSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            r1 = float(reg.major.to(radunit).value)
            r2 = float(reg.minor.to(radunit).value)
            ang = float(reg.angle.to('deg').value)
            output += ds9_strings['ellipse'].format(**locals())
        elif isinstance(reg, shapes.ellipse.EllipsePixelRegion):
            x = reg.center.x
            y = reg.center.y
            r1 = reg.major
            r2 = reg.minor
            ang = reg.angle
            output += ds9_strings['ellipse'].format(**locals())
        elif isinstance(reg, shapes.polygon.PolygonSkyRegion):
            v = reg.vertices.transform_to(frame)
            coords = [(x.to('deg').value, y.to('deg').value) for x, y in
                      zip(v.spherical.lon, v.spherical.lat)]
            val = "{:" + fmt + "}"
            temp = [val.format(x) for _ in coords for x in _]
            c = ",".join(temp)
            output += ds9_strings['polygon'].format(**locals())
        elif isinstance(reg, shapes.polygon.PolygonPixelRegion):
            v = reg.vertices
            coords = [(x, y) for x, y in zip(v.x, v.y)]
            val = "{:" + fmt + "}"
            temp = [val.format(x) for _ in coords for x in _]
            c = ",".join(temp)
            output += ds9_strings['polygon'].format(**locals())

    return output


def write_ds9(regions, filename='ds9.reg', coordsys='fk5'):
    """Convert list of regions to ds9 string and write to file.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects
    filename : str
        Filename
    coordsys : {TODO}
        Coordinate system
    """
    output = ds9_objects_to_string(regions, coordsys)
    with open(filename, 'w') as fh:
        fh.write(output)

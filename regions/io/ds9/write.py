# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy import coordinates
from ... import shapes

__all__ = [
    'write_ds9',
    'ds9_objects_to_string',
]

coordsys_name_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                 coordinates.frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic'  # needs expert attention TODO


def ds9_objects_to_string(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    """Convert list of regions to ds9 region strings.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects

    Returns
    -------
    region_string : str
        ds9 region string
    fmt : str
        A python string format defining the output precision.  Default is .6f,
        which is accurate to 0.0036 arcseconds.


    Examples
    --------
    TODO
    """

    ids = {
        shapes.circle.CircleSkyRegion: 'skycircle',
        shapes.circle.CirclePixelRegion: 'pixcircle',
        shapes.annulus.CircleAnnulusSkyRegion: 'skycircle',
        shapes.annulus.CircleAnnulusPixelRegion: 'pixcircle',
        shapes.ellipse.EllipseSkyRegion: 'skyellipse',
        shapes.ellipse.EllipsePixelRegion: 'pixellipse',
        shapes.rectangle.RectangleSkyRegion: 'skyrectangle',
        shapes.rectangle.RectanglePixelRegion: 'pixrectangle',
        shapes.polygon.PolygonSkyRegion: 'skypolygon',
        shapes.polygon.PolygonPixelRegion: 'pixpolygon',
        shapes.line.LineSkyRegion: 'skyline',
        shapes.line.LinePixelRegion: 'pixline',
    }

    if radunit == 'arcsec':
        if coordsys in coordsys_name_mapping.keys():
            radunitstr = '"'
        else:
            raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
    else:
        radunitstr = ''

    ds9_strings = {
        'circle': 'circle({x:' + fmt + '},{y:' + fmt + '},{r:' + fmt + '}' + radunitstr + ')',
        'annulus': 'annulus({x:' + fmt + '},{y:' + fmt + '},{inner:' + fmt + '}' + radunitstr + ',{outer:' + fmt + '}' + radunitstr + ')',
        'ellipse': 'ellipse({x:' + fmt + '},{y:' + fmt + '},{r1:' + fmt + '}' + radunitstr + ',{r2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})',
        'rectangle': 'box({x:' + fmt + '},{y:' + fmt + '},{h1:' + fmt + '}' + radunitstr + ',{h2:' + fmt + '}' + radunitstr + ',{ang:' + fmt + '})',
        'polygon': 'polygon({c})',
        'point': 'point({{x:{fmt}}}, {{y:{fmt}}})'.format(fmt=fmt),
        'line': 'line({{x1:{fmt}}}, {{y1:{fmt}}}, {{x2:{fmt}}}, {{y2:{fmt}}}, )'.format(fmt=fmt),
    }
    for key in ds9_strings:
        # include must be a prefix "-" or ""
        ds9_strings[key] = "{include}" + ds9_strings[key]

    output = '# Region file format: DS9 astropy/regions\n'
    output += '{}\n'.format(coordsys)

    # convert coordsys string to coordsys object
    if coordsys in coordsys_name_mapping:
        frame = coordinates.frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
    else:
        # for pixel/image/physical frames
        frame = None

    for reg in regions:

        # default if unspecified is that include is True, which means we
        # prepend nothing
        include = reg.meta.get('include', '')

        meta_str = " ".join("{0}={{{1}}}".format(key,val) for key,val in
                            reg.meta.items() if key not in ('include','tag'))
        if 'tag' in reg.meta:
            meta_str += " "+" ".join(["tag={{{0}}}".format(tag)
                                      for tag in reg.meta['tag']])
        if 'comment' in reg.meta:
            meta_str += " " + reg.meta['comment']

        if isinstance(reg, shapes.circle.CircleSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            r = float(reg.radius.to(radunit).value)
            line = ds9_strings['circle'].format(**locals())
        elif isinstance(reg, shapes.circle.CirclePixelRegion):
            x = reg.center.x
            y = reg.center.y
            r = reg.radius
            line = ds9_strings['circle'].format(**locals())
        elif isinstance(reg, shapes.annulus.CircleAnnulusSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            inner = float(reg.inner_radius.to(radunit).value)
            outer = float(reg.outer_radius.to(radunit).value)
            line = ds9_strings['annulus'].format(**locals())
        elif isinstance(reg, shapes.annulus.CircleAnnulusPixelRegion):
            x = reg.center.x
            y = reg.center.y
            inner = reg.inner_radius
            outer = reg.outer_radius
            line = ds9_strings['annulus'].format(**locals())
        elif isinstance(reg, shapes.point.PointSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            line = ds9_strings['point'].format(**locals())
        elif isinstance(reg, shapes.point.PointPixelRegion):
            x = reg.center.x
            y = reg.center.y
            line = ds9_strings['point'].format(**locals())
        elif isinstance(reg, shapes.line.LineSkyRegion):
            x1 = float(reg.start.transform_to(frame).spherical.lon.to('deg').value)
            y1 = float(reg.start.transform_to(frame).spherical.lat.to('deg').value)
            x2 = float(reg.end.transform_to(frame).spherical.lon.to('deg').value)
            y2 = float(reg.end.transform_to(frame).spherical.lat.to('deg').value)
            line = ds9_strings['line'].format(**locals())
        elif isinstance(reg, shapes.line.LinePixelRegion):
            x1 = reg.start.x
            y1 = reg.start.y
            x2 = reg.end.x
            y2 = reg.end.y
            line = ds9_strings['line'].format(**locals())
        elif isinstance(reg, shapes.ellipse.EllipseSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            r1 = float(reg.width.to(radunit).value)
            r2 = float(reg.height.to(radunit).value)
            ang = float(reg.angle.to('deg').value)
            line = ds9_strings['ellipse'].format(**locals())
        elif isinstance(reg, shapes.rectangle.RectangleSkyRegion):
            x = float(reg.center.transform_to(frame).spherical.lon.to('deg').value)
            y = float(reg.center.transform_to(frame).spherical.lat.to('deg').value)
            h1 = float(reg.width.to(radunit).value)
            h2 = float(reg.height.to(radunit).value)
            ang = float(reg.angle.to('deg').value)
            line = ds9_strings['rectangle'].format(**locals())
        elif isinstance(reg, shapes.ellipse.EllipsePixelRegion):
            x = reg.center.x
            y = reg.center.y
            r1 = reg.width
            r2 = reg.height
            ang = reg.angle
            line = ds9_strings['ellipse'].format(**locals())
        elif isinstance(reg, shapes.rectangle.RectanglePixelRegion):
            x = reg.center.x
            y = reg.center.y
            h1 = reg.width
            h2 = reg.height
            ang = reg.angle
            line = ds9_strings['rectangle'].format(**locals())
        elif isinstance(reg, shapes.polygon.PolygonSkyRegion):
            v = reg.vertices.transform_to(frame)
            coords = [(x.to('deg').value, y.to('deg').value) for x, y in
                      zip(v.spherical.lon, v.spherical.lat)]
            val = "{:" + fmt + "}"
            temp = [val.format(x) for _ in coords for x in _]
            c = ",".join(temp)
            line = ds9_strings['polygon'].format(**locals())
        elif isinstance(reg, shapes.polygon.PolygonPixelRegion):
            v = reg.vertices
            coords = [(x, y) for x, y in zip(v.x, v.y)]
            val = "{:" + fmt + "}"
            temp = [val.format(x) for _ in coords for x in _]
            c = ",".join(temp)
            line = ds9_strings['polygon'].format(**locals())
        else:
            raise ValueError("Cannot write {0}".format(reg))

        output += "{0} # {1} \n".format(line, meta_str)

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

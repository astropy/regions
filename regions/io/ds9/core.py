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
from ... import shapes
from ...core import PixCoord

__all__ = [
    'Shape',
    'ShapeList',
]

class ShapeList(list):
    """List of Shape
    """
    def to_region(self):
        regions = list()
        for shape in self:
            # Skip elliptical annuli for now
            if shape.region_type == 'ellipse' and len(shape.coord) > 5:
                from .read import DS9RegionParserWarning
                msg = 'Skipping elliptical annulus {}'.format(shape)
                warn(msg, DS9RegionParserWarning)
                continue
            log.debug(shape)
            region = shape.to_region()
            log.debug(region)
            regions.append(region)
        return regions

class Shape(object):
    """Helper class to represent a DS9
    
    This serves as intermediate step in the parsing process.

    Parameters
    ---------
    coordsys : str
        Coordinate system
    region_type : str
        Region type
    coord : list of `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        Coordinates
    meta : dict
        Meta attributed
    composite : bool
        Composite region
    include : bool
        Include/exclude region
    """
    shape_to_sky_region = dict(
        circle=shapes.CircleSkyRegion,
        ellipse=shapes.EllipseSkyRegion,
        box=shapes.RectangleSkyRegion,
        polygon=shapes.PolygonSkyRegion,
        annulus=shapes.CircleAnnulusSkyRegion,
        line=shapes.LineSkyRegion,
        point=shapes.PointSkyRegion,
    )
    shape_to_pixel_region = dict(
        circle=shapes.CirclePixelRegion,
        ellipse=shapes.EllipsePixelRegion,
        box=shapes.RectanglePixelRegion,
        polygon=shapes.PolygonPixelRegion,
        annulus=shapes.CircleAnnulusPixelRegion,
        line=shapes.LinePixelRegion,
        point=shapes.PointPixelRegion,
    )

    def __init__(self, coordsys, region_type, coord, meta, composite, include):
        self.coordsys = coordsys
        self.region_type=region_type
        self.coord=coord
        self.meta=meta
        self.composite=composite
        self.include=include

    def __str__(self):
        ss = self.__class__.__name__
        ss += '\nCoord sys : {}'.format(self.coordsys)
        ss += '\nRegion type : {}'.format(self.region_type)
        ss += '\nCoord: {}'.format(self.coord)
        ss += '\nMeta: {}'.format(self.meta)
        ss += '\nComposite: {}'.format(self.composite)
        ss += '\nInclude: {}'.format(self.include)
        ss += '\n'
        return ss

    def convert_coords(self):
        """Process list of coordinates

        This mainly seaches for tuple of coordinates in the coordinate list and
        creates a SkyCoord or PixCoord object from them if appropriate for a
        given region type. This involves again some coordinate transformation,
        so this step could be moved to the parsing process
        """
        from .read import DS9Parser
        if self.coordsys in DS9Parser.coordsys_mapping:
            coords = self._convert_sky_coords()
        else:
            coords = self._convert_pix_coords()

        if self.region_type == 'line':
            coords = [coords[0][0], coords[0][1]]

        return coords

    def _convert_sky_coords(self):
        """Convert to sky coords"""
        parsed_angles = [(x, y)
                         for x, y in zip(self.coord[:-1:2], self.coord[1::2])
                         if (isinstance(x, coordinates.Angle) and
                             isinstance(y, coordinates.Angle))
                        ]
        frame = coordinates.frame_transform_graph.lookup_name(self.coordsys)

        lon, lat = zip(*parsed_angles)
        if hasattr(lon, '__len__') and hasattr(lon, '__lat__') and len(lon) == 1 and len(lat==1):
            # force entries to be scalar if they are length-1
            lon, lat = u.Quantity(lon[0]), u.Quantity(lat[0])
        else:
            # otherwise, they are vector quantitites
            lon, lat = u.Quantity(lon), u.Quantity(lat)
        sphcoords = coordinates.UnitSphericalRepresentation(lon, lat)
        coords = [frame(sphcoords)]

        if self.region_type != 'polygon':
            coords += self.coord[len(coords * 2):]

        return coords

    def _convert_pix_coords(self):
        """Convert to pix coords"""
        if self.region_type in ['polygon', 'line']:
            # have to special-case polygon in the phys coord case
            # b/c can't typecheck when iterating as in sky coord case
            coords = [PixCoord(self.coord[0::2], self.coord[1::2])]
        else:
            temp = [_.value for _ in self.coord]
            coord = PixCoord(temp[0], temp[1])
            coords = [coord] + temp[2:]

        return coords

    def to_region(self):
        """Convert to region object
        """
        coords = self.convert_coords()
        log.debug(coords)
        viz_keywords = ['color', 'dashed', 'width', 'point', 'font']

        if isinstance(coords[0], BaseCoordinateFrame):
            reg = self.shape_to_sky_region[self.region_type](*coords)
        elif isinstance(coords[0], PixCoord):
            reg = self.shape_to_pixel_region[self.region_type](*coords)
        else:
            raise DS9RegionParserError("No central coordinate")

        reg.visual = {key: self.meta[key] for key in self.meta.keys() if key in viz_keywords}
        reg.meta = {key: self.meta[key] for key in self.meta.keys() if key not in viz_keywords}
        return reg

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function
from warnings import warn
import string

from astropy import units as u
from astropy import coordinates
from astropy.coordinates import BaseCoordinateFrame, Angle
from astropy import log

from .. import shapes
from ..core import PixCoord, SkyRegion, RegionMeta, RegionVisual
from .ds9.core import DS9RegionParserWarning, DS9RegionParserError
from .crtf.core import CRTFRegionParserWarning, CRTFRegionParserError

__all__ = ['ShapeList', 'Shape', 'to_shape_list', 'to_crtf_meta', 'to_ds9_meta']

regions_attributes = dict(circle=['center', 'radius'],
                          ellipse=['center', 'width', 'height', 'angle'],
                          rectangle=['center', 'width', 'height', 'angle'],
                          polygon=['vertices'],
                          annulus=['center', 'inner_radius', 'outer_radius'],
                          line=['start', 'end'],
                          point=['center']
                          )

# This helps to map the region names in the respective format to the ones available in this package
reg_mapping = {'DS9': {x: x for x in regions_attributes},
               'CRTF': {x: x for x in regions_attributes}}
reg_mapping['DS9']['box'] = 'rectangle'
reg_mapping['CRTF']['rotbox'] = 'rectangle'
reg_mapping['CRTF']['poly'] = 'polygon'


# valid astropy coordinate frames in their respective formats.
valid_coordsys = {'DS9': ['image', 'physical', 'fk4', 'fk5', 'icrs', 'galactic',
                          'geocentrictrueecliptic', 'wcs'],
                  'CRTF': ['image', 'fk5', 'fk4', 'galactic',
                           'geocentrictrueecliptic', 'supergalactic', 'icrs']
                  }
valid_coordsys['DS9'] += ['wcs{}'.format(x) for x in string.ascii_lowercase]

# Maps astropy's coordinate frame names with their respective name in the file format.
coordsys_mapping = {'DS9': {x: x for x in valid_coordsys['DS9']},
                    'CRTF': {x: x for x in valid_coordsys['CRTF']}
                    }
coordsys_mapping['CRTF']['geocentrictrueecliptic'] = 'ecliptic'
coordsys_mapping['CRTF']['fk5'] = 'j2000'
coordsys_mapping['CRTF']['fk4'] = 'b1950'
coordsys_mapping['CRTF']['supergalactic'] = 'supergal'

coordsys_mapping['DS9']['geocentrictrueecliptic'] = 'ecliptic'


class ShapeList(list):
    """
    List of Shape
    """
    def to_regions(self):
        regions = list()
        for shape in self:
            # Skip elliptical annulus for now
            if shape.region_type == 'ellipse' and len(shape.coord) > 5:
                msg = 'Skipping elliptical annulus {}'.format(shape)
                warn(msg, DS9RegionParserWarning)
                continue
            if shape.region_type in ['box', 'symbol', 'text'] and shape.format_type == 'CRTF':
                msg = 'Skipping {} {}'.format(shape.region_type, shape)
                warn(msg, CRTFRegionParserWarning)
                continue
            log.debug(shape)
            region = shape.to_region()
            log.debug(region)
            regions.append(region)
        return regions

    def to_crtf(self):
        raise NotImplementedError('This method intends to convert ShapeList to CRTF strings')

    def to_ds9(self, coordsys='fk5', fmt='.6f', radunit='deg'):
        """
        Converts a list of ``regions.Shape`` objects to ds9 region strings.

        Parameters
        ----------
        coordsys : str
            This overrides the coordinate system frame for all regions.

        fmt : str
            A python string format defining the output precision.
            Default is .6f, which is accurate to 0.0036 arcseconds.

        radunit : str
            This denotes the unit of the radius.

        Returns
        -------
        region_string : str
            ds9 region string

        Examples
        --------
        TODO
        """

        ds9_strings = {
            'circle': '{0}circle({1:FMT},{2:FMT},{3:FMT}RAD)',
            'annulus': '{0}annulus({1:FMT},{2:FMT},{3:FMT}RAD,{4:FMT}RAD)',
            'ellipse': '{0}ellipse({1:FMT},{2:FMT},{3:FMT}RAD,{4:FMT}RAD,{5:FMT})',
            'rectangle': '{0}box({1:FMT},{2:FMT},{3:FMT}RAD,{4:FMT}RAD,{5:FMT})',
            'polygon': '{0}polygon({1})',
            'point': '{0}point({1:FMT},{2:FMT})',
            'line': '{0}line({1:FMT},{2:FMT},{3:FMT},{4:FMT})'
                      }

        output = '# Region file format: DS9 astropy/regions\n'

        if radunit == 'arcsec':
            # what's this for?
            if coordsys in coordsys_mapping['DS9'].values():
                radunitstr = '"'
            else:
                raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else:
            radunitstr = ''

        for key, val in ds9_strings.items():
            ds9_strings[key] = val.replace("FMT", fmt).replace("RAD", radunitstr)

        output += '{}\n'.format(coordsys)

        for shape in self:

            shape.check_ds9()

            # if unspecified, include is True.
            include = "-" if shape.meta.get('include') in (False, '-') else ""

            meta_str = " ".join("{0}={1}".format(key, val) for key, val in
                                shape.meta.items() if key not in ('include', 'tag', 'comment'))
            if 'tag' in shape.meta:
                meta_str += " " + " ".join(["tag={0}".format(tag) for tag in shape.meta['tag']])
            if 'comment' in shape.meta:
                meta_str += " " + shape.meta['comment']

            coord = []

            if coordsys not in ['image', 'physical']:
                for val in shape.coord:
                    if isinstance(val, Angle):
                        coord.append(float(val.value))
                    else:
                        if radunit == '' or None:
                            coord.append(float(val.value))
                        else:
                            coord.append(float(val.to(radunit).value))
                if shape.region_type in ['ellipse', 'rectangle'] and len(shape.coord) % 2 == 1:
                    coord[-1] = float(shape.coord[-1].to('deg').value)
            else:
                for val in shape.coord:
                    if isinstance(val, u.Quantity):
                        coord.append(float(val.value))
                    else:
                        coord.append(float(val))
                if shape.region_type in ['polygon', 'line']:
                    coord = [x+1 for x in coord]
                else:
                    coord[0] += 1
                    coord[1] += 1

            if shape.region_type == 'polygon':
                val = "{0:" + fmt + "}"
                temp = [val.format(x) for x in coord]
                coord = ",".join(temp)
                line = ds9_strings['polygon'].format(include, coord)
            else:
                line = ds9_strings[shape.region_type].format(include, *coord)

            if meta_str.strip():
                output += "{0} # {1} \n".format(line, meta_str)
            else:
                output += "{0}\n".format(line)

        return output


class Shape(object):
    """
    Helper class to represent a DS9/CRTF Region.

    This serves as intermediate step in the parsing process.

    Parameters
    ----------
    format_type : str
        File Format type

    coordsys : str
        Astropy Coordinate system frame used in the region.

    region_type : str
        Type of the region (as defined in this package).

    coord : list of `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        Coordinates

    meta : dict
        Meta attributes

    composite : bool
        Composite region

    include : bool
        Include/exclude region
    """

    shape_to_sky_region = dict(circle=shapes.CircleSkyRegion,
                               ellipse=shapes.EllipseSkyRegion,
                               rectangle=shapes.RectangleSkyRegion,
                               polygon=shapes.PolygonSkyRegion,
                               annulus=shapes.CircleAnnulusSkyRegion,
                               line=shapes.LineSkyRegion,
                               point=shapes.PointSkyRegion
                               )

    shape_to_pixel_region = dict(circle=shapes.CirclePixelRegion,
                                 ellipse=shapes.EllipsePixelRegion,
                                 rectangle=shapes.RectanglePixelRegion,
                                 polygon=shapes.PolygonPixelRegion,
                                 annulus=shapes.CircleAnnulusPixelRegion,
                                 line=shapes.LinePixelRegion,
                                 point=shapes.PointPixelRegion
                                 )

    error = {'DS9': DS9RegionParserError, 'CRTF': CRTFRegionParserError}
    warning = {'DS9': DS9RegionParserWarning, 'CRTF': CRTFRegionParserWarning}

    def __init__(self, format_type, coordsys, region_type, coord, meta, composite, include):

        from . import CRTFRegionParser, DS9Parser
        self.parser = {'DS9': DS9Parser, 'CRTF': CRTFRegionParser}

        self._format_type = format_type
        self._coordsys = coordsys
        self._region_type = region_type
        self.coord = coord
        self.meta = meta
        self.composite = composite
        self.include = include
        self._validate(self.format_type)

    @property
    def format_type(self):
        return self._format_type

    @property
    def coordsys(self):
        return self._coordsys

    @coordsys.setter
    def coordsys(self, value):
        self._coordsys = value.lower()
        self._validate(self.format_type)

    @property
    def region_type(self):
        return self._region_type

    @region_type.setter
    def region_type(self, value):
        self._region_type = value.lower()
        self._validate(self.format_type)

    def __str__(self):
        ss = self.__class__.__name__
        ss += '\nFormat Type : {}'.format(self.format_type)
        if self.format_type == 'CRTF':
            ss += '\nType : {}'.format(self.meta.get('type', 'reg'))
        ss += '\nCoord sys : {}'.format(self.coordsys)
        ss += '\nRegion type : {}'.format(self.region_type)
        if self.region_type == 'symbol':
            ss += '\nSymbol : {}'.format(self.meta['symbol'])
        if self.region_type == 'text':
            ss += '\nText : {}'.format(self.meta['string'])
        ss += '\nCoord: {}'.format(self.coord)
        ss += '\nMeta: {}'.format(self.meta)
        ss += '\nComposite: {}'.format(self.composite)
        ss += '\nInclude: {}'.format(self.include)
        ss += '\n'
        return ss

    def convert_coords(self):
        """
        Process list of coordinates

        This mainly searches for tuple of coordinates in the coordinate list and
        creates a SkyCoord or PixCoord object from them if appropriate for a
        given region type. This involves again some coordinate transformation,
        so this step could be moved to the parsing process
        """
        if self.coordsys in self.parser[self.format_type].coordsys_mapping:
            coords = self._convert_sky_coords()
        else:
            coords = self._convert_pix_coords()

        if self.region_type == 'line':
            coords = [coords[0][0], coords[0][1]]

        return coords

    def _convert_sky_coords(self):
        """
        Convert to sky coordinates
        """
        parsed_angles = [(x, y)
                         for x, y in zip(self.coord[:-1:2], self.coord[1::2])
                         if (isinstance(x, coordinates.Angle) and isinstance(y, coordinates.Angle))
                         ]
        frame = coordinates.frame_transform_graph.lookup_name(self.coordsys)

        lon, lat = zip(*parsed_angles)
        if hasattr(lon, '__len__') and hasattr(lat, '__len__') and len(lon) == 1 and len(lat) == 1:
            # force entries to be scalar if they are length-1
            lon, lat = u.Quantity(lon[0]), u.Quantity(lat[0])
        else:
            # otherwise, they are vector quantities
            lon, lat = u.Quantity(lon), u.Quantity(lat)
        sphcoords = coordinates.UnitSphericalRepresentation(lon, lat)
        coords = [frame(sphcoords)]

        if self.region_type != 'polygon':
            coords += self.coord[len(coords * 2):]

        return coords

    def _convert_pix_coords(self):
        """
        Convert to pixel coordinates, `regions.PixCoord`
        """
        if self.region_type in ['polygon', 'line', 'poly']:
            # have to special-case polygon in the phys coord case
            # b/c can't typecheck when iterating as in sky coord case
            coords = [PixCoord(self.coord[0::2], self.coord[1::2])]
        else:
            temp = [_.value for _ in self.coord]
            coord = PixCoord(temp[0], temp[1])
            coords = [coord] + temp[2:]

        return coords

    def to_region(self):
        """
        Converts to region, ``regions.Region`` object
        """

        coords = self.convert_coords()
        log.debug(coords)
        viz_keywords = ['color', 'dash', 'dashlist', 'width', 'font', 'symsize',
                        'symsize', 'fontsize', 'fontstyle', 'usetex',
                        'labelpos', 'labeloff', 'linewidth', 'linestyle']

        if isinstance(coords[0], BaseCoordinateFrame):
            reg = self.shape_to_sky_region[self.region_type](*coords)
        elif isinstance(coords[0], PixCoord):
            reg = self.shape_to_pixel_region[self.region_type](*coords)
        else:
            self._raise_error("No central coordinate")

        reg.visual = RegionVisual()
        reg.meta = RegionMeta()
        for key in self.meta.keys():
            if key in viz_keywords:
                reg.visual[key] = self.meta[key]
            else:
                reg.meta[key] = self.meta[key]
        reg.meta['include'] = self.include
        return reg

    def _raise_error(self, msg):
        raise self.error[self.format_type](msg)

    def check_crtf(self):
        """
        Checks for CRTF compatibility.
        """
        self._validate('CRTF')

    def check_ds9(self):
        """
        Checks for DS9 compatibility.
        """
        self._validate('DS9')

    def _validate(self, format_type):
        """
        Checks whether all the attributes of this object is valid according to the given format.
        """
        if format_type not in ['CRTF', 'DS9']:
            raise ValueError("'{0}' is not available as io".format(format_type))
        if self.region_type not in regions_attributes.keys():
            raise ValueError("'{0}' is not a valid region type in this package".format(self.region_type, format_type))
        if self.coordsys not in valid_coordsys[format_type]:
            raise ValueError("'{0}' is not a valid coordinate reference frame in astropy".
                             format(self.coordsys, format_type))


def to_shape_list(region_list, format_type='DS9', coordinate_system='fk5'):
    """
    Converts a list of regions into a `regions.ShapeList` object.

    Parameters
    ----------
    region_list: python list
        Lists of `regions.Region` objects

    format_type: str ('DS9' or 'CRTF')
        The format type of the Shape object. Default is 'DS9'.

    coordinate_system: str
        The astropy coordinate system frame in which all the coordinates present
        in the `region_list` will be converted. Default is 'fk5'.

    Returns
    -------
    shape_list: `regions.ShapeList` object
        list of `regions.Shape` objects.
    """

    shape_list = ShapeList()

    for region in region_list:

        coord = []
        reg_type = str(type(region)).split(".")[2]

        for val in regions_attributes[reg_type]:
            coord.append(getattr(region, val))

        if reg_type == 'polygon':
            coord = [x for x in region.vertices]

        if coordinate_system:
            coordsys = coordinate_system
        else:
            if isinstance(region, SkyRegion):
                coordsys = coord[0].name
            else:
                coordsys = 'image'

        frame = coordinates.frame_transform_graph.lookup_name(coordsys)

        new_coord = []
        for val in coord:
            if isinstance(val, Angle) or isinstance(val, u.Quantity) or isinstance(val, float):
                new_coord.append(val)
            elif isinstance(val, PixCoord):
                new_coord.append(u.Quantity(val.x, u.dimensionless_unscaled))
                new_coord.append(u.Quantity(val.y, u.dimensionless_unscaled))
            else:
                new_coord.append(Angle(val.transform_to(frame).spherical.lon))
                new_coord.append(Angle(val.transform_to(frame).spherical.lat))

        if format_type == 'DS9':
            meta = to_ds9_meta(region.meta, region.visual)
        else:
            meta = to_crtf_meta(region.meta, region.visual)

        shape_list.append(Shape(format_type, coordsys, reg_type, new_coord, meta, False,
                                region.meta.get('include', False)))

    return shape_list


def to_ds9_meta(region_meta, region_visual):
    """
    Makes the meta data DS9 compatible by filtering and mapping the valid keys

    Parameters
    ----------
    region_meta: dict
        meta attribute of a `regions.Region` object

    region_visual : dict
        visual attribute of a `regions.Region` object.

    Returns
    -------
    meta : dict
        DS9 compatible meta dictionary

    """

    # meta keys allowed in DS9.
    valid_keys = ['label', 'symbol', 'include', 'tag', 'line', 'comment',
                  'name', 'select', 'highlite', 'fixed',
                  'edit', 'move', 'rotate', 'delete', 'source', 'background']

    # visual keys allowed in DS9
    valid_keys += ['color', 'dash', 'linewidth', 'font', 'dashlist', 'fill']

    # mapped to actual names in DS9
    key_mappings = {'symbol': 'point', 'label': 'text', 'linewidth': 'width'}

    return _to_io_meta(region_meta, region_visual, valid_keys, key_mappings)


def to_crtf_meta(region_meta, region_visual):
    """
    Makes the meta data CRTF compatible by filtering and mapping the valid keys

    Parameters
    ----------
    region_meta: dict
        meta attribute of a `regions.Region` object

    region_visual : dict
        visual attribute of a `regions.Region` object

    Returns
    -------
    meta : dict
        CRTF compatible meta dictionary

    """
    # please refer : https://casaguides.nrao.edu/index.php/CASA_Region_Format

    # meta keys allowed in CRTF
    valid_keys = ['label', 'symbol', 'include', 'frame', 'range', 'veltype',
                  'restfreq', 'coord']

    # visual keys allowed in CRTF
    valid_keys += ['color', 'width', 'font', 'symthick', 'symsize', 'fontsize',
                   'fontstyle', 'usetex', 'labelpos', 'labeloff', 'linewidth',
                   'linestyle']

    key_mappings = {}

    return _to_io_meta(region_meta, region_visual, valid_keys, key_mappings)


def _to_io_meta(region_meta, region_visual, valid_keys, key_mappings):
    """
    This is used to make meta data compatible with a specific io
    by filtering and mapping to it's valid keys

    Parameters
    ----------
    region_meta: dict
        meta attribute of a `regions.Region` object

    region_visual : dict
        visual attribute of a `regions.Region` object.

    valid_keys : python list
        Contains all the valid keys of a particular file format.

    key_mappings : python dict
        Maps to the actual name of the key in the format.

    Returns
    -------
    meta : dict
        io compatible meta dictionary according to valid_keys and key_mappings

    """

    meta = dict()

    for key in region_meta:
        if key in valid_keys:
            meta[key_mappings.get(key, key)] = region_meta[key]
    for key in region_visual:
        if key in valid_keys:
            meta[key_mappings.get(key, key)] = region_visual[key]

    return meta


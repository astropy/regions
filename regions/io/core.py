# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function
from warnings import warn
import string
import numbers
import numpy as np

from astropy import units as u
from astropy import coordinates
from astropy.coordinates import Angle, SkyCoord
from astropy import log
from astropy.table import Table

from .. import shapes
from ..core import PixCoord, SkyRegion
from ..core.attributes import RegionMeta, RegionVisual
from .ds9.core import DS9RegionParserWarning, valid_symbols_ds9
from .crtf.core import CRTFRegionParserWarning

__all__ = ['ShapeList', 'Shape', 'to_shape_list', 'to_crtf_meta', 'to_ds9_meta']


class RegionConversionError(ValueError):
    """
    A generic error class for Shape to Region conversion.
    """


regions_attributes = dict(circle=['center', 'radius'],
                          ellipse=['center', 'width', 'height', 'angle'],
                          rectangle=['center', 'width', 'height', 'angle'],
                          polygon=['vertices'],
                          circleannulus=['center', 'inner_radius', 'outer_radius'],
                          ellipseannulus=['center', 'inner_width',
                                          'inner_height', 'outer_width',
                                          'outer_height', 'angle'],
                          line=['start', 'end'],
                          point=['center'],
                          text=['center']
                          )
regions_attributes['rectangleannulus'] = regions_attributes['ellipseannulus']

# This helps to map the region names in the respective format to the ones
# available in this package
reg_mapping = {'DS9': {x: x for x in regions_attributes},
               'CRTF': {x: x for x in regions_attributes},
               'FITS_REGION': {x: x for x in regions_attributes}}
reg_mapping['DS9']['box'] = 'rectangle'
reg_mapping['CRTF']['rotbox'] = 'rectangle'
reg_mapping['CRTF']['box'] = 'rectangle'
reg_mapping['CRTF']['centerbox'] = 'rectangle'
reg_mapping['CRTF']['poly'] = 'polygon'
reg_mapping['CRTF']['symbol'] = 'point'
reg_mapping['CRTF']['text'] = 'text'
reg_mapping['CRTF']['annulus'] = 'circleannulus'
reg_mapping['DS9']['text'] = 'text'
reg_mapping['DS9']['annulus'] = 'circleannulus'
reg_mapping['FITS_REGION']['annulus'] = 'circleannulus'
reg_mapping['FITS_REGION']['box'] = 'rectangle'
reg_mapping['FITS_REGION']['rotbox'] = 'rectangle'
reg_mapping['FITS_REGION']['elliptannulus'] = 'ellipseannulus'

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
    """List of Shape"""
    def to_regions(self):
        regions = list()
        for shape in self:
            # Skip elliptical annulus for now
            if shape.region_type == 'ellipse' and len(shape.coord) > 5:
                msg = 'Skipping elliptical annulus {}'.format(shape)
                warn(msg, DS9RegionParserWarning)
                continue
            if shape.region_type in ['box'] and shape.format_type == 'CRTF':
                msg = 'Skipping {} {}'.format(shape.region_type, shape)
                warn(msg, CRTFRegionParserWarning)
                continue
            log.debug(shape)
            region = shape.to_region()
            log.debug(region)
            regions.append(region)
        return regions

    def to_crtf(self,  coordsys='fk5', fmt='.6f', radunit='deg'):
        """
        Converts a list of ``regions.Shape`` objects to crtf region strings.

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
            crtf region string

        Examples
        --------
        TODO
        """

        crtf_strings = {
            'circle': '{0}circle[[{1:FMT}deg, {2:FMT}deg], {3:FMT}RAD]',
            'circleannulus': '{0}annulus[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}RAD, {4:FMT}RAD]]',
            'ellipse': '{0}ellipse[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}RAD, {4:FMT}RAD], {5:FMT}deg]',
            'rectangle': '{0}rotbox[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}RAD, {4:FMT}RAD], {5:FMT}deg]',
            'polygon': '{0}poly[{1}]',
            'point': '{0}point[[{1:FMT}deg, {2:FMT}deg]]',
            'symbol': '{0}symbol[[{1:FMT}deg, {2:FMT}deg], {symbol}]',
            'text': '{0}text[[{1:FMT}deg, {2:FMT}deg], \'{text}\']',
            'line': '{0}line[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}deg, {4:FMT}deg]]'
                        }

        output = '#CRTF\n'

        if radunit == 'arcsec':
            # what's this for?
            if coordsys in coordsys_mapping['CRTF'].values():
                radunitstr = '"'
            else:
                raise ValueError(
                    'Radius unit arcsec not valid for coordsys {}'.format(
                        coordsys))
        else:
            radunitstr = radunit

        for key, val in crtf_strings.items():
            crtf_strings[key] = val.replace("FMT", fmt).replace("RAD",
                                                                radunitstr)

        # CASA does not support global coordinate specification, even though the
        # documentation for the specification explicitly states that it does.
        # output += 'global coord={}\n'.format(coordsys)

        for shape in self:

            shape.check_crtf()
            shape.meta = to_crtf_meta(shape.meta)

            # if unspecified, include is True.
            # Despite the specification, CASA does *not* support a preceding
            # "+".  If you want a region included, leave the opening character
            # blank.
            include = "-" if shape.include in (False, '-') else ""
            include += "ann " if shape.meta.get('type', 'reg') == 'ann' else ""

            if shape.meta.get('label', "") != "":
                shape.meta['label'] = "'{}'".format(shape.meta['label'])
            meta_str = ", ".join("{0}={1}".format(key, val) for key, val in
                                shape.meta.items() if
                                key not in ('include', 'comment', 'symbol',
                                            'coord', 'text', 'range', 'corr',
                                            'type'))

            # the first item should be the coordinates, since CASA cannot
            # recognize a region without an inline coordinate specification
            # It can be, but does not need to be, comma-separated at the start
            meta_str = "coord={0}, ".format(coordsys.upper()) + meta_str

            if 'comment' in shape.meta:
                meta_str += ", " + shape.meta['comment']
            if 'range' in shape.meta:
                shape.meta['range'] = [str(str(x).replace(" ", "")) for x in
                                       shape.meta['range']]
                meta_str += ", range={}".format(shape.meta['range']).replace("'", "")
            if 'corr' in shape.meta:
                meta_str += ", corr={}".format(shape.meta['corr']).replace("'", "")

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
            else:
                for val in shape.coord:
                    if isinstance(val, u.Quantity):
                        coord.append(float(val.value))
                    else:
                        coord.append(float(val))

            if shape.region_type in ['ellipse', 'rectangle'] and len(shape.coord) % 2 == 1:
                coord[-1] = float(shape.coord[-1].to('deg').value)

            if shape.region_type == 'polygon':
                val = '[{0:' + fmt + '}deg, {1:' + fmt + '}deg]'
                temp = [val.format(x, y) for x, y in zip(coord[::2], coord[1::2])]
                coord = ", ".join(temp)
                line = crtf_strings['polygon'].format(include, coord)
            elif shape.region_type == 'point':
                if 'symbol' in shape.meta:
                    line = crtf_strings['symbol'].format(include, *coord,
                                                         symbol=shape.meta['symbol'])
                else:
                    line = crtf_strings['point'].format(include, *coord)
            elif shape.region_type == 'ellipse':
                coord[2:] = [x / 2 for x in coord[2:]]
                if len(coord) % 2 == 1:
                    coord[-1] *= 2
                line = crtf_strings['ellipse'].format(include, *coord)
            elif shape.region_type == 'text':
                line = crtf_strings['text'].format(include, *coord, text=shape.meta['text'])
            else:
                line = crtf_strings[shape.region_type].format(include, *coord)

            if meta_str.strip():
                output += "{0}, {1}\n".format(line, meta_str)
            else:
                output += "{0}\n".format(line)

        return output

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
        valid_symbols_reverse = {y: x for x, y in valid_symbols_ds9.items()}

        ds9_strings = {
            'circle': '{0}circle({1:FMT},{2:FMT},{3:FMT}RAD)',
            'circleannulus': '{0}annulus({1:FMT},{2:FMT},{3:FMT}RAD,{4:FMT}RAD)',
            'ellipse': '{0}ellipse({1:FMT},{2:FMT},{3:FMT}RAD,{4:FMT}RAD,{5:FMT})',
            'rectangle': '{0}box({1:FMT},{2:FMT},{3:FMT}RAD,{4:FMT}RAD,{5:FMT})',
            'polygon': '{0}polygon({1})',
            'point': '{0}point({1:FMT},{2:FMT})',
            'line': '{0}line({1:FMT},{2:FMT},{3:FMT},{4:FMT})',
            'text': '{0}text({1:FMT},{2:FMT})'
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
            shape.meta = to_ds9_meta(shape.meta)

            # if unspecified, include is True.
            include = "-" if shape.include in (False, '-') else ""

            if 'point' in shape.meta:
                shape.meta['point'] = valid_symbols_reverse[shape.meta['point']]

            if 'symsize' in shape.meta:
                shape.meta['point'] += " {}".format(shape.meta.pop('symsize'))

            meta_str = " ".join("{0}={1}".format(key, val) for key, val in
                                shape.meta.items() if key not in ('include', 'tag', 'comment', 'font', 'text'))
            if 'tag' in shape.meta:
                meta_str += " " + " ".join(["tag={0}".format(tag) for tag in shape.meta['tag']])
            if 'font' in shape.meta:
                meta_str += " " + 'font="{0}"'.format(shape.meta['font'])
            if shape.meta.get('text', '') != '':
                meta_str += " " + 'text={' + shape.meta['text'] + '}'
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
            elif shape.region_type == 'ellipse':
                coord[2:] = [x / 2 for x in coord[2:]]
                if len(coord) % 2 == 1:
                    coord[-1] *= 2
                line = ds9_strings['ellipse'].format(include, *coord)
            else:
                line = ds9_strings[shape.region_type].format(include, *coord)

            if meta_str.strip():
                output += "{0} # {1}\n".format(line, meta_str)
            else:
                output += "{0}\n".format(line)

        return output

    def to_fits(self):
        """
        Converts a `~regions.ShapeList` to a `~astropy.table.Table` object.
        """

        max_length_coord = 1
        coord_x = []
        coord_y = []
        shapes = []
        radius = []
        rotangle_deg = []
        components = []

        reg_reverse_mapping = {value: key for key, value in
                               reg_mapping['FITS_REGION'].items()}
        reg_reverse_mapping['rectangle'] = 'ROTBOX'
        reg_reverse_mapping['circleannulus'] = 'ANNULUS'
        reg_reverse_mapping['ellipseannulus'] = 'ELLIPTANNULUS'

        for num, shape in enumerate(self):
            shapes.append(reg_reverse_mapping[shape.region_type])
            if shape.region_type == 'polygon':
                max_length_coord = max(len(shape.coord)/2, max_length_coord)
                coord = [x.value for x in shape.coord]
                coord_x.append(coord[::2])
                coord_y.append(coord[1::2])
                radius.append(0)
                rotangle_deg.append(0)
            else:
                coord_x.append(shape.coord[0].value)
                coord_y.append(shape.coord[1].value)
                if shape.region_type in ['circle', 'circleannulus', 'point']:
                    radius.append([float(val) for val in shape.coord[2:]])
                    rotangle_deg.append(0)
                else:
                    radius.append([float(x) for x in shape.coord[2:-1]])
                    rotangle_deg.append(shape.coord[-1].to('deg').value)

            tag = shape.meta.get('tag', '')
            if tag.isdigit():
                components.append(int(tag))
            else:
                components.append(num + 1)

        # padding every value with zeros at the end to make sure that all values
        # in the column have same length.
        for i in range(len(self)):
            if np.isscalar(coord_x[i]):
                coord_x[i] = np.array([coord_x[i]])
            if np.isscalar(coord_y[i]):
                coord_y[i] = np.array([coord_y[i]])
            if np.isscalar(radius[i]):
                radius[i] = np.array([radius[i]])

            coord_x[i] = np.pad(coord_x[i], (0, int(max_length_coord - len(coord_x[i]))),
                                'constant', constant_values=(0, 0))
            coord_y[i] = np.pad(coord_y[i], (0, int(max_length_coord - len(coord_y[i]))),
                                'constant', constant_values=(0, 0))
            radius[i] = np.pad(radius[i], (0, 4 - len(radius[i])), 'constant',
                               constant_values=(0, 0))

        table = Table([coord_x, coord_y, shapes, radius, rotangle_deg, components],
                       names=('X', 'Y', 'SHAPE', 'R', 'ROTANG', 'COMPONENT'))
        table['X'].unit = 'pix'
        table['Y'].unit = 'pix'
        table['R'].unit = 'pix'
        table['ROTANG'].unit = 'deg'

        return table


class Shape(object):
    """
    Helper class to represent a DS9/CRTF Region.

    This serves as intermediate step in the parsing process.

    Parameters
    ----------
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
                               circleannulus=shapes.CircleAnnulusSkyRegion,
                               ellipseannulus=shapes.EllipseAnnulusSkyRegion,
                               rectangleannulus=shapes.RectangleAnnulusSkyRegion,
                               line=shapes.LineSkyRegion,
                               point=shapes.PointSkyRegion,
                               text=shapes.TextSkyRegion
                               )

    shape_to_pixel_region = dict(circle=shapes.CirclePixelRegion,
                                 ellipse=shapes.EllipsePixelRegion,
                                 rectangle=shapes.RectanglePixelRegion,
                                 polygon=shapes.PolygonPixelRegion,
                                 circleannulus=shapes.CircleAnnulusPixelRegion,
                                 ellipseannulus=shapes.EllipseAnnulusPixelRegion,
                                 rectangleannulus=shapes.RectangleAnnulusPixelRegion,
                                 line=shapes.LinePixelRegion,
                                 point=shapes.PointPixelRegion,
                                 text=shapes.TextPixelRegion
                                 )

    error = RegionConversionError

    def __init__(self, coordsys, region_type, coord, meta, composite, include):
        self._coordsys = coordsys
        self._region_type = region_type
        self.coord = coord
        self.meta = meta
        self.composite = composite
        self.include = include

    @property
    def coordsys(self):
        return self._coordsys

    @coordsys.setter
    def coordsys(self, value):
        self._coordsys = value.lower()
        self._validate()

    @property
    def region_type(self):
        return self._region_type

    @region_type.setter
    def region_type(self, value):
        self._region_type = value.lower()
        self._validate()

    def __str__(self):
        ss = self.__class__.__name__
        ss += '\nType : {}'.format(self.meta.get('type', 'reg'))
        ss += '\nCoord sys : {}'.format(self.coordsys)
        ss += '\nRegion type : {}'.format(self.region_type)
        if self.region_type == 'symbol':
            ss += '\nSymbol : {}'.format(self.meta['symbol'])
        if self.region_type == 'text':
            ss += '\nText : {}'.format(self.meta['text'])
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
        if self.coordsys in ['image', 'physical']:
            coords = self._convert_pix_coords()
        else:
            coords = self._convert_sky_coords()

        if self.region_type == 'line':
            coords = [coords[0][0], coords[0][1]]

        if self.region_type == 'text':
            coords.append(self.meta['text'])

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
        coords = [SkyCoord(frame(sphcoords))]

        if self.region_type != 'polygon':
            coords += self.coord[len(coords * 2):]

        return coords

    def _convert_pix_coords(self):
        """
        Convert to pixel coordinates, `regions.PixCoord`
        """
        if self.region_type in ['polygon', 'line']:
            # have to special-case polygon in the phys coord case
            # b/c can't typecheck when iterating as in sky coord case
            coords = [PixCoord(self.coord[0::2], self.coord[1::2])]
        else:
            temp = [_.value for _ in self.coord]
            coord = PixCoord(temp[0], temp[1])
            coords = [coord] + temp[2:]

        # The angle remains as a quantity object.
        # Modulus check makes sure that it works for ellipse/rectangle annulus
        if self.region_type in ['ellipse', 'rectangle'] and len(coords) % 2 == 0:
            coords[-1] = self.coord[-1]

        return coords

    def to_region(self):
        """
        Converts to region, ``regions.Region`` object
        """

        coords = self.convert_coords()
        log.debug(coords)
        viz_keywords = ['color', 'dash', 'dashlist', 'width', 'font', 'symsize',
                        'symbol', 'symsize', 'fontsize', 'fontstyle', 'usetex',
                        'labelpos', 'labeloff', 'linewidth', 'linestyle',
                        'point', 'textangle', 'fontweight']

        if isinstance(coords[0], SkyCoord):
            reg = self.shape_to_sky_region[self.region_type](*coords)
        elif isinstance(coords[0], PixCoord):
            reg = self.shape_to_pixel_region[self.region_type](*coords)
        else:
            self._raise_error("No central coordinate")

        reg.visual = RegionVisual()
        reg.meta = RegionMeta()

        # both 'text' and 'label' should be set to the same value, where we
        # default to the 'text' value since that is the one used by ds9 regions
        label = self.meta.get('text',
                              self.meta.get('label', ""))
        if label != '':
            reg.meta['label'] = label
        for key in self.meta:
            if key in viz_keywords:
                reg.visual[key] = self.meta[key]
            else:
                reg.meta[key] = self.meta[key]
        reg.meta['include'] = self.include
        return reg

    def _raise_error(self, msg):
        raise self.error(msg)

    def check_crtf(self):
        """
        Checks for CRTF compatibility.
        """
        if self.region_type not in regions_attributes:
            raise ValueError("'{0}' is not a valid region type in this package"
                             "supported by CRTF".format(self.region_type))

        if self.coordsys not in valid_coordsys['CRTF']:
            raise ValueError("'{0}' is not a valid coordinate reference frame in "
                             "astropy supported by CRTF".format(self.coordsys))

    def check_ds9(self):
        """
        Checks for DS9 compatibility.
        """
        if self.region_type not in regions_attributes:
            raise ValueError("'{0}' is not a valid region type in this package"
                             "supported by DS9".format(self.region_type))

        if self.coordsys not in valid_coordsys['DS9']:
            raise ValueError("'{0}' is not a valid coordinate reference frame "
                             "in astropy supported by DS9".format(self.coordsys))

    def _validate(self):
        """
        Checks whether all the attributes of this object is valid.
        """
        if self.region_type not in regions_attributes:
            raise ValueError("'{0}' is not a valid region type in this package"
                             .format(self.region_type))

        if self.coordsys not in valid_coordsys['DS9'] + valid_coordsys['CRTF']:
            raise ValueError("'{0}' is not a valid coordinate reference frame "
                             "in astropy".format(self.coordsys))


def to_shape_list(region_list, coordinate_system='fk5'):
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
        if isinstance(region, SkyRegion):
            reg_type = region.__class__.__name__[:-9].lower()
        else:
            reg_type = region.__class__.__name__[:-11].lower()

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
            if isinstance(val, Angle) or isinstance(val, u.Quantity) or isinstance(val, numbers.Number):
                new_coord.append(val)
            elif isinstance(val, PixCoord):
                new_coord.append(u.Quantity(val.x, u.dimensionless_unscaled))
                new_coord.append(u.Quantity(val.y, u.dimensionless_unscaled))
            else:
                new_coord.append(Angle(val.transform_to(frame).spherical.lon))
                new_coord.append(Angle(val.transform_to(frame).spherical.lat))

        meta = dict(region.meta)
        meta.update(region.visual)

        if reg_type == 'text':
            meta['text'] = meta.get('text', meta.pop('label', ''))

        include = region.meta.pop('include', True)

        shape_list.append(Shape(coordsys, reg_type, new_coord, meta, False,
                                include))

    return shape_list


def to_ds9_meta(shape_meta):
    """
    Makes the meta data DS9 compatible by filtering and mapping the valid keys

    Parameters
    ----------
    shape_meta: dict
        meta attribute of a `regions.Shape` object

    Returns
    -------
    meta : dict
        DS9 compatible meta dictionary
    """

    # meta keys allowed in DS9.
    valid_keys = ['symbol', 'include', 'tag', 'line', 'comment',
                  'name', 'select', 'highlite', 'fixed', 'label', 'text',
                  'edit', 'move', 'rotate', 'delete', 'source', 'background']

    # visual keys allowed in DS9
    valid_keys += ['color', 'dash', 'linewidth', 'font', 'dashlist',
                   'fill', 'textangle', 'symsize']

    # mapped to actual names in DS9
    key_mappings = {'symbol': 'point', 'linewidth': 'width', 'label': 'text'}

    meta = _to_io_meta(shape_meta, valid_keys, key_mappings)

    if 'font' in meta:
        meta['font'] += " {0} {1} {2}".format(shape_meta.get('fontsize', 12),
                                              shape_meta.get('fontstyle', 'normal'),
                                              shape_meta.get('fontweight', 'roman'))

    return meta


def to_crtf_meta(shape_meta):
    """
    Makes the meta data CRTF compatible by filtering and mapping the valid keys

    Parameters
    ----------
    shape_meta: dict
        meta attribute of a `regions.Region` object

    Returns
    -------
    meta : dict
        CRTF compatible meta dictionary

    """
    # please refer : https://casaguides.nrao.edu/index.php/CASA_Region_Format

    # meta keys allowed in CRTF
    valid_keys = ['label', 'include', 'frame', 'range', 'veltype',
                  'restfreq', 'coord', 'type', 'text', 'corr']

    # visual keys allowed in CRTF
    valid_keys += ['color', 'width', 'font', 'symthick', 'symsize', 'fontsize',
                   'fontstyle', 'usetex', 'labelpos', 'labeloff', 'linewidth',
                   'linestyle', 'symbol']

    key_mappings = {}

    return _to_io_meta(shape_meta, valid_keys, key_mappings)


def _to_io_meta(shape_meta, valid_keys, key_mappings):
    """
    This is used to make meta data compatible with a specific io
    by filtering and mapping to it's valid keys

    Parameters
    ----------
    shape_meta: dict
        meta attribute of a `regions.Region` object
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

    for key in shape_meta:
        if key in valid_keys:
            meta[key_mappings.get(key, key)] = shape_meta[key]

    return meta

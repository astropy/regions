# Licensed under a 3-clause BSD style license - see LICENSE.rst
from warnings import warn
import string
import numbers
import numpy as np

from astropy import units as u
from astropy import coordinates
from astropy.coordinates import Angle, SkyCoord
from astropy.io.registry import UnifiedReadWriteMethod
from astropy import log
from astropy.table import Table

from .. import shapes
from ..core import PixCoord, SkyRegion
from ..core.attributes import RegionMeta, RegionVisual
from .connect import ShapeListRead, ShapeListWrite
from .ds9.core import DS9RegionParserWarning, valid_symbols_ds9
from .crtf.core import CRTFRegionParserWarning

__all__ = ['ShapeList', 'Shape', 'to_shape_list', 'to_crtf_meta',
           'to_ds9_meta']


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
valid_coordsys['DS9'] += [f'wcs{x}' for x in string.ascii_lowercase]

# Maps astropy's coordinate frame names with their respective name in the file format.
coordsys_mapping = {'DS9': {x: x for x in valid_coordsys['DS9']},
                    'CRTF': {x: x.upper() for x in valid_coordsys['CRTF']}
                    }
# Coordinate frame names must be uppercase following the CASA CRTF syntax
coordsys_mapping['CRTF']['geocentrictrueecliptic'] = 'ECLIPTIC'
coordsys_mapping['CRTF']['fk5'] = 'J2000'
coordsys_mapping['CRTF']['fk4'] = 'B1950'
coordsys_mapping['CRTF']['supergalactic'] = 'SUPERGAL'

coordsys_mapping['DS9']['geocentrictrueecliptic'] = 'ECLIPTIC'


class ShapeList(list):
    """
    A class to hold a list of `~regions.Shape` objects.
    """

    # Unified I/O read and write methods from .connect
    read = UnifiedReadWriteMethod(ShapeListRead)
    write = UnifiedReadWriteMethod(ShapeListWrite)

    def to_regions(self):
        """
        Convert to a list of `~regions.Region` objects.
        """
        regions = []
        for shape in self:
            # Skip elliptical annulus for now
            if shape.region_type == 'ellipse' and len(shape.coord) > 5:
                msg = f'Skipping elliptical annulus {shape}'
                warn(msg, DS9RegionParserWarning)
                continue
            if shape.region_type in ['box'] and shape.format_type == 'CRTF':
                msg = f'Skipping {shape.region_type} {shape}'
                warn(msg, CRTFRegionParserWarning)
                continue
            log.debug(shape)
            region = shape.to_region()
            log.debug(region)
            regions.append(region)
        return regions

    def to_crtf(self, coordsys='fk5', fmt='.6f', radunit='deg'):
        """
        Convert to CRTF (CASA Region Text Format) region strings.

        Parameters
        ----------
        coordsys : str
            An Astropy coordinate system that overrides the coordinate
            system frame for all regions.

        fmt : str
            A python string format defining the output precision.
            Default is '.6f', which is accurate to 0.0036 arcseconds.

        radunit : str
            The unit of the radius.

        Returns
        -------
        region_string : str
            A CRFT region string.
        """
        crtf_strings = {
            'circle': '{0}circle[[{1:FMT}deg, {2:FMT}deg], {3:FMT}RAD]',
            'circleannulus': '{0}annulus[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}RAD, {4:FMT}RAD]]',
            # Make sure that width goes to minor axis and height to major axis
            'ellipse': '{0}ellipse[[{1:FMT}deg, {2:FMT}deg], [{4:FMT}RAD, {3:FMT}RAD], {5:FMT}deg]',
            'rectangle': '{0}rotbox[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}RAD, {4:FMT}RAD], {5:FMT}deg]',
            'polygon': '{0}poly[{1}]',
            'point': '{0}point[[{1:FMT}deg, {2:FMT}deg]]',
            'symbol': '{0}symbol[[{1:FMT}deg, {2:FMT}deg], {symbol}]',
            'text': '{0}text[[{1:FMT}deg, {2:FMT}deg], \'{text}\']',
            'line': '{0}line[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}deg, {4:FMT}deg]]'
                        }

        output = '#CRTFv0\n'

        if radunit == 'arcsec':
            # arcseconds are allowed for all but image coordinates
            if coordsys.lower() not in ('image',):
                radunitstr = '"'
            else:
                raise ValueError(
                    f'Radius unit arcsec not valid for coordsys {coordsys}')
        else:
            radunitstr = radunit

        for key, val in crtf_strings.items():
            crtf_strings[key] = val.replace("FMT", fmt).replace("RAD",
                                                                radunitstr)

        # CASA does not support global coordinate specification, even though the
        # documentation for the specification explicitly states that it does.
        output += f"global coord={coordsys_mapping['CRTF'][coordsys.lower()]}\n"

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
                shape.meta['label'] = f"'{shape.meta['label']}'"
            meta_str = ", ".join(f"{key}={val}" for key, val in
                                 shape.meta.items() if
                                 key not in ('include', 'comment', 'symbol',
                                             'coord', 'text', 'range', 'corr',
                                             'type'))

            # the first item should be the coordinates, since CASA cannot
            # recognize a region without an inline coordinate specification
            # It can be, but does not need to be, comma-separated at the start
            shape_coordsys = getattr(shape, 'coordsys')
            if shape_coordsys.lower() != coordsys.lower():
                if meta_str.strip():
                    meta_str = f"coord={coordsys_mapping['CRTF'][coordsys.lower()]}, " + meta_str
                else:
                    # if there is no metadata at all (above), the trailing comma is incorrect
                    meta_str = f"coord={coordsys_mapping['CRTF'][coordsys.lower()]}"

            if 'comment' in shape.meta:
                meta_str += ", " + shape.meta['comment']
            if 'range' in shape.meta:
                shape.meta['range'] = [str(str(x).replace(" ", "")) for x in
                                       shape.meta['range']]
                meta_str += f", range={shape.meta['range']}".replace("'", "")
            if 'corr' in shape.meta:
                meta_str += f", corr={shape.meta['corr']}".replace("'", "")

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
                vals = [f'[{x:{fmt}}deg, {y:{fmt}}deg]'
                        for x, y in zip(coord[::2], coord[1::2])]
                coord = ", ".join(vals)
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
                output += f"{line}, {meta_str}\n"
            else:
                output += f"{line}\n"

        return output

    def to_ds9(self, coordsys='fk5', fmt='.6f', radunit='deg'):
        """
        Convert to DS9 region strings.

        Parameters
        ----------
        coordsys : str
            An Astropy coordinate system that overrides the coordinate
            system frame for all regions.

        fmt : str
            A python string format defining the output precision.
            Default is '.6f', which is accurate to 0.0036 arcseconds.

        radunit : str
            The unit of the radius.

        Returns
        -------
        region_string : str
            A DS9 region string.
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
                raise ValueError(f'Radius unit arcsec not valid for coordsys {coordsys}')
        else:
            radunitstr = ''

        for key, val in ds9_strings.items():
            ds9_strings[key] = val.replace("FMT", fmt).replace("RAD", radunitstr)

        output += f'{coordsys}\n'

        for shape in self:

            shape.check_ds9()
            shape.meta = to_ds9_meta(shape.meta)

            # if unspecified, include is True.
            include = "-" if shape.include in (False, '-') else ""

            if 'point' in shape.meta:
                shape.meta['point'] = valid_symbols_reverse[shape.meta['point']]

            if 'symsize' in shape.meta:
                shape.meta['point'] += f" {shape.meta.pop('symsize')}"

            meta_str = " ".join(f"{key}={val}" for key, val in
                                shape.meta.items() if key not in ('include', 'tag', 'comment', 'font', 'text'))
            if 'tag' in shape.meta:
                meta_str += " " + " ".join([f"tag={tag}" for tag in shape.meta['tag']])
            if 'font' in shape.meta:
                meta_str += " " + f"font=\"{shape.meta['font']}\""
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
                coord = ",".join([f'{x:{fmt}}' for x in coord])
                line = ds9_strings['polygon'].format(include, coord)
            elif shape.region_type == 'ellipse':
                coord[2:] = [x / 2 for x in coord[2:]]
                if len(coord) % 2 == 1:
                    coord[-1] *= 2
                line = ds9_strings['ellipse'].format(include, *coord)
            else:
                line = ds9_strings[shape.region_type].format(include, *coord)

            if meta_str.strip():
                output += f"{line} # {meta_str}\n"
            else:
                output += f"{line}\n"

        return output

    def to_fits(self):
        """
        Convert to a `~astropy.table.Table` object.
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


class Shape:
    """
    Helper class to represent a DS9/CRTF Region.

    This serves as intermediate step in the parsing process.

    Parameters
    ----------
    coordsys : str
        An Astropy Coordinate system frame used in the region.

    region_type : str
        The type of the region (as defined in this package).

    coord : list of `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The region coordinates.

    meta : dict
        The meta attributes.

    composite : bool
        Flag indicting wheter the region is a Composite region.

    include : bool
        Flag indicating where to include or exclude the region.
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
        ss += f"\nType: {self.meta.get('type', 'reg')}"
        ss += f'\nCoord sys: {self.coordsys}'
        ss += f'\nRegion type: {self.region_type}'
        if self.region_type == 'symbol':
            ss += f"\nSymbol: {self.meta['symbol']}"
        if self.region_type == 'text':
            ss += f"\nText: {self.meta['text']}"
        ss += f'\nMeta: {self.meta}'
        ss += f'\nComposite: {self.composite}'
        ss += f'\nInclude: {self.include}'
        ss += '\n'
        return ss

    def convert_coords(self):
        """
        Process a list of coordinates.

        This mainly searches for a tuple of coordinates in the
        coordinate list and creates a SkyCoord or PixCoord object from
        them if appropriate for a given region type. This involves again
        some coordinate transformation, so this step could be moved to
        the parsing process.
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
        Convert to sky coordinates.
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
        Convert to pixel coordinates.
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
        Convert to a `~regions.Region` object.
        """
        coords = self.convert_coords()
        log.debug(coords)
        viz_keywords = ['color', 'dash', 'dashlist', 'width', 'font', 'symsize',
                        'symbol', 'symsize', 'fontsize', 'fontstyle', 'usetex',
                        'labelpos', 'labeloff', 'linewidth', 'linestyle',
                        'point', 'textangle', 'fontweight', 'symthick']

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
        Check for CRTF compatibility.
        """
        if self.region_type not in regions_attributes:
            raise ValueError(f"'{self.region_type}' is not a valid region "
                             "type")

        if self.coordsys not in valid_coordsys['CRTF']:
            raise ValueError(f"'{self.coordsys}' is not a valid coordinate "
                             "reference frame")

    def check_ds9(self):
        """
        Check for DS9 compatibility.
        """
        if self.region_type not in regions_attributes:
            raise ValueError(f"'{self.region_type}' is not a valid region "
                             "type")

        if self.coordsys not in valid_coordsys['DS9']:
            raise ValueError(f"'{self.coordsys}' is not a valid coordinate "
                             "reference frame")

    def _validate(self):
        """
        Check whether all the attributes of this object is valid.
        """
        if self.region_type not in regions_attributes:
            raise ValueError(f"'{self.region_type}' is not a valid region "
                             "type")

        if self.coordsys not in valid_coordsys['DS9'] + valid_coordsys['CRTF']:
            raise ValueError(f"'{self.coordsys}' is not a valid coordinate "
                             "reference frame")


def to_shape_list(region_list, coordinate_system='fk5'):
    """
    Convert a list of regions into a `~regions.ShapeList` object.

    Parameters
    ----------
    region_list : list
        A list of `regions.Region` objects.

    coordinate_system : str, optional
        The Astropy coordinate system frame in which all the coordinates
        present in the ``region_list`` will be converted. Default is
        'fk5'.

    Returns
    -------
    shape_list: `regions.ShapeList`
        A `~regions.ShapeList` object.
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
    Make the metadata DS9 compatible by filtering and mapping the valid
    keys.

    Parameters
    ----------
    shape_meta : dict
        The meta attribute of a `regions.Shape` object.

    Returns
    -------
    meta : dict
        A DS9 compatible meta dictionary.
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
        meta['font'] += (f" {shape_meta.get('fontsize', 12)} "
                         f"{shape_meta.get('fontstyle', 'normal')} "
                         f"{shape_meta.get('fontweight', 'roman')}")

    return meta


def to_crtf_meta(shape_meta):
    """
    Make the metadata CRTF compatible by filtering and mapping the valid
    keys.

    Parameters
    ----------
    shape_meta : dict
        A meta attribute of a `regions.Region` object.

    Returns
    -------
    meta : dict
        A CRTF compatible meta dictionary.
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
    Make metadata compatible with a specific IO by filtering and mapping
    to its valid keys.

    Parameters
    ----------
    shape_meta : dict
        A meta attribute of a `regions.Region` object.

    valid_keys : list
        The valid keys of a particular file format.

    key_mappings : dict
        A dictionary mapping of the actual name of the key in the
        format.

    Returns
    -------
    meta : dict
        An IO compatible meta dictionary according to ``valid_keys`` and
        ``key_mappings``.
    """
    meta = dict()

    for key in shape_meta:
        if key in valid_keys:
            meta[key_mappings.get(key, key)] = shape_meta[key]

    return meta

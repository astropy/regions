# Licensed under a 3-clause BSD style license - see LICENSE.rst

from warnings import warn

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.utils import deprecated
from astropy.wcs import WCS
import numpy as np

from ...core.registry import RegionsRegistry
from ...core import RegionList
from ..core import Shape, ShapeList, reg_mapping
from .core import (FITSRegionParserError, FITSRegionParserWarning,
                   language_spec)

__all__ = ['FITSRegionParser', 'read_fits']


class FITSRegionParser:
    """
    Parses a FITS Region table.

    Parameters
    ----------
    table : `~astropy.table.Table`
        A FITS region table.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.FITSRegionParserError`. 'warn' will raise a
        `~regions.FITSRegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).

    Examples
    --------
    >>> from regions import FITSRegionParser
    >>> from astropy.table import Table
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> filename = get_pkg_data_filename('data/fits_region.fits',
    ...                                  package='regions.io.fits.tests')
    >>> table = Table.read(filename)
    >>> parser = FITSRegionParser(table)
    >>> shapes = parser.shapes
    >>> regions = shapes.to_regions()
    >>> regions[5]
    <PointPixelRegion(center=PixCoord(x=341.0, y=345.0))>
    """

    valid_columns = ['X', 'Y', 'SHAPE', 'COMPONENT', 'R', 'ROTANG']

    def __init__(self, table, errors='strict'):
        if errors not in ('strict', 'ignore', 'warn'):
            raise ValueError('errors must be one of "strict", "ignore", or '
                             '"warn"')
        if not isinstance(table, Table):
            raise TypeError('The table should be an astropy table object')
        self.table = table
        self.errors = errors
        self.unit = {}
        self._shapes = {}

        self.parse_table()

    @property
    def shapes(self):
        shape_list = ShapeList()
        components = list(self._shapes.keys())
        components.sort()
        for component in components:
            shape_list += ShapeList(self._shapes[component])
        return shape_list

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, FITSRegionParserWarning)
        elif self.errors == 'strict':
            raise FITSRegionParserError(msg)

    def parse_table(self):
        for col in self.table.colnames:
            if col not in self.valid_columns:
                self._raise_error('This table has an invalid column name: '
                                  f'"{col}"')
            else:
                self.unit[col] = self.table[col].unit

        if 'COMPONENT' not in self.table.colnames:
            self.table['COMPONENT'] = np.ones([1, len(self.table)])

        for row in self.table:
            reg = _FITSRegionRowParser(row, self.unit, self.errors)
            component, shape = reg.parse()
            if component in self.shapes:
                self._shapes[component].append(shape)
            else:
                self._shapes[component] = [shape]


class _FITSRegionRowParser():
    """
    Parse a single row of a FITS region table.

    Parameters
    ----------
    row : `~astropy.table.Row`
        A row of the region table that is to be parsed.

    unit : dict
        The units of each column in the row.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.FITSRegionParserError`. 'warn' will raise a
        `~regions.FITSRegionParserWarning`, and 'ignore' will do
        nothing (i.e., be silent).
    """

    def __init__(self, row, unit, errors='strict'):
        self.errors = errors
        self.unit = unit
        self.row = row

        region_type = self._get_col_value('SHAPE0', 'POINT')[0]

        # Default region is POINT
        if region_type == '':
            region_type = 'POINT'

        if region_type[0] == '!':
            self.include = False
            region_type = region_type[1:]
        else:
            self.include = True

        region_type = region_type.strip().upper()

        if region_type in language_spec:
            self.region_type = region_type
        else:
            self._raise_error(f'"{region_type}" is not a valid FITS Region '
                              'type')

        self.component = str(self.row['COMPONENT'])

    def _get_col_value(self, colname, default=None):
        index = None

        if colname[-1].isdigit():
            index = int(colname[-1])
            colname = colname[: -1]

        try:
            val, unit = self.row[colname], self.unit[colname]

            if np.isscalar(val):
                val = np.array(val).reshape(1, )
            if index is not None:
                if index < len(val) and val[index] is not None:
                    return val[index], unit
                else:
                    raise ValueError(f'The column "{colname}" must have more '
                                     f'than {index} value for the region '
                                     f'{self.region_type}')
            else:
                return val, unit
        except KeyError:
            if default is None:
                self._raise_error(f'The column "{colname}" is missing in '
                                  'the table')
            else:
                return default

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, FITSRegionParserWarning)
        elif self.errors == 'strict':
            raise FITSRegionParserError(msg)

    def parse(self):
        coords = []

        for x in language_spec[self.region_type]:
            y, unit = self._get_col_value(x)
            coords.append(self._parse_value(y, unit))

        meta = {'tag': self.component}

        if self.region_type == 'POLYGON':
            coords_new = []
            for x, y in zip(coords[0], coords[1]):
                coords_new += [x, y]
            coords = coords_new
        elif self.region_type == 'BOX':
            # Add a 0-degree rotation to turn it into ROTBOX
            coords.append(0.0 * u.deg)

        region_type = self.region_type.lower()
        if region_type in reg_mapping['FITS_REGION']:
            region_type = reg_mapping['FITS_REGION'][region_type]
        else:
            self._raise_error(f'"{self.region_type}" is currently not '
                              'supported')

        return self.component, Shape('physical', region_type,
                                     coords, meta, False, False)

    def _parse_value(self, val, unit):
        units = dict(pix=u.dimensionless_unscaled,
                     deg=u.deg,
                     rad=u.rad,
                     )

        if unit is not None:
            return val * units.get(str(unit), unit)
        else:
            self._raise_error(f'The unit: {unit} is invalid')


@RegionsRegistry.register('RegionList', 'read', 'fits')
def read_fits(filename, errors='strict', cache=False):
    """
    Read a FITS region file, converting a FITS regions table to a list
    of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The file path.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.FITSRegionParserError`. 'warn' will raise a
        `~regions.FITSRegionParserWarning`, and 'ignore' will do
        nothing (i.e., be silent).

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.

    Examples
    --------
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from regions import read_fits
    >>> file_read = get_pkg_data_filename('data/fits_region.fits',
    ...                                   package='regions.io.fits.tests')
    >>> regions = read_fits(file_read)
    """
    regions = []

    with fits.open(filename) as hdul:
        for hdu in hdul:
            if hdu.name == 'REGION':
                table = Table.read(hdu)
                wcs = WCS(hdu.header, keysel=['image', 'binary', 'pixel'])
                parser = FITSRegionParser(table, errors)
                regions_list = parser.shapes.to_regions()
                for reg in regions_list:
                    regions.append(reg.to_sky(wcs))

    return RegionList(regions)


@deprecated('0.5', alternative='read_fits')
def read_fits_region(filename, errors='strict'):
    """
    Read a FITS region file, converting a FITS regions table to a list
    of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The file path.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.FITSRegionParserError`. 'warn' will raise a
        `~regions.FITSRegionParserWarning`, and 'ignore' will do
        nothing (i.e., be silent).

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.
    """
    return read_fits(filename, errors=errors)

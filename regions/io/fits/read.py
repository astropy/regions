# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

from warnings import warn

from astropy import units as u
import numpy as np

from .core import FITSRegionParserError, FITSRegionParserWarning, language_spec
from ..core import Shape, ShapeList, reg_mapping

__all__ = ['FITSRegionParser']


class FITSRegionParser(object):
    """
    table: `~astropy.table.Table` object
        A fits region table
    errors : ``warn``, ``ignore``, ``strict``
        The error handling scheme to use for handling parsing errors.
        The default is 'strict', which will raise a ``CRTFRegionParserError``.
        `warn`` will raise a ``CRTFRegionParserWarning``, and ``ignore`` will do nothing
        (i.e., be silent).
    """

    valid_columns = ['X', 'Y', 'SHAPE', 'COMPONENT', 'R', 'ROTANG']

    def __init__(self, table, errors='strict'):
        self.errors = errors
        self.table = table
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
                self._raise_error('This table has an invalid col name: {}'
                                  .format(col))
            else:
                self.unit[col] = self.table[col].unit

        if 'COMPONENT' not in self.table.colnames:
            self.table['COMPONENT'] = np.ones([1, len(self.table)])

        for row in self.table:
            reg = FITSRegionRowParser(row, self.unit, self.errors)
            component, shape = reg.parse()
            if component in self.shapes:
                self._shapes[component].append(shape)
            else:
                self._shapes[component] = [shape]


class FITSRegionRowParser():

    def __init__(self, row, unit, errors='strict'):
        self.errors = errors
        self.unit = unit
        self.row = row

        region_type = self._get_col_value('SHAPE0', 'POINT')[0]
        if region_type[0] == '!':
            self.include = False
            region_type = region_type[1:]
        else:
            self.include = True

        region_type = region_type.strip()

        if region_type in language_spec:
            self.region_type = region_type
        else:
            self._raise_error("{} is not a valid FITS Region type"
                              .format(region_type))

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
                if index < len(val) and val[index] != 0:
                    print("here")
                    return val[index], unit
                else:
                    raise ValueError("The {0} must have more than {1} value for the "
                                     "region {2}".format(colname, index,
                                                        self.region_type))
            else:
                return val, unit
        except KeyError:
            if default is None:
                self._raise_error("The {0} is missing in the table"
                                  .format(colname))
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
            print(x, y)
            coords.append(self._parse_value(y, unit))

        meta = {'tag': self.component}

        return self.component, Shape('physical',
                                     reg_mapping['FITS_REGION'][self.region_type.lower()],
                                     coords, meta, False, False)

    def _parse_value(self, val, unit):

        units = dict(pix=u.dimensionless_unscaled,
                     deg=u.deg,
                     rad=u.rad,
                     )

        if unit is not None:
            return val * units.get(str(unit), unit)
        else:
            self._raise_error("The unit: {} is invalid".format(unit))

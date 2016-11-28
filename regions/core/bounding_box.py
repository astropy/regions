# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = ['BoundingBox']


class BoundingBox(object):
    """
    Rectangular bounding box in integer pixel indices (not floats).

    Parameters
    ----------
    ixmin, ixmax, iymin, iymax : int
        The bounding indices.
        Note that the upper values (iymax and ixmax) are
        exclusive as for normal slices in Python.

    Examples
    --------

    >>> from regions import BoundingBox
    # Constructing a BoundingBox like this is cryptic:
    >>> bbox = BoundingBox(1, 10, 2, 20)
    # Better use keyword arguments for readability:
    >>> bbox = BoundingBox(ixmin=1, ixmax=10, iymin=2, iymax=20)
    >>> bbox  # nice repr, useful for interactive work
    BoundingBox(ixmin=1, ixmax=10, iymin=2, iymax=20)
    # sometimes it's useful to check if two bboxes are the same
    >>> bbox == BoundingBox(ixmin=1, ixmax=10, iymin=2, iymax=20)
    True
    >>> bbox == BoundingBox(ixmin=99, ixmax=10, iymin=2, iymax=20)
    False
    # Shape and slices can be useful when working with cutout numpy arrays
    >>> bbox.shape  # numpy order: (y, x)
    (18, 9)
    >>> bbox.slices  # numpy order: (y, x)
    (slice(2, 20, None), slice(1, 10, None))
    # Extent is useful when plotting bbox with matplotlib
    >>> bbox.extent  # matplotlib order: (x, y)
    (0.5, 9.5, 1.5, 19.5)
    >>> print(bbox.as_patch())
    Rectangle(0.5,1.5;9x18)
    """

    def __init__(self, ixmin, ixmax, iymin, iymax):
        self.ixmin = ixmin
        self.ixmax = ixmax
        self.iymin = iymin
        self.iymax = iymax

    def __eq__(self, other):
        if not isinstance(other, BoundingBox):
            raise TypeError('Can only compare BoundingBox to other BoundingBox')

        return (
            (self.ixmin == other.ixmin) and
            (self.ixmax == other.ixmax) and
            (self.iymin == other.iymin) and
            (self.iymax == other.iymax)
        )

    def __repr__(self):
        data = self.__dict__
        data['name'] = self.__class__.__name__
        fmt = '{name}(ixmin={ixmin}, ixmax={ixmax}, iymin={iymin}, iymax={iymax})'
        return fmt.format(**data)

    @property
    def shape(self):
        """
        The shape of the bounding box.

        Numpy axis order `(y, x)`.
        """
        return self.iymax - self.iymin, self.ixmax - self.ixmin

    @property
    def slices(self):
        """
        The bounding box as a pair of `slice` objects.

        Numpy axis order `(y, x)`. Can be used to index into Numpy arrays.
        """
        return (
            slice(self.iymin, self.iymax),
            slice(self.ixmin, self.ixmax),
        )

    @property
    def extent(self):
        """
        The 'extent' of the mask, i.e the bounding box from the bottom left
        corner of the lower left pixel to the upper right corner of the upper
        right pixel.

        This can be used for example when plotting using Matplotlib.
        """
        return (
            self.ixmin - 0.5,
            self.ixmax - 0.5,
            self.iymin - 0.5,
            self.iymax - 0.5,
        )

    def as_patch(self, **kwargs):
        """
        Return a Matplotlib patch that represents the bounding box.

        TODO: show full code example how to add it to a plot
        """
        from matplotlib.patches import Rectangle
        return Rectangle(
            xy=(self.ixmin - 0.5, self.iymin - 0.5),
            width=self.ixmax - self.ixmin,
            height=self.iymax - self.iymin,
            **kwargs
        )

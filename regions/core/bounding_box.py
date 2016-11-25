class BoundingBox(object):
    """
    A simple class used to represent a bounding box in integer pixel
    coordinates.

    Parameters
    ----------
    ixmin, ixmax, iymin, iymax : int
        The bounding indices. Note that the upper values (iymax and ixmax) are
        exlusive as for normal slices in Python.
    """

    def __init__(self, ixmin, ixmax, iymin, iymax):
        self.ixmin = ixmin
        self.ixmax = ixmax
        self.iymin = iymin
        self.iymax = iymax

    @property
    def shape(self):
        """
        The shape of the bounding box
        """
        return self.iymax - self.iymin, self.ixmax - self.ixmin

    @property
    def slices(self):
        """
        The bounding box as a pair of `slice` objects.
        """
        return (slice(self.iymin, self.iymax),
                slice(self.ixmin, self.ixmax))

    @property
    def extent(self):
        """
        The 'extent' of the mask, i.e the bounding box from the bottom left
        corner of the lower left pixel to the upper right corner of the upper
        right pixel. This can be used for example when plotting using Matplotlib.
        """
        return self.ixmin - 0.5, self.ixmax - 0.5, self.iymin - 0.5, self.iymax - 0.5

    def as_patch(self, **kwargs):
        """
        Return a Matplotlib patch that represents the bounding box.
        """
        from matplotlib.patches import Rectangle
        return Rectangle((self.ixmin - 0.5, self.iymin - 0.5),
                         self.ixmax - self.ixmin,
                         self.iymax - self.iymin, **kwargs)

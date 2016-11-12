class BoundingBox(object):
    """
    A simple class used to represent a bounding box in integer pixel
    coordinates.

    Parameters
    ----------
    jmin, jmax, imin, imax : int
        The bounding indices. Note that the upper values (jmax and imax) are
        exlusive as for normal slices in Python.
    """

    def __init__(self, jmin, jmax, imin, imax):
        self.jmin = jmin
        self.jmax = jmax
        self.imin = imin
        self.imax = imax

    @property
    def shape(self):
        """
        The shape of the bounding box
        """
        return self.jmax - self.jmin, self.imax - self.imin

    @property
    def slices(self):
        """
        The bounding box as a pair of `slice` objects.
        """
        return (slice(self.jmin, self.jmax),
                slice(self.imin, self.imax))

__all__ = ['Mask']

class Mask(object):
    """
    A class that can be used to represent a mask that can then be applied to
    an array.

    Rather than require the full mask for the array to be specified, this class
    allows you to represent the mask using a cutout and an origin, which will be
    far more efficient when small regions are masked.

    Parameters
    ----------
    mask : `numpy.ndarray`
        The mask values. This should be a floating-point array where a value
        of 0 indicates no overlap and a value of 1 means full overlap.
    origin : tuple of int, optional
        The position of the lower left pixel. If not specified, the mask needs
        to be the same shape as the array it is then applied to.
    """

    def __init__(self, mask, origin=None):
        self.mask = mask
        self.origin = origin

    def apply(self, array):
        if self.origin is None:
            if array.shape != self.mask.shape:
                raise ValueError("No origin was specified for mask and mask shape does not match array shape")
            else:
                return array * self.mask
        else:
            pass

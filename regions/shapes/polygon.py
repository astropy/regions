# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord

from ..core import PixelRegion, SkyRegion, RegionMask, BoundingBox, PixCoord
from .._geometry import polygonal_overlap_grid
from .._geometry.pnpoly import points_in_polygon
from ..core.attributes import OneDPix, OneDSky, RegionMeta, RegionVisual

__all__ = ['PolygonPixelRegion', 'PolygonSkyRegion']


class PolygonPixelRegion(PixelRegion):
    """
    A polygon in pixel coordinates.

    Parameters
    ----------
    vertices : `~regions.PixCoord`
        The vertices of the polygon
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.

    Examples
    --------

    .. plot::
        :include-source:

        import numpy as np
        from astropy.coordinates import Angle
        from regions import PixCoord, PolygonPixelRegion
        import matplotlib.pyplot as plt

        x, y = [45, 45, 55, 60], [75, 70, 65, 75]
        fig, ax = plt.subplots(1, 1)

        vertices = PixCoord(x=x, y=y)
        reg = PolygonPixelRegion(vertices=vertices)
        patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.xlim(30, 80)
        plt.ylim(50, 80)
        ax.set_aspect('equal')
    """
    _params = ('vertices',)
    vertices = OneDPix('vertices')

    def __init__(self, vertices, meta=None, visual=None):
        self.vertices = vertices
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @staticmethod
    def _find_slope(point_1, point_2):
        """
        Finds the slope of a line between point_1 and point_2

        Parameters
        ----------
        point_1: PixCoord(x, y) -> contains pixel coordinates of point_1
        point_2: PixCoord(x, y) -> contains pixel coordinates of point_2

        Returns
        -------
        float: slope value of the line with

        * Returns "None" if the line is perpendicular (Zero division)
        """

        del_y = point_2.y - point_1.y
        del_x = point_2.x - point_1.x

        if float(del_x) == 0.0:
            slope = None
        else:
            slope = del_y / del_x

        return slope

    @staticmethod
    def _area_under_line(point_1, point_2):
        """
        Finds the area under a line. (The area between the line and the x-axis.)

        Parameters
        ----------
        point_1: PixCoord(x, y) -> contains pixel coordinates of point_1
        point_2: PixCoord(x, y) -> contains pixel coordinates of point_2

        Returns
        -------
        float: Area value
        """

        area = abs(point_2.x - point_1.x) * (point_2.y + point_1.y) * 0.5

        return area

    @staticmethod
    def _find_intersection(line_1, line_2):
        """
        Finds the intersection coordinates of two given lines.
        If the intersection point goes beyond the limits of the lines, returns (None, None)

        line equation:
        y = m*x + n

        (m -> slope)
        (n -> constant)

        m and n values will be None in case of zero division.
        (If the line is perpendicular to the x-axis there will be zero division.)

        ----------
        line_1
        line_2

        Returns
        -------
        float, float: x and y coordinate of intersection point.
        """

        # a "near zero" value for comparing float values
        sensitivity = 1E-8

        # line: 1, point: 1
        l1x1, l1y1 = line_1[0].x, line_1[0].y
        # line: 1, point: 2
        l1x2, l1y2 = line_1[1].x, line_1[1].y

        # slope (m) and constant (n) values for line 1;
        if abs(float(l1x2 - l1x1)) > sensitivity:
            m1 = (l1y2 - l1y1) / (l1x2 - l1x1)
            n1 = -((l1y2 - l1y1) / (l1x2 - l1x1)) * l1x2 + l1y2
        else:
            m1 = None
            n1 = None

        # line: 2, point: 1
        l2x1, l2y1 = line_2[0].x, line_2[0].y
        # line: 2, point: 2
        l2x2, l2y2 = line_2[1].x, line_2[1].y

        # slope (m) and constant (n) values for line 2;
        if abs(float(l2x2 - l2x1)) > sensitivity:
            m2 = (l2y2 - l2y1) / (l2x2 - l2x1)
            n2 = -((l2y2 - l2y1) / (l2x2 - l2x1)) * l2x2 + l2y2
        else:
            m2 = None
            n2 = None

        # initial values for x and y coordinates of intersection point
        x, y = None, None

        # if m1 and m2 are not perpendicular;
        if (m1 is not None) and (m2 is not None):

            # if they are not parallel;
            if abs(float(m1 - m2)) > sensitivity:
                # intersection coordinates would be:
                x = (n2 - n1) / (m1 - m2)
                y = m1 * x + n1

        # if m1 is perpendicular and m2 is not;
        if m1 is None and m2 is not None:
            # In this case, the equation of the line 1 will be like this:
            # x = constant
            # So the constant x value will also be the x value of intersection point.
            x = l1x1
            y = m2 * x + n2

        # if m2 is perpendicular and m1 is not;
        if m2 is None and m1 is not None:
            # In this case, the equation of the line 2 will be like this:
            # x = constant
            # So the constant x value will also be the x value of intersection point.
            x = l2x1
            y = m1 * x + n1

        if x is not None and y is not None:
            # These conditions ensures the (x,y) intersection point doesn't go beyond the line limits.
            condition_1 = min(l1x1, l1x2) < x < max(l1x1, l1x2)
            condition_2 = min(l1y1, l1y2) < y < max(l1y1, l1y2)
            condition_3 = min(l2x1, l2x2) < x < max(l2x1, l2x2)
            condition_4 = min(l2y1, l2y2) < y < max(l2y1, l2y2)

            if (condition_1 and condition_2) and (condition_3 and condition_4):
                return x, y
            else:
                return None, None
        else:
            return None, None

    def _split_intersection(self, intersection_point, line_indexes):
        """
        This method splits the intersecting lines into two parts from the intersection point
        and creates a new polygon with these additional points.
        Returns the new x and y values as separate lists.

        Parameters
        ----------
        intersection_point: (x, y) values of the intersection point
        line_indexes: index values of intersecting lines. (Example: [1,3] -> 1. and 3. lines are intersecting.)

        Returns:
        -------
        list(float), list(float): list of x (float) values and list of y (float) values

        Example:
        -------

        points_of_my_polygon: [p1, p2, p3, p4, p5, p6, p7, p8]

        If the [p2,p3] line and  the [p5,p6] line are intersecting at point p';

        -Step 1
            Add p' between p2-p3 and p5-p6
            [p1, p2, p', p3, p4, p5, p', p6, p7, p8]

        -Step 2
            Reverse the order of points between two p'
            [p1, p2, p', p5, p4, p3, p', p6, p7, p8]


        In case of multiple intersections, the method will be called separately for each.

        """
        x_list = self.vertices.x
        y_list = self.vertices.y
        p_count = len(x_list)

        x_i, y_i = intersection_point
        index_1, index_2 = line_indexes

        x_left, x_mid, x_right = x_list[:index_1 + 1], x_list[index_1 + 1:int(index_2 + 1 % p_count)], x_list[int(
            index_2 + 1 % p_count):]
        y_left, y_mid, y_right = y_list[:index_1 + 1], y_list[index_1 + 1:int(index_2 + 1 % p_count)], y_list[int(
            index_2 + 1 % p_count):]

        new_x_list = np.array([*x_left, x_i, *x_mid[::-1], x_i, *x_right])
        new_y_list = np.array([*y_left, y_i, *y_mid[::-1], y_i, *y_right])

        return new_x_list, new_y_list

    @property
    def area(self):
        """
        Calculates the area of a polygon.
        The polygon can be concave or convex with multiple self intersecting lines.

        Returns:
        -------
        float: Area of the polygon
        """

        # The vertices will be changed if there is an intersection.
        # The original values should be stored before the process.
        initial_x_vertices = self.vertices.x.copy()
        initial_y_vertices = self.vertices.y.copy()

        # Searching for intersections
        search_done = False
        while not search_done:
            self.index_counter = 0

            for i in range(self.index_counter, len(self.vertices)):
                p1 = self.vertices[i]
                p2 = self.vertices[int((i + 1) % len(self.vertices))]

                for j in range(i + 1, len(self.vertices)):
                    p3 = self.vertices[j]
                    p4 = self.vertices[int((j + 1) % len(self.vertices))]

                    x, y = self._find_intersection([p1, p2], [p3, p4])

                    if x is not None and y is not None:
                        self.vertices.x, self.vertices.y = self._split_intersection((x, y), (i, j))
                        self.index_counter = i
                        break

                if (i + 1) == len(self.vertices):
                    search_done = True
                    del self.index_counter
                    break

        # start calculating area
        total_area = 0.0
        p_count = len(self.vertices)
        for i in range(p_count):
            p1 = self.vertices[i]
            p2 = self.vertices[int((i + 1) % p_count)]

            slope = self._find_slope(p1, p2)

            if slope is not None:
                if p2.x > p1.x:
                    total_area += self._area_under_line(p1, p2)
                else:
                    total_area -= self._area_under_line(p1, p2)

        total_area = abs(total_area)

        # When the calculation finishes, original vertices are put back.
        self.vertices.x = initial_x_vertices
        self.vertices.y = initial_y_vertices

        return total_area

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, 'pixcoord')
        x = np.atleast_1d(np.asarray(pixcoord.x, dtype=float))
        y = np.atleast_1d(np.asarray(pixcoord.y, dtype=float))
        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        shape = x.shape
        mask = points_in_polygon(x.flatten(), y.flatten(), vx, vy).astype(bool)
        in_poly = mask.reshape(shape)
        if self.meta.get('include', True):
            return in_poly
        else:
            return np.logical_not(in_poly)

    def to_sky(self, wcs):
        vertices_sky = pixel_to_skycoord(self.vertices.x, self.vertices.y, wcs)
        return PolygonSkyRegion(vertices=vertices_sky)

    @property
    def bounding_box(self):
        xmin = self.vertices.x.min()
        xmax = self.vertices.x.max()
        ymin = self.vertices.y.min()
        ymax = self.vertices.y.max()
        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):

        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5
        xmax = float(bbox.ixmax) - 0.5
        ymin = float(bbox.iymin) - 0.5
        ymax = float(bbox.iymax) - 0.5

        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        fraction = polygonal_overlap_grid(
            xmin, xmax, ymin, ymax,
            nx, ny, vx, vy, use_exact, subpixels,
        )

        return RegionMask(fraction, bbox=bbox)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Matplotlib patch object for this region (`matplotlib.patches.Polygon`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs : `dict`
            All keywords that a `~matplotlib.patches.Polygon` object accepts

        Returns
        -------
        patch : `~matplotlib.patches.Polygon`
            Matplotlib polygon patch
        """
        from matplotlib.patches import Polygon
        xy = np.vstack([self.vertices.x - origin[0],
                        self.vertices.y - origin[1]]).transpose()

        mpl_params = self.mpl_properties_default('patch')
        mpl_params.update(kwargs)

        return Polygon(xy=xy, **mpl_params)

    def rotate(self, center, angle):
        """Make a rotated region.

        Rotates counter-clockwise for positive ``angle``.

        Parameters
        ----------
        center : `PixCoord`
            Rotation center point
        angle : `~astropy.coordinates.Angle`
            Rotation angle

        Returns
        -------
        region : `PolygonPixelRegion`
            Rotated region (an independent copy)
        """
        vertices = self.vertices.rotate(center, angle)
        return self.copy(vertices=vertices)


class PolygonSkyRegion(SkyRegion):
    """
    A polygon defined using vertices in sky coordinates.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the polygon
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """
    _params = ('vertices',)
    vertices = OneDSky('vertices')

    def __init__(self, vertices, meta=None, visual=None):
        self.vertices = vertices
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        x, y = skycoord_to_pixel(self.vertices, wcs)
        vertices_pix = PixCoord(x, y)
        return PolygonPixelRegion(vertices_pix)

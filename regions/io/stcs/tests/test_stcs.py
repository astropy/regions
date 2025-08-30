# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord

from regions.core import PixCoord, Regions
from regions.io.stcs.connect import is_stcs
from regions.io.stcs.core import (STCSParserError, parse_coordinate_frame,
                                  parse_numbers, parse_unit, validate_stcs_string)
from regions.io.stcs.read import _parse_stcs
from regions.io.stcs.write import _serialize_stcs
from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
                            EllipsePixelRegion, EllipseSkyRegion,
                            PolygonPixelRegion, PolygonSkyRegion,
                            PointPixelRegion, PointSkyRegion,
                            RectanglePixelRegion, RectangleSkyRegion)


class TestSTCSCore:
    """Test core STC-S functionality."""

    def test_validate_stcs_string(self):
        """Test STC-S string validation."""
        # Valid strings
        assert validate_stcs_string("Circle ICRS BARYCENTER 180.0 10.0 0.5")
        assert validate_stcs_string("Position FK5 GEOCENTER 45.0 -30.0")
        assert validate_stcs_string("polygon icrs barycenter 0 0 1 1 2 2")

        # Invalid strings
        assert not validate_stcs_string("not an stcs string")
        assert not validate_stcs_string("")
        assert not validate_stcs_string(None)

    def test_parse_coordinate_frame(self):
        """Test coordinate frame parsing."""
        frame, refpos = parse_coordinate_frame("Circle ICRS BARYCENTER 180.0 10.0 0.5")
        assert frame == 'ICRS'
        assert refpos == 'BARYCENTER'

        frame, refpos = parse_coordinate_frame("Position FK5 GEOCENTER 45.0 -30.0")
        assert frame == 'FK5'
        assert refpos == 'GEOCENTER'

        # Test defaults
        frame, refpos = parse_coordinate_frame("Circle 180.0 10.0 0.5")
        assert frame == 'ICRS'
        assert refpos == 'UNKNOWN'

    def test_parse_unit(self):
        """Test unit parsing."""
        unit = parse_unit("Circle ICRS unit deg 180.0 10.0 0.5")
        assert unit == 'degree'

        unit = parse_unit("Circle ICRS unit arcsec 180.0 10.0 0.5")
        assert unit == 'arcsec'

        # Test default
        unit = parse_unit("Circle ICRS 180.0 10.0 0.5", frame='ICRS')
        assert unit == 'degree'

        unit = parse_unit("Circle IMAGE 180.0 10.0 0.5", frame='IMAGE')
        assert unit == 'pixel'

    def test_parse_numbers(self):
        """Test number parsing."""
        numbers = parse_numbers("180.0 10.0 0.5")
        assert numbers == [180.0, 10.0, 0.5]

        numbers = parse_numbers("1.5e-3 -45.0 90")
        assert numbers == [1.5e-3, -45.0, 90.0]

        numbers = parse_numbers("")
        assert numbers == []


class TestSTCSConnect:
    """Test STC-S file identification."""

    def test_is_stcs_by_extension(self):
        """Test file identification by extension."""
        assert is_stcs('write', 'test.stcs')
        assert is_stcs('write', 'test.stc')
        assert not is_stcs('write', 'test.ds9')
        assert not is_stcs('write', 'test.txt')

        assert is_stcs('read', 'test.stcs')
        assert is_stcs('read', 'test.stc')
        assert is_stcs('read', 'test.stcs.txt')


class TestSTCSParsing:
    """Test STC-S parsing functionality."""

    def test_parse_circle_sky(self):
        """Test parsing sky circle regions."""
        stcs_str = "Circle ICRS BARYCENTER 180.0 10.0 0.5"
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 1
        region = regions[0]
        assert isinstance(region, CircleSkyRegion)
        assert region.center.ra.degree == 180.0
        assert region.center.dec.degree == 10.0
        assert region.radius.degree == 0.5

    def test_parse_circle_pixel(self):
        """Test parsing pixel circle regions."""
        stcs_str = "Circle IMAGE UNKNOWN 100.5 200.3 15.0"
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 1
        region = regions[0]
        assert isinstance(region, CirclePixelRegion)
        assert region.center.x == 100.5
        assert region.center.y == 200.3
        assert region.radius == 15.0

    def test_parse_ellipse_sky(self):
        """Test parsing sky ellipse regions."""
        stcs_str = "Ellipse ICRS BARYCENTER 180.0 10.0 0.5 0.3 45.0"
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 1
        region = regions[0]
        assert isinstance(region, EllipseSkyRegion)
        assert region.center.ra.degree == 180.0
        assert region.center.dec.degree == 10.0
        assert region.width.degree == 1.0  # 2 * semi_major
        assert region.height.degree == 0.6  # 2 * semi_minor
        assert region.angle.degree == 45.0

    def test_parse_box_sky(self):
        """Test parsing sky box regions."""
        stcs_str = "Box ICRS BARYCENTER 180.0 10.0 1.0 0.5 30.0"
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 1
        region = regions[0]
        assert isinstance(region, RectangleSkyRegion)
        assert region.center.ra.degree == 180.0
        assert region.center.dec.degree == 10.0
        assert region.width.degree == 1.0
        assert region.height.degree == 0.5
        assert region.angle.degree == 30.0

    def test_parse_polygon_sky(self):
        """Test parsing sky polygon regions."""
        stcs_str = "Polygon ICRS BARYCENTER 179.0 9.0 181.0 9.0 181.0 11.0 179.0 11.0"
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 1
        region = regions[0]
        assert isinstance(region, PolygonSkyRegion)
        assert len(region.vertices) == 4
        assert region.vertices[0].ra.degree == 179.0
        assert region.vertices[0].dec.degree == 9.0

    def test_parse_position_sky(self):
        """Test parsing sky position regions."""
        stcs_str = "Position ICRS BARYCENTER 180.0 10.0"
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 1
        region = regions[0]
        assert isinstance(region, PointSkyRegion)
        assert region.center.ra.degree == 180.0
        assert region.center.dec.degree == 10.0

    def test_parse_multiple_regions(self):
        """Test parsing multiple regions."""
        stcs_str = """
        Circle ICRS BARYCENTER 180.0 10.0 0.5
        Position FK5 GEOCENTER 45.0 -30.0
        """
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 2
        assert isinstance(regions[0], CircleSkyRegion)
        assert isinstance(regions[1], PointSkyRegion)

    def test_parse_with_comments(self):
        """Test parsing with comment lines."""
        stcs_str = """
        # This is a comment
        Circle ICRS BARYCENTER 180.0 10.0 0.5
        # Another comment
        Position FK5 GEOCENTER 45.0 -30.0
        """
        regions = _parse_stcs(stcs_str)

        assert len(regions) == 2

    def test_parse_invalid_stcs(self):
        """Test parsing invalid STC-S strings."""
        with pytest.raises(STCSParserError):
            _parse_stcs("not a valid stcs string")

        with pytest.raises(STCSParserError):
            _parse_stcs("")

    def test_parse_insufficient_parameters(self):
        """Test parsing with insufficient parameters."""
        # Circle needs at least 3 parameters
        with pytest.raises(STCSParserError):
            _parse_stcs("Circle ICRS BARYCENTER 180.0")

        # Ellipse needs at least 5 parameters
        with pytest.raises(STCSParserError):
            _parse_stcs("Ellipse ICRS BARYCENTER 180.0 10.0")


class TestSTCSSerialization:
    """Test STC-S serialization functionality."""

    def test_serialize_circle_sky(self):
        """Test serializing sky circle regions."""
        center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
        region = CircleSkyRegion(center=center, radius=0.5 * u.degree)

        stcs_str = _serialize_stcs(region)
        expected = "Circle ICRS BARYCENTER 180 10 0.5"
        assert stcs_str == expected

    def test_serialize_circle_pixel(self):
        """Test serializing pixel circle regions."""
        center = PixCoord(100.5, 200.3)
        region = CirclePixelRegion(center=center, radius=15.0)

        stcs_str = _serialize_stcs(region)
        expected = "Circle IMAGE UNKNOWN 100.5 200.3 15"
        assert stcs_str == expected

    def test_serialize_ellipse_sky(self):
        """Test serializing sky ellipse regions."""
        center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
        region = EllipseSkyRegion(center=center, width=1.0 * u.degree,
                                 height=0.6 * u.degree, angle=45.0 * u.degree)

        stcs_str = _serialize_stcs(region)
        expected = "Ellipse ICRS BARYCENTER 180 10 0.5 0.3 45"
        assert stcs_str == expected

    def test_serialize_box_sky(self):
        """Test serializing sky box regions."""
        center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
        region = RectangleSkyRegion(center=center, width=1.0 * u.degree,
                                   height=0.5 * u.degree, angle=30.0 * u.degree)

        stcs_str = _serialize_stcs(region)
        expected = "Box ICRS BARYCENTER 180 10 1 0.5 30"
        assert stcs_str == expected

    def test_serialize_polygon_sky(self):
        """Test serializing sky polygon regions."""
        vertices = SkyCoord([179.0, 181.0, 181.0, 179.0],
                           [9.0, 9.0, 11.0, 11.0],
                           unit='degree', frame='icrs')
        region = PolygonSkyRegion(vertices=vertices)

        stcs_str = _serialize_stcs(region)
        expected = "Polygon ICRS BARYCENTER 179 9 181 9 181 11 179 11"
        assert stcs_str == expected

    def test_serialize_position_sky(self):
        """Test serializing sky position regions."""
        center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
        region = PointSkyRegion(center=center)

        stcs_str = _serialize_stcs(region)
        expected = "Position ICRS BARYCENTER 180 10"
        assert stcs_str == expected

    def test_serialize_multiple_regions(self):
        """Test serializing multiple regions."""
        center1 = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
        region1 = CircleSkyRegion(center=center1, radius=0.5 * u.degree)

        center2 = SkyCoord(45.0, -30.0, unit='degree', frame='fk5')
        region2 = PointSkyRegion(center=center2)

        regions = Regions([region1, region2])
        stcs_str = _serialize_stcs(regions)

        lines = stcs_str.strip().split('\n')
        assert len(lines) == 2
        assert "Circle ICRS BARYCENTER 180 10 0.5" in lines[0]
        assert "Position FK5 BARYCENTER 45 -30" in lines[1]


class TestSTCSRoundTrip:
    """Test round-trip conversion between regions and STC-S."""

    def test_circle_roundtrip(self):
        """Test circle region round-trip conversion."""
        original_stcs = "Circle ICRS BARYCENTER 180.0 10.0 0.5"
        regions = _parse_stcs(original_stcs)
        serialized = _serialize_stcs(regions)

        # Parse again to verify consistency
        regions2 = _parse_stcs(serialized)

        assert len(regions) == len(regions2) == 1
        r1, r2 = regions[0], regions2[0]

        assert isinstance(r1, CircleSkyRegion)
        assert isinstance(r2, CircleSkyRegion)
        assert abs(r1.center.ra.degree - r2.center.ra.degree) < 1e-10
        assert abs(r1.center.dec.degree - r2.center.dec.degree) < 1e-10
        assert abs(r1.radius.degree - r2.radius.degree) < 1e-10

    def test_ellipse_roundtrip(self):
        """Test ellipse region round-trip conversion."""
        original_stcs = "Ellipse ICRS BARYCENTER 180.0 10.0 0.5 0.3 45.0"
        regions = _parse_stcs(original_stcs)
        serialized = _serialize_stcs(regions)
        regions2 = _parse_stcs(serialized)

        assert len(regions) == len(regions2) == 1
        r1, r2 = regions[0], regions2[0]

        assert isinstance(r1, EllipseSkyRegion)
        assert isinstance(r2, EllipseSkyRegion)
        assert abs(r1.center.ra.degree - r2.center.ra.degree) < 1e-10
        assert abs(r1.center.dec.degree - r2.center.dec.degree) < 1e-10
        assert abs(r1.width.degree - r2.width.degree) < 1e-10
        assert abs(r1.height.degree - r2.height.degree) < 1e-10
        assert abs(r1.angle.degree - r2.angle.degree) < 1e-10

    def test_pixel_roundtrip(self):
        """Test pixel region round-trip conversion."""
        original_stcs = "Circle IMAGE UNKNOWN 100.5 200.3 15.0"
        regions = _parse_stcs(original_stcs)
        serialized = _serialize_stcs(regions)
        regions2 = _parse_stcs(serialized)

        assert len(regions) == len(regions2) == 1
        r1, r2 = regions[0], regions2[0]

        assert isinstance(r1, CirclePixelRegion)
        assert isinstance(r2, CirclePixelRegion)
        assert abs(r1.center.x - r2.center.x) < 1e-10
        assert abs(r1.center.y - r2.center.y) < 1e-10
        assert abs(r1.radius - r2.radius) < 1e-10

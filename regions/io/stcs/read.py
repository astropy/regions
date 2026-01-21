# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_readable_fileobj
from astropy.utils.exceptions import AstropyUserWarning

from regions.core import PixCoord, RegionMeta, Regions, RegionVisual
from regions.core.registry import RegionsRegistry
from regions.io.stcs.core import (STCS_PATTERNS, STCSParserError,
                                  parse_coordinate_frame, parse_numbers,
                                  parse_unit, stcs_frame_map,
                                  stcs_shape_to_region, validate_stcs_string)

__all__ = []


@RegionsRegistry.register(Regions, 'read', 'stcs')
def _read_stcs(filename, cache=False):
    """
    Read an STC-S region file as a list of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The filename of the file to access.

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : `regions.Regions`
        A `Regions` object containing a list of `~regions.Region`
        objects.
    """
    with get_readable_fileobj(filename, cache=cache) as fh:
        region_string = fh.read()
        return _parse_stcs(region_string)


@RegionsRegistry.register(Regions, 'parse', 'stcs')
def _parse_stcs(region_str):
    """
    Parse an STC-S region string to `~regions.Region` objects.

    Parameters
    ----------
    region_str : str
        STC-S region string.

    Returns
    -------
    regions : `regions.Regions`
        A `Regions` object containing a list of `~regions.Region`
        objects.
    """
    if not validate_stcs_string(region_str):
        raise STCSParserError("Invalid STC-S string format")

    regions = []

    # Split the string into lines and process each line
    lines = [line.strip() for line in region_str.split('\n') if line.strip()]

    for line in lines:
        # Skip comments
        if line.startswith('#'):
            continue

        try:
            region = _parse_stcs_line(line)
            if region is not None:
                regions.append(region)
        except Exception as e:
            warnings.warn(f"Failed to parse STC-S line '{line}': {e}",
                         AstropyUserWarning)
            continue

    return Regions(regions)


def _parse_stcs_line(line):
    """
    Parse a single STC-S line into a Region object.

    Parameters
    ----------
    line : str
        A single STC-S line.

    Returns
    -------
    region : `~regions.Region` or None
        Parsed region object or None if parsing failed.
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    # Parse coordinate frame and reference position
    frame, refpos = parse_coordinate_frame(line)
    unit = parse_unit(line, frame)

    # Try to match different shape patterns
    for shape_name, pattern in STCS_PATTERNS.items():
        if shape_name in ['frame', 'refpos', 'unit', 'number']:
            continue

        match = pattern.search(line)
        if match:
            try:
                return _parse_shape(shape_name, match, frame, refpos, unit)
            except Exception as e:
                raise STCSParserError(f"Failed to parse {shape_name}: {e}")

    raise STCSParserError(f"No recognized shape pattern found in: {line}")


def _parse_shape(shape_name, match, frame, refpos, unit):
    """
    Parse a specific shape from regex match.

    Parameters
    ----------
    shape_name : str
        Name of the shape (circle, ellipse, etc.).
    match : re.Match
        Regex match object.
    frame : str
        Coordinate frame.
    refpos : str
        Reference position.
    unit : str
        Coordinate unit.

    Returns
    -------
    region : `~regions.Region`
        Parsed region object.
    """
    # Extract coordinate and parameter data from the match
    coord_frame_str = match.group(1).strip()
    numbers_str = match.group(2).strip()

    # Parse numbers from the string
    numbers = parse_numbers(numbers_str)

    # Create metadata - only use valid RegionMeta keys
    meta = RegionMeta()
    meta['frame'] = frame  # Use the standard 'frame' key
    # Store STC-S specific info in comment for round-trip compatibility
    meta['comment'] = f'stcs_refpos={refpos} stcs_unit={unit}'

    # Determine if this is pixel or sky coordinates
    is_pixel = (frame.upper() == 'IMAGE')
    coord_type = 'pixel' if is_pixel else 'sky'

    # Get the appropriate region class (convert shape name to title case for mapping)
    shape_title = shape_name.title()
    if shape_title not in stcs_shape_to_region[coord_type]:
        raise STCSParserError(f"Unsupported shape: {shape_name}")

    region_class = stcs_shape_to_region[coord_type][shape_title]

    # Parse based on shape type
    if shape_name == 'circle':
        return _parse_circle(region_class, numbers, frame, unit, meta, is_pixel)
    elif shape_name == 'ellipse':
        return _parse_ellipse(region_class, numbers, frame, unit, meta, is_pixel)
    elif shape_name == 'box':
        return _parse_box(region_class, numbers, frame, unit, meta, is_pixel)
    elif shape_name == 'polygon':
        return _parse_polygon(region_class, numbers, frame, unit, meta, is_pixel)
    elif shape_name == 'position':
        return _parse_position(region_class, numbers, frame, unit, meta, is_pixel)
    else:
        raise STCSParserError(f"Parser not implemented for shape: {shape_name}")


def _parse_circle(region_class, numbers, frame, unit, meta, is_pixel):
    """Parse Circle shape."""
    if len(numbers) < 3:
        raise STCSParserError("Circle requires at least 3 parameters: lon, lat, radius")

    lon, lat, radius = numbers[0], numbers[1], numbers[2]

    if is_pixel:
        center = PixCoord(lon, lat)
        region = region_class(center=center, radius=radius)
    else:
        astropy_frame = stcs_frame_map.get(frame, 'icrs')
        center = SkyCoord(lon * u.Unit(unit), lat * u.Unit(unit), frame=astropy_frame)
        region = region_class(center=center, radius=radius * u.Unit(unit))

    # Set metadata after creation
    region.meta = meta
    return region


def _parse_ellipse(region_class, numbers, frame, unit, meta, is_pixel):
    """Parse Ellipse shape."""
    if len(numbers) < 5:
        raise STCSParserError("Ellipse requires at least 5 parameters: lon, lat, semi_major, semi_minor, angle")

    lon, lat, semi_major, semi_minor, angle = numbers[0], numbers[1], numbers[2], numbers[3], numbers[4]

    if is_pixel:
        center = PixCoord(lon, lat)
        region = region_class(center=center, width=2*semi_major, height=2*semi_minor,
                            angle=angle * u.degree)
    else:
        astropy_frame = stcs_frame_map.get(frame, 'icrs')
        center = SkyCoord(lon * u.Unit(unit), lat * u.Unit(unit), frame=astropy_frame)
        region = region_class(center=center, width=2*semi_major * u.Unit(unit),
                            height=2*semi_minor * u.Unit(unit),
                            angle=angle * u.degree)

    # Set metadata after creation
    region.meta = meta
    return region


def _parse_box(region_class, numbers, frame, unit, meta, is_pixel):
    """Parse Box shape."""
    if len(numbers) < 5:
        raise STCSParserError("Box requires at least 5 parameters: lon, lat, width, height, angle")

    lon, lat, width, height, angle = numbers[0], numbers[1], numbers[2], numbers[3], numbers[4]

    if is_pixel:
        center = PixCoord(lon, lat)
        region = region_class(center=center, width=width, height=height,
                            angle=angle * u.degree)
    else:
        astropy_frame = stcs_frame_map.get(frame, 'icrs')
        center = SkyCoord(lon * u.Unit(unit), lat * u.Unit(unit), frame=astropy_frame)
        region = region_class(center=center, width=width * u.Unit(unit),
                            height=height * u.Unit(unit),
                            angle=angle * u.degree)

    # Set metadata after creation
    region.meta = meta
    return region


def _parse_polygon(region_class, numbers, frame, unit, meta, is_pixel):
    """Parse Polygon shape."""
    if len(numbers) < 6 or len(numbers) % 2 != 0:
        raise STCSParserError("Polygon requires an even number of coordinates (at least 6)")

    if is_pixel:
        vertices = PixCoord([numbers[i] for i in range(0, len(numbers), 2)],
                           [numbers[i] for i in range(1, len(numbers), 2)])
        region = region_class(vertices=vertices)
    else:
        astropy_frame = stcs_frame_map.get(frame, 'icrs')
        lon_coords = [numbers[i] * u.Unit(unit) for i in range(0, len(numbers), 2)]
        lat_coords = [numbers[i] * u.Unit(unit) for i in range(1, len(numbers), 2)]
        vertices = SkyCoord(lon_coords, lat_coords, frame=astropy_frame)
        region = region_class(vertices=vertices)

    # Set metadata after creation
    region.meta = meta
    return region


def _parse_position(region_class, numbers, frame, unit, meta, is_pixel):
    """Parse Position shape."""
    if len(numbers) < 2:
        raise STCSParserError("Position requires at least 2 parameters: lon, lat")

    lon, lat = numbers[0], numbers[1]

    if is_pixel:
        center = PixCoord(lon, lat)
        region = region_class(center=center)
    else:
        astropy_frame = stcs_frame_map.get(frame, 'icrs')
        center = SkyCoord(lon * u.Unit(unit), lat * u.Unit(unit), frame=astropy_frame)
        region = region_class(center=center)

    # Set metadata after creation
    region.meta = meta
    return region

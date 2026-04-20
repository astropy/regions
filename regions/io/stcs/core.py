# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools
import re

from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
                            EllipsePixelRegion, EllipseSkyRegion,
                            PolygonPixelRegion, PolygonSkyRegion,
                            PointPixelRegion, PointSkyRegion,
                            RectanglePixelRegion, RectangleSkyRegion)

__all__ = []


# STC-S coordinate frame mappings to astropy coordinate frames
stcs_frame_map = {
    'ICRS': 'icrs',
    'FK5': 'fk5',
    'FK4': 'fk4',
    'GALACTIC': 'galactic',
    'ECLIPTIC': 'barycentricmeanecliptic',
    'IMAGE': 'image',
    'GEOCENTER': 'icrs',  # Default for geocentric
    'BARYCENTER': 'icrs', # Default for barycentric
    'TOPOCENTER': 'icrs', # Default for topocentric
}

# Reference position mappings
stcs_refpos_map = {
    'GEOCENTER': 'geocentric',
    'BARYCENTER': 'barycentric',
    'TOPOCENTER': 'topocentric',
    'HELIOCENTER': 'heliocentric',
    'LSR': 'lsr',
    'LSRK': 'lsrk',
    'LSRD': 'lsrd',
    'GALACTIC_CENTER': 'galactic_center',
    'LOCAL_GROUP_CENTER': 'local_group_center',
    'UNKNOWN': 'unknown',
    'RELOCATABLE': 'relocatable',
}

# STC-S shape mappings to region classes
stcs_shape_to_region = {
    'pixel': {
        'Circle': CirclePixelRegion,
        'Ellipse': EllipsePixelRegion,
        'Box': RectanglePixelRegion,
        'Polygon': PolygonPixelRegion,
        'Position': PointPixelRegion,
    },
    'sky': {
        'Circle': CircleSkyRegion,
        'Ellipse': EllipseSkyRegion,
        'Box': RectangleSkyRegion,
        'Polygon': PolygonSkyRegion,
        'Position': PointSkyRegion,
    }
}

# Region class to STC-S shape mappings
region_to_stcs_shape = {
    'CirclePixelRegion': 'Circle',
    'CircleSkyRegion': 'Circle',
    'EllipsePixelRegion': 'Ellipse',
    'EllipseSkyRegion': 'Ellipse',
    'RectanglePixelRegion': 'Box',
    'RectangleSkyRegion': 'Box',
    'PolygonPixelRegion': 'Polygon',
    'PolygonSkyRegion': 'Polygon',
    'PointPixelRegion': 'Position',
    'PointSkyRegion': 'Position',
}

# STC-S unit mappings
stcs_unit_map = {
    'deg': 'degree',
    'degree': 'degree',
    'degrees': 'degree',
    'rad': 'radian',
    'radian': 'radian',
    'radians': 'radian',
    'arcmin': 'arcmin',
    'arcsec': 'arcsec',
    'mas': 'mas',
    'pixel': 'pixel',
    'pix': 'pixel',
}

# Default units for different coordinate systems
default_units = {
    'ICRS': 'degree',
    'FK5': 'degree',
    'FK4': 'degree',
    'GALACTIC': 'degree',
    'ECLIPTIC': 'degree',
    'IMAGE': 'pixel',
}

# Regular expressions for parsing STC-S strings
STCS_PATTERNS = {
    # Basic shape patterns
    'circle': re.compile(r'Circle\s+([\w\s]+?)\s+([\d\.\s\-\+e]+)', re.IGNORECASE),
    'ellipse': re.compile(r'Ellipse\s+([\w\s]+?)\s+([\d\.\s\-\+e]+)', re.IGNORECASE),
    'box': re.compile(r'Box\s+([\w\s]+?)\s+([\d\.\s\-\+e]+)', re.IGNORECASE),
    'polygon': re.compile(r'Polygon\s+([\w\s]+?)\s+([\d\.\s\-\+e]+)', re.IGNORECASE),
    'position': re.compile(r'Position\s+([\w\s]+?)\s+([\d\.\s\-\+e]+)', re.IGNORECASE),

    # Coordinate system patterns
    'frame': re.compile(r'\b(ICRS|FK5|FK4|GALACTIC|ECLIPTIC|IMAGE)\b', re.IGNORECASE),
    'refpos': re.compile(r'\b(GEOCENTER|BARYCENTER|TOPOCENTER|HELIOCENTER|LSR|LSRK|LSRD)\b', re.IGNORECASE),

    # Unit patterns
    'unit': re.compile(r'\bunit\s+(deg|degree|degrees|rad|radian|radians|arcmin|arcsec|mas|pixel|pix)\b', re.IGNORECASE),

    # Number patterns
    'number': re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?'),
}

# Template for STC-S serialization formats
stcs_templates = {
    'Circle': 'Circle {frame} {refpos} {center_lon} {center_lat} {radius}',
    'Ellipse': 'Ellipse {frame} {refpos} {center_lon} {center_lat} {semi_major} {semi_minor} {angle}',
    'Box': 'Box {frame} {refpos} {center_lon} {center_lat} {width} {height} {angle}',
    'Polygon': 'Polygon {frame} {refpos} {vertices}',
    'Position': 'Position {frame} {refpos} {lon} {lat}',
}


class STCSParserError(Exception):
    """
    A custom exception for STC-S parsing errors.
    """
    pass


def parse_coordinate_frame(stcs_string):
    """
    Parse coordinate frame and reference position from STC-S string.

    Parameters
    ----------
    stcs_string : str
        The STC-S string to parse.

    Returns
    -------
    frame : str
        The coordinate frame (default: 'ICRS').
    refpos : str
        The reference position (default: 'UNKNOWN').
    """
    frame_match = STCS_PATTERNS['frame'].search(stcs_string)
    frame = frame_match.group(1).upper() if frame_match else 'ICRS'

    refpos_match = STCS_PATTERNS['refpos'].search(stcs_string)
    refpos = refpos_match.group(1).upper() if refpos_match else 'UNKNOWN'

    return frame, refpos


def parse_unit(stcs_string, frame='ICRS'):
    """
    Parse unit from STC-S string.

    Parameters
    ----------
    stcs_string : str
        The STC-S string to parse.
    frame : str
        The coordinate frame (for default units).

    Returns
    -------
    unit : str
        The unit string (default based on frame).
    """
    unit_match = STCS_PATTERNS['unit'].search(stcs_string)
    if unit_match:
        return stcs_unit_map.get(unit_match.group(1).lower(), 'degree')
    else:
        return default_units.get(frame, 'degree')


def parse_numbers(number_string):
    """
    Parse a string of space-separated numbers.

    Parameters
    ----------
    number_string : str
        String containing space-separated numbers.

    Returns
    -------
    numbers : list of float
        List of parsed numbers.
    """
    numbers = STCS_PATTERNS['number'].findall(number_string)
    return [float(num) for num in numbers]


def format_coordinate(value, precision=8):
    """
    Format a coordinate value for STC-S output.

    Parameters
    ----------
    value : float
        The coordinate value to format.
    precision : int
        Number of decimal places.

    Returns
    -------
    formatted : str
        Formatted coordinate string.
    """
    return f"{value:.{precision}f}".rstrip('0').rstrip('.')


def validate_stcs_string(stcs_string):
    """
    Basic validation of an STC-S string.

    Parameters
    ----------
    stcs_string : str
        The STC-S string to validate.

    Returns
    -------
    is_valid : bool
        True if the string appears to be valid STC-S.
    """
    if not isinstance(stcs_string, str):
        return False

    # Check for at least one shape keyword
    shape_keywords = ['Circle', 'Ellipse', 'Box', 'Polygon', 'Position']
    return any(keyword.lower() in stcs_string.lower() for keyword in shape_keywords)

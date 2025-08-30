#!/usr/bin/env python3
"""
Example usage of the STC-S I/O module for astropy-regions.

This script demonstrates how to read and write STC-S region files
using the astropy-regions package.
"""

import astropy.units as u
from astropy.coordinates import SkyCoord

from regions.core import PixCoord, Regions
from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
                            EllipseSkyRegion, PolygonSkyRegion,
                            PointSkyRegion)


def example_read_stcs():
    """Example of reading STC-S files."""
    print("=== Reading STC-S Files ===")

    # Read from file (this would work if the file exists)
    # regions = Regions.read('example.stcs', format='stcs')

    # Parse STC-S string directly
    stcs_string = """
    # Example STC-S regions
    Circle ICRS BARYCENTER 180.0 10.0 0.5
    Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0
    Position FK5 GEOCENTER 85.0 -15.0
    """

    regions = Regions.parse(stcs_string, format='stcs')

    print(f"Parsed {len(regions)} regions:")
    for i, region in enumerate(regions):
        print(f"  {i+1}. {region}")

    return regions


def example_write_stcs():
    """Example of writing STC-S files."""
    print("\n=== Writing STC-S Files ===")

    # Create some example regions
    regions = []

    # Sky circle
    center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
    circle = CircleSkyRegion(center=center, radius=0.5 * u.degree)
    regions.append(circle)

    # Sky ellipse
    center = SkyCoord(150.0, -20.0, unit='degree', frame='icrs')
    ellipse = EllipseSkyRegion(center=center, width=2.0 * u.degree,
                              height=1.0 * u.degree, angle=45.0 * u.degree)
    regions.append(ellipse)

    # Sky polygon
    vertices = SkyCoord([45.0, 50.0, 50.0, 45.0],
                       [45.0, 45.0, 50.0, 50.0],
                       unit='degree', frame='icrs')
    polygon = PolygonSkyRegion(vertices=vertices)
    regions.append(polygon)

    # Point region
    center = SkyCoord(85.0, -15.0, unit='degree', frame='fk5')
    point = PointSkyRegion(center=center)
    regions.append(point)

    # Pixel circle
    center = PixCoord(100.5, 200.3)
    pixel_circle = CirclePixelRegion(center=center, radius=15.0)
    regions.append(pixel_circle)

    regions_obj = Regions(regions)

    # Serialize to STC-S string
    stcs_string = regions_obj.serialize(format='stcs')
    print("Serialized STC-S:")
    print(stcs_string)

    # Write to file (uncomment to actually write)
    # regions_obj.write('output.stcs', format='stcs', overwrite=True)

    return stcs_string


def example_round_trip():
    """Example of round-trip conversion."""
    print("\n=== Round-trip Conversion ===")

    original_stcs = """Circle ICRS BARYCENTER 180.0 10.0 0.5
Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0
Position FK5 GEOCENTER 85.0 -15.0"""

    print("Original STC-S:")
    print(original_stcs)

    # Parse
    regions = Regions.parse(original_stcs, format='stcs')
    print(f"\nParsed {len(regions)} regions")

    # Serialize back
    serialized = regions.serialize(format='stcs')
    print("\nSerialized back to STC-S:")
    print(serialized)

    # Parse again to verify
    regions2 = Regions.parse(serialized, format='stcs')
    print(f"\nRe-parsed {len(regions2)} regions - round-trip successful!")


if __name__ == '__main__':
    print("STC-S I/O Module Examples")
    print("=" * 40)

    try:
        example_read_stcs()
        example_write_stcs()
        example_round_trip()

        print("\n" + "=" * 40)
        print("All examples completed successfully!")

    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()

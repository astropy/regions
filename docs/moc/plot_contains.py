import numpy as np

from astropy.utils.data import get_pkg_data_filename
from regions import MOCSkyRegion, WCS
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
# Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
filename = get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions')
moc = MOCSkyRegion.from_fits(filename)

# Generate random sky coordinates
ra = np.random.randint(low=0, high=360, size=1000) * u.deg
dec = np.random.randint(low=-90, high=90, size=1000) * u.deg
# Get the mask of the sky coordinates contained in the MOCSkyRegion instance
inside_mask = moc.contains(ra, dec)

coords_inside = SkyCoord(ra[inside_mask], dec[inside_mask])
coords_outside = SkyCoord(ra[~inside_mask], dec[~inside_mask])

# Plot the MOC using matplotlib
import matplotlib.pyplot as plt
fig = plt.figure(111, figsize=(10, 10))
# Define a WCS as a context
with WCS(fig, 
        fov=100 * u.deg,
        center=SkyCoord(0, 0, unit='deg', frame='icrs'),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Cast the MOC sky region to a MOC pixel region
    reg = moc.to_pixel(wcs)
    # Call the plot method of MOC pixel region
    reg.plot(ax=ax, alpha=0.5, fill=True, linewidth=0, color='r')

    # Project the sky coordinates to the image system
    from astropy.wcs.utils import skycoord_to_pixel
    x_in, y_in = skycoord_to_pixel(coords_inside, wcs=wcs)
    x_out, y_out = skycoord_to_pixel(coords_outside, wcs=wcs)

    plt.scatter(x_out, y_out, s=16, c='black', alpha=0.5, marker='^', zorder=2, label='outside')
    plt.scatter(x_in, y_in, s=16, c='green', alpha=0.5, marker='^', zorder=3, label='inside')

plt.legend()
plt.title("MOCSkyRegion of GALEX")
plt.xlabel('ra')
plt.ylabel('dec')
plt.grid(color="black", linestyle="dotted")
plt.show()
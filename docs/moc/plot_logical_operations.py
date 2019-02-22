from astropy.utils.data import get_pkg_data_filename
from regions import MOCSkyRegion, WCS
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
# Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))
sdss9 = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-SDSS9-r.fits', package='regions'))

# Compute the union and intersection
intersection = galex.intersection(sdss9)
union = galex.union(sdss9)
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
    union.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, linewidth=0, color='r', label="Union")
    union.border(ax=ax, wcs=wcs, color='r')
    intersection.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, linewidth=0, color='g', label="Intersection")
    intersection.border(ax=ax, wcs=wcs, color='g')
    ax.legend()

plt.title("Logical operations between GALEX and SDSS9")
plt.xlabel('ra')
plt.ylabel('dec')
plt.grid(color="black", linestyle="dotted")
plt.show()
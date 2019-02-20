# Example illustrating how to:
# - Load a MOC from a FITS file
# - Plot a MOC defining a WCS

from astropy.utils.data import get_pkg_data_filename
from regions import MOCSkyRegion, WCS
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
# Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
filename = get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions')
moc = MOCSkyRegion.from_fits(filename)
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
    reg.plot(ax=ax, alpha=0.5, fill=True, linewidth=0.5, color='r')

plt.title("MOCSkyRegion of GALEX")
plt.xlabel('ra')
plt.ylabel('dec')
plt.grid(color="black", linestyle="dotted")
plt.show()
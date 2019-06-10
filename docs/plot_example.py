"""Example how to plot sky regions on a sky image.
"""
from astropy.coordinates import SkyCoord, Angle
from regions import (
    make_example_dataset,
    CircleSkyRegion,
)
import matplotlib.pyplot as plt

config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
dataset = make_example_dataset(data='simulated', config=config)
wcs = dataset.wcs

fig = plt.figure()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)

ax.set_xlim(-0.5, dataset.config['shape'][1] - 0.5)
ax.set_ylim(-0.5, dataset.config['shape'][0] - 0.5)

ax.imshow(dataset.image.data, cmap='gray', vmin=0, vmax=1,
          interpolation='nearest', origin='lower')

for source in dataset.source_table:
    # Plot a sky circle around each source
    center = SkyCoord(source['GLON'], source['GLAT'], unit='deg', frame='galactic')
    radius = Angle(20, 'deg')
    region = CircleSkyRegion(center=center, radius=radius)
    pix_region = region.to_pixel(wcs=wcs)

    pix_region.plot(ax=ax, edgecolor='yellow', facecolor='yellow', alpha=0.5, lw=3)

"""
Example script illustrating compound regions.
"""
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion, make_example_dataset

# load example dataset to get skymap
config = dict(crval=(0, 0),
              crpix=(180, 90),
              cdelt=(-1, 1),
              shape=(180, 360))

dataset = make_example_dataset(data='simulated', config=config)
wcs = dataset.wcs

# remove sources
dataset.image.data = np.zeros_like(dataset.image.data)

# define 2 sky circles
circle1 = CircleSkyRegion(
    center=SkyCoord(20, 0, unit='deg', frame='galactic'),
    radius=Angle('30 deg')
)

circle2 = CircleSkyRegion(
    center=SkyCoord(50, 45, unit='deg', frame='galactic'),
    radius=Angle('30 deg'),
)

# define skycoords
lon = np.arange(-180, 181, 10)
lat = np.arange(-90, 91, 10)
coords = np.array(np.meshgrid(lon, lat)).T.reshape(-1, 2)
skycoords = SkyCoord(coords, unit='deg', frame='galactic')

# get events in AND and XOR
compound_and = circle1 & circle2
compound_xor = circle1 ^ circle2

mask_and = compound_and.contains(skycoords, wcs)
skycoords_and = skycoords[mask_and]
mask_xor = compound_xor.contains(skycoords, wcs)
skycoords_xor = skycoords[mask_xor]

# plot
fig = plt.figure()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs, aspect='equal')

ax.scatter(skycoords.l.value, skycoords.b.value, label='all',
           transform=ax.get_transform('galactic'))
ax.scatter(skycoords_xor.l.value, skycoords_xor.b.value, color='orange',
           label='xor', transform=ax.get_transform('galactic'))
ax.scatter(skycoords_and.l.value, skycoords_and.b.value, color='magenta',
           label='and', transform=ax.get_transform('galactic'))

circle1.to_pixel(wcs=wcs).plot(ax=ax, edgecolor='green', facecolor='none', alpha=0.8, lw=3)
circle2.to_pixel(wcs=wcs).plot(ax=ax, edgecolor='red', facecolor='none', alpha=0.8, lw=3)

ax.legend(loc='lower right')

ax.set_xlim(-0.5, dataset.config['shape'][1] - 0.5)
ax.set_ylim(-0.5, dataset.config['shape'][0] - 0.5)

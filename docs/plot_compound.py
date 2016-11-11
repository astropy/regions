import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion, make_example_dataset

# load example dataset to get skymap
config = dict(crval=(0, 0),
              crpix=(8, 8),
              cdelt=(-1, 1),
              shape=(16, 16))

dataset = make_example_dataset(data='simulated', config=config)
wcs = dataset.wcs

# remove sources
dataset.image.data = np.zeros_like(dataset.image.data)

# define 2 sky circles
circle1 = CircleSkyRegion(
    center=SkyCoord(1,2, unit='deg', frame='galactic'),
    radius=Angle('5 deg')
)

circle2 = CircleSkyRegion(
    center=SkyCoord(-4,3, unit='deg', frame='galactic'),
    radius=Angle('3 deg'),
)


# define skycoords
lon = np.concatenate([np.arange(352, 360), np.arange(0,8)])
lat = np.arange(-7,7)
coords = np.array(np.meshgrid(lon, lat)).T.reshape(-1,2)
skycoords = SkyCoord(coords, unit='deg', frame='galactic')


# get events in AND and XOR
compound_and = circle1 & circle2
compound_xor = circle1 ^ circle2

mask_and = compound_and.contains(skycoords)
skycoords_and = skycoords[mask_and]
mask_xor = compound_xor.contains(skycoords)
skycoords_xor = skycoords[mask_xor]

# plot
fig = plt.figure()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)

ax.imshow(dataset.image.data, cmap='gray', vmin=0, vmax=1,
          interpolation='nearest', origin='lower')

ax.scatter(skycoords.l, skycoords.b, label='all events',
           transform=ax.get_transform('galactic'))
ax.scatter(skycoords_xor.l, skycoords_xor.b, color='orange', label='XOR',
          transform=ax.get_transform('galactic'))
ax.scatter(skycoords_and.l, skycoords_and.b, color='magenta', label='AND',
          transform=ax.get_transform('galactic'))

circle1.plot(ax, edgecolor='yellow', facecolor='none', alpha=0.8, lw=3)
circle2.plot(ax, edgecolor='cyan', facecolor='none', alpha=0.8, lw=3)

ax.legend(loc='lower right')

ax.set_xlim(-0.5, dataset.config['shape'][1] - 0.5)
ax.set_ylim(-0.5, dataset.config['shape'][0] - 0.5)

plt.show()

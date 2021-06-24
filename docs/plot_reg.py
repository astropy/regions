# Licensed under a 3-clause BSD style license - see LICENSE.rst

from matplotlib import pyplot as plt

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

from regions import Regions


image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
print(image_file)
image_data = fits.getdata(image_file, ext=0, memmap=False)

ax = plt.gca()
plt.imshow(image_data, cmap='gray')

print(ax.get_xlim(), ax.get_ylim())
ax.set_ylim([-0.5, 892.5])
regs = Regions.read(get_pkg_data_filename('data/plot_image.reg',
                                          package='regions.io.ds9.tests'),
                    format='ds9')

for i, reg in enumerate(regs):
    reg.plot(ax=ax)

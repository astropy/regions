import os

from matplotlib import pyplot as plt

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

from regions import read_ds9

image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
print(image_file)
image_data = fits.getdata(image_file, ext=0)

ax = plt.gca()
plt.imshow(image_data, cmap='gray')

print(ax.get_xlim(), ax.get_ylim())
ax.set_ylim([-0.5, 892.5])
regs = read_ds9(os.path.join(os.path.dirname(__file__), 'plot_image.reg'))

for i, reg in enumerate(regs):
    reg.plot(ax=ax)

plt.show()

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy import units as u


from matplotlib import pyplot as plt

from regions import *

image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
print(image_file)
image_data = fits.getdata(image_file, ext=0)

ax = plt.gca()
plt.imshow(image_data, cmap='gray')

regs = []
print(ax.get_xlim(), ax.get_ylim())
ax.set_ylim([-0.5, 892.5])
colors = ['red', 'blue', 'orange', 'brown', 'yellow', 'violet']

regs.append(CirclePixelRegion(PixCoord(500, 600), 50))
# regs.append(EllipseAnnulusPixelRegion(PixCoord(200, 700), 50, 60, 80, 100, 45 * u.deg ))
regs.append(LinePixelRegion(PixCoord(400, 200), PixCoord(500, 100)))
# regs.append(RectangleAnnulusPixelRegion(PixCoord(300, 400), 50, 60, 80, 100, -45 * u.deg))
regs.append(PointPixelRegion(PixCoord(100, 800), visual=RegionVisual(symbol='+')))
regs.append(PointPixelRegion(PixCoord(100, 100), visual=RegionVisual(symbol='x')))
regs.append(PointPixelRegion(PixCoord(100, 300), visual=RegionVisual(symbol='*')))
regs.append(PointPixelRegion(PixCoord(100, 400), visual=RegionVisual(symbol='x')))
regs.append(PointPixelRegion(PixCoord(100, 200)))
regs.append(EllipsePixelRegion(PixCoord(300, 750), 90, 60, 75 * u.deg))
regs.append(TextPixelRegion(PixCoord(150, 550), "Text", visual=RegionVisual(textangle='45',font='algeria')))
regs.append(CircleAnnulusPixelRegion(PixCoord(650, 300), 60, 90))
regs.append(PolygonPixelRegion(PixCoord([450, 450, 550, 600], [750, 700, 650, 750])))
regs.append(RectanglePixelRegion(PixCoord(400, 400), 100, 80))

for i, reg in enumerate(regs):
    reg.visual['color'] = colors[i%6]
    reg.plot(ax)

plt.show()

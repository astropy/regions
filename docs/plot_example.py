"""Make an example plot for regions.
"""
from regions import make_example_dataset
import matplotlib.pyplot as plt

config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
dataset = make_example_dataset(data='simulated', config=config)
print(dataset.image.data.max())

plt.imshow(dataset.image.data, cmap='gray', vmin=0, vmax=1, interpolation='none', origin='bottom')

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import numpy as np

fig = plt.figure(figsize=(8, 4))
gs = gridspec.GridSpec(1, 2, width_ratios=[20, 1])

ax_img = fig.add_subplot(gs[0])
im = ax_img.imshow(np.random.rand(10, 10), cmap='viridis')

# Manually place a tight colorbar axis
cax = fig.add_axes([ax_img.get_position().x1 + 0.01,  # left
                    ax_img.get_position().y0,         # bottom
                    0.015,                            # width
                    ax_img.get_position().height])    # height

fig.colorbar(im, cax=cax)
plt.show()

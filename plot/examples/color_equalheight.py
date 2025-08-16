# ----------------------------------------------------------------------------
# Title:   Scientific Visualisation - Python & Matplotlib
# Author:  Nicolas P. Rougier
# License: BSD
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

p = plt.rcParams
p["figure.figsize"] = 7, 7
# p["font.sans-serif"] = ["Roboto Condensed"]
p["font.weight"] = "light"
p["ytick.minor.visible"] = True
p["xtick.minor.visible"] = True

# Some data
n = 64
X, Z = np.meshgrid(
    np.linspace(-0.5 + 0.5 / n, +0.5 - 0.5 / n, n),
    np.linspace(-0.5 + 0.5 / n, +0.5 - 0.5 / n, n),
)
Y = 0.75 * np.exp(-10 * (X ** 2 + Z ** 2))


def f(x, y):
    return (1 - x / 2 + x ** 5 + y ** 3) * np.exp(-(x ** 2) - y ** 2)


x = np.linspace(-3, 3, 100)
y = np.linspace(-3, 3, 100)
X, Y = np.meshgrid(x, y)
Z = 0.5 * f(X, Y)

fig = plt.figure()
nrows, ncols = 1, 1
gspec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=fig)

ax = plt.subplot(gspec[0, 0], aspect=1)
ax.set_xlim(0, 1)
ax.set_xticks(np.linspace(0, 1, 4 + 1))
ax.set_xlabel("X Label")
ax.set_ylim(0, 1)
ax.set_yticks(np.linspace(0, 1, 4 + 1))
ax.set_ylabel("Y Label")
ax.set_title("Title", weight=500)
I = ax.imshow(Z, extent=[0, 1, 0, 1])

# Get position of imshow Axes and place colorbar tightly next to it
bbox = ax.get_position()
cax = fig.add_axes([
    bbox.x1 + 0.01,   # left
    bbox.y0,          # bottom
    0.02,             # width
    bbox.height       # height
])
fig.colorbar(I, cax=cax)

plt.show()

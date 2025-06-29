# ----------------------------------------------------------------------------
# Title:   Scientific Visualisation - Python & Matplotlib
# Author:  Nicolas P. Rougier
# License: BSD
# ----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


fig = plt.figure(figsize=(5.5, 4))
ax = plt.subplot(1, 1, 1, xlim=[-1, +1], xticks=[], ylim=[-1, +1], yticks=[], aspect=1)

np.random.seed(123)
X = np.random.normal(0, 0.35, 1000)
Y = np.random.normal(0, 0.35, 1000)

ax.scatter(X, Y, edgecolor="None", facecolor="C1", alpha=0.5)

I = np.random.choice(len(X), size=5, replace=False)
Px, Py = X[I], Y[I]
I = np.argsort(Y[I])[::-1]
Px, Py = Px[I], Py[I]

ax.scatter(Px, Py, edgecolor="black", facecolor="white", zorder=20)
ax.scatter(Px, Py, edgecolor="None", facecolor="C1", alpha=0.5, zorder=30)

y, dy = 0.25, 0.125
style = "arc,angleA=-0,angleB=0,armA=-100,armB=0,rad=0"
for i in range(len(I)):
    text = ax.annotate(
        "Point " + chr(ord("A") + i),
        xy=(Px[i], Py[i]),
        xycoords="data",
        xytext=(1.25, y - i * dy),
        textcoords="data",
        arrowprops=dict(
            arrowstyle="->",
            color="black",
            linewidth=0.75,
            shrinkA=20,
            shrinkB=5,
            patchA=None,
            patchB=None,
            connectionstyle=style,
        ),
    )
    text.arrow_patch.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )


plt.tight_layout()
plt.show()

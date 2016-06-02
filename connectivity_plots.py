import numpy as np
import matplotlib.pyplot as plt

import ebola_small_network as eb

mats = []
for p in np.logspace(-1, 0, 4):
    s = eb.SmallWorld(20, 1, p)
    connect = np.empty((20, 20))
    for i in range(20):
        for j in range(20):
            connect[i,j] = len(s.lattice[(i,j)].getNearestNeighbors())
    mats += [connect]
         
         
# Plot each slice as an independent subplot
fig, axes = plt.subplots(nrows=2, ncols=2)
for mat, ax in zip(mats, axes.flat):
    # The vmin and vmax arguments specify the color limits
    im = ax.imshow(mat, vmin=0, vmax=8)

# Make an axis for the colorbar on the right side
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
fig.colorbar(im, cax=cax)

plt.show()
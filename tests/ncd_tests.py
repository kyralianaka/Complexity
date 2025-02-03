# Testing functions in the ncd subpackage
from complexity.ncd.NCD import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

arr = np.random.binomial(1, 0.5, size=(16, 2000))

# NCD & NCD_clusters test
mat = NCD(arr, "lz")
hm = NCD_clusters(mat)

fig, ax = plt.subplots()
im = ax.matshow(hm, cmap="GnBu")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
plt.colorbar(im, cax=cax, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.setp(ax, xticks=[], yticks=[])
plt.show()

# NCD_pairwise test
x = np.random.binomial(1, 0.5, 1000)
y = np.random.binomial(1, 0.5, 1000)

pairwise_dists = dict.fromkeys([gzip, bz2, "gzip_padded", snappy, "ppm", "ppmc", "lz"])
for i in list(pairwise_dists.keys()):
    print(NCD_pairwise(x, y, i))

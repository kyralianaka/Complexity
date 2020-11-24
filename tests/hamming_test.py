from complexity.Hamming import Hamming
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#Binomial random array that should have large pairwise hamming distances
rand_sarr = np.random.binomial(1,0.5,size=(20,100))
rand_hm = Hamming(rand_sarr)

#Similar sequences that should have small pairwise hamming distances
s1 = np.tile([1,0,0,1,0,0,0,1,0,0],(10,10))
s2 = np.tile([1,1,0,0,1,0,0,0,1,0],(10,10))
sim_sarr = np.vstack((s1,s2))
sim_hm = Hamming(sim_sarr)

#Plot them 
fig, axs = plt.subplots(1,2)
axs[0].imshow(rand_hm,cmap='binary',vmin=0,vmax=1)
axs[1].imshow(sim_hm,cmap='binary',vmin=0,vmax=1)
plt.colorbar(cm.ScalarMappable(cmap='binary'),ax=axs.ravel().tolist())

#A little formatting
plt.setp(axs,xticks=[],yticks=[])
axs[0].title.set_text('Random Sequences')
axs[1].title.set_text('Similar Sequences')

plt.show()
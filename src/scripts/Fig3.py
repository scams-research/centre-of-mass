import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as matplotlib
from Variables import font_size
import paths

k = np.loadtxt(paths.data / 'sampling.txt')

# Plotting
matplotlib.rcParams.update({'font.size': font_size})
fig, axis = plt.subplots(figsize = (4.03, 4.03))

num_bins = 100
credible_intervals = [[16, 84], [2.5, 97.5], [0.15, 99.85]]

axis.set_aspect('equal')
axis.scatter(k[:, 6], k[:, 1], alpha = 0.4, color = '#0173B2', rasterized = True )
axis.plot(k[:, 6], k[:, 3], color = 'k')
axis.set_xlabel('Intrinsic mean')
axis.set_ylabel('Centre of mass')
# axis.set_xlim([0,1])
# axis.set_ylim([0,1])

plt.savefig(paths.figures / 'Fig3.pdf', bbox_inches='tight', pad_inches=0)
plt.close()

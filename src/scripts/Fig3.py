import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as matplotlib
from Variables import font_size
import paths

k = np.loadtxt(paths.data / 'sampling.txt')

# Binning Bai and Breen error for visulation
N_bins = 100
bins = np.logspace(-1, 2, N_bins)
binned_indices = np.digitize(k[:, 4], bins)

error_bb_array = np.empty([len(bins), 2000000])
error_bb_array[:] = np.nan
err = (k[:, 1] - k[:, 2]) / k[:, 5]
for z in range(0, N_bins):
    length = len(np.where(binned_indices == z)[0])
    if length != 0:
        error_bb_array[z, 0:length] = err[np.where(binned_indices == z)[0]]

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
axis.set_xlim([0,1])
axis.set_ylim([0,1])


plt.savefig(paths.figures / 'Fig3.pdf', bbox_inches='tight', pad_inches=0)
plt.close()

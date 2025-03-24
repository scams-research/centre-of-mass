import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import uniform
import matplotlib as matplotlib
from Functions import pi_centrer
from Variables import font_size
import paths

matplotlib.rcParams.update({'font.size': font_size})

k = np.loadtxt(paths.data / 'sampling.txt')
x_max = 1
x_min = 0

# Removing errors from incorrect boundary wrappings
k = k[np.where((k[:,6] > 0.02) & (k[:,6] < 0.98))]

# Binning Bai and Breen error for visulation 
N_bins = 100
bins = np.logspace(-1, 2, N_bins)
binned_indices = np.digitize(k[:, 4], bins)

error_bb_array = np.empty([len(bins),2000000])
error_bb_array[:] = np.nan
err = (k[:, 1] - k[:, 6])/k[:,5]
for z in range(0,N_bins):
    length = len(np.where(binned_indices == z)[0])
    if length != 0:
        error_bb_array[z,0:length] = err[np.where(binned_indices == z)[0]]

credible_intervals = [[16, 84], [2.5, 97.5], [0.15, 99.85]]
alpha = [0.6, 0.4, 0.2]
fig, axis = plt.subplots(figsize = (4.03, 3.3))

for i, ci in enumerate(credible_intervals):
    axis.fill_between(bins,
                    *np.nanpercentile(error_bb_array,ci, axis=1),
                    alpha=alpha[i],
                    color='#0173B2',
                    label= f'${i+1}\sigma$ - B&B method',
                    lw=0)
    

axis.set_xlabel('Asymmetry')
axis.set_ylabel("Centre of mass error / Span of particles")
axis.set_yticks([-0.1, 0., 0.1])
axis.set_xticks([0, 25, 50, 75, 100])
plt.tight_layout()

plt.savefig(paths.figures / 'Fig5.pdf', bbox_inches='tight', pad_inches=0)
plt.close()

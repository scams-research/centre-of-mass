import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import uniform
import matplotlib as matplotlib
from Functions import pi_centrer
import paths

k = np.loadtxt(paths.data / 'sampling.txt')

# Binning Bai and Breen error for visulation 
N_bins = 100
bins = np.logspace(-1, 2, N_bins)
binned_indices = np.digitize(k[:, 4], bins)

error_bb_array = np.empty([len(bins),2000000])
error_bb_array[:] = np.nan
err = (k[:, 1] - k[:, 2])/k[:,5]
for z in range(0,N_bins):
    length = len(np.where(binned_indices == z)[0])
    if length != 0:
        error_bb_array[z,0:length] = err[np.where(binned_indices == z)[0]]


credible_intervals = [[16, 84], [2.5, 97.5], [0.15, 99.85]]
alpha = [0.6, 0.4, 0.2]
matplotlib.rcParams.update({'font.size': 7})
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
plt.tight_layout()

plt.savefig(paths.figures / 'Fig5.pdf', bbox_inches='tight', pad_inches=0)
plt.close()

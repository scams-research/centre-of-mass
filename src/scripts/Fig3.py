import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import paths

k = np.loadtxt(paths.data / 'Fig3.npy')

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

credible_intervals = [[16, 84], [2.5, 97.5], [0.15, 99.85]]
alpha = [0.6, 0.4, 0.2]
fig, axis = plt.subplots(figsize=(4.03, 3.3))

axis.plot(k[:, 4], (k[:, 3] - k[:, 2]) / k[:, 5],
          label='This Work',
          color='orange',
          solid_capstyle='round')

for i, ci in enumerate(credible_intervals):
    axis.fill_between(bins,
                      *np.nanpercentile(error_bb_array, ci, axis=1),
                      alpha=alpha[i],
                      color='#0173B2',
                      label=f'${i+1}\sigma$ - B&B method',
                      lw=0)

axis.set_yticks([-0.08, -0.04, 0, 0.04, 0.08])
axis.set_xticks([0, 50, 100])
axis.set_ylim(-0.082, 0.082)
axis.set_xlabel('Asymmetry')
axis.set_ylabel('Normalised centre of mass error')
plt.tight_layout()
plt.savefig(paths.figures / 'Fig3.pdf', bbox_inches='tight', pad_inches=0)
plt.close()

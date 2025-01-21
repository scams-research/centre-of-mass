import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import uniform
import matplotlib as matplotlib
from Functions import pi_centrer

# Main
m = int(2 ** 24) 
k = np.zeros((m, 6))
x_max = 1
x_min = 0

# Variables to define wrapping
minimum_uniform = 0.65


for i in tqdm(range(m)):
    rng = np.random.default_rng(i)
    N = rng.integers(3, 512, size=1)
    n = rng.uniform(0.1, 0.45, size=1)
    particles = rng.uniform(minimum_uniform, minimum_uniform+n, size=N)
    y = particles
    particle_span = particles.max() - particles.min()
    y_pbc = y % x_max
    k[i, 0] = N
    k[i, 1] = pi_centrer(y_pbc, np.ones_like(y_pbc), x_min, x_max)
    k[i, 2] = np.average(y)
    offset = (k[i, 1] - ((x_max - x_min) / 2))
    moved = (y_pbc - offset) % (x_max - x_min)
    k[i, 3] = np.average(moved) + offset
    offset = (k[i, 1] - ((x_max - x_min) / 2))
    moved = (y - offset) % (x_max - x_min)
    # aymmetry metric
    y_hist, x_hist = np.histogram(moved, bins=np.linspace(x_min, x_max, 101), density=True)
    asymmetry = y_hist[:int(y_hist.size / 2)] - y_hist[int(y_hist.size / 2):][::-1]
    k[i, 4] = np.sum(np.abs(asymmetry))
    k[i,5] = particle_span

# Unwrapping Centre of mass if required to allow for comparison with true Centre of mass
k[:, 1][np.where(k[:, 1] < (x_max - x_min) * 0.5)] += x_max
k[:, 2][np.where(k[:, 2] < (x_max - x_min) * 0.5)] += x_max
k[:, 3][np.where(k[:, 3] < (x_max - x_min) * 0.5)] += x_max


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

axis.plot(k[:, 4], (k[:, 3] - k[:, 2])/k[:,5], label='This Work', color = 'orange',solid_capstyle='round')

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
plt.show()

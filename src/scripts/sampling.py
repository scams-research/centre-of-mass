import numpy as np
from tqdm import tqdm
from scipy.stats import uniform
from Functions import pi_centrer, intrinsic_mean
import paths

# Main
m = int(2**10)
k = np.zeros((m, 7))
x_max = 1
x_min = 0

for i in tqdm(range(m)):
    rng = np.random.default_rng(i)
    minimum_uniform = rng.uniform(0,1,size=1)
    N = rng.integers(3, 512, size=1)
    n = rng.uniform(0.1, 0.59, size=1)
    particles = rng.uniform(minimum_uniform, minimum_uniform + n, size=N)
    y = particles
    particle_span = particles.max() - particles.min()
    y_pbc = y % x_max
    k[i, 0] = N
    k[i, 1] = pi_centrer(y_pbc, np.ones_like(y_pbc), x_min, x_max)
    k[i, 2] = np.average(y)
    offset = (k[i, 1] - ((x_max - x_min) / 2))
    moved = (y_pbc - offset) % (x_max - x_min)
    k[i, 3] = (np.average(moved) + offset) % 1
    offset = (k[i, 1] - ((x_max - x_min) / 2))
    moved = (y - offset) % (x_max - x_min)
    # aymmetry metric
    y_hist, x_hist = np.histogram(moved,
                                  bins=np.linspace(x_min, x_max, 101),
                                  density=True)
    asymmetry = y_hist[:int(y_hist.size / 2)] - y_hist[int(y_hist.size /
                                                           2):][::-1]
    k[i, 4] = np.sum(np.abs(asymmetry))
    k[i, 5] = particle_span
    k[i, 6] = intrinsic_mean(y_pbc,np.ones_like(y_pbc))

# Unwrapping Centre of mass if required to allow for comparison with true Centre of mass
# k[:, 1][np.where(k[:, 1] < (x_max - x_min) * 0.5)] += x_max
# k[:, 2][np.where(k[:, 2] < (x_max - x_min) * 0.5)] += x_max
# k[:, 3][np.where(k[:, 3] < (x_max - x_min) * 0.5)] += x_max

np.savetxt(paths.data / 'sampling.txt', k)

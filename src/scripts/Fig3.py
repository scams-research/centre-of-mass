import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from Functions import pi_centrer,frechet_statistic,intrinsic_mean
from Variables import font_size
import paths


# Main
m = int(2 ** 19) 
k = np.zeros((m, 4))
x_max = 1
x_min = 0


k = np.zeros((m, 4))

for i in range(m):
    rng = np.random.default_rng(i)
    minimum_uniform = rng.uniform(0,1,size=1)
    N = rng.integers(3, 512, size=1)
    n = rng.uniform(0.1, 0.59, size=1)
    particles = rng.uniform(minimum_uniform, minimum_uniform+n, size=N)
    y = particles
    particle_span = particles.max() - particles.min()
    y_pbc = y % x_max
    k[i, 0] = N
    k[i, 1] = pi_centrer(y_pbc, np.ones_like(y_pbc), x_min, x_max)
    k[i, 2] = intrinsic_mean(y_pbc,np.ones_like(y_pbc))
    offset = (k[i, 1] - ((x_max - x_min) / 2))
    moved = (y_pbc - offset) % (x_max - x_min)
    k[i, 3] = (np.average(moved) + offset) % 1

# Plotting
matplotlib.rcParams.update({'font.size': font_size})
fig, axis = plt.subplots(figsize = (4.03, 4.03))

num_bins = 100
credible_intervals = [[16, 84], [2.5, 97.5], [0.15, 99.85]]
alpha = [0.6, 0.4, 0.2]
axis.set_aspect('equal')
axis.scatter(k[:,2],k[:,1],alpha = 0.4, color = '#0173B2', rasterized = True )
axis.plot(k[:,2],k[:,3]%1, color = 'k')
axis.set_xlabel('Intrinsic mean')
axis.set_ylabel('Centre of mass')
axis.set_xlim([0,1])
axis.set_ylim([0,1])


plt.savefig(paths.figures / 'Fig3.pdf', bbox_inches='tight', pad_inches=0)
plt.close()
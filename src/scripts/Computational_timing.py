"""
To re run the computational timing for your own machine, install the geomstats package
via pip and run this script from the top directory of the git repository.
"""

import numpy as np
import time
from geomstats.geometry.hypersphere import Hypersphere
from geomstats.learning.frechet_mean import CircleMean

from Functions import pi_centrer,geom_stats_intrinsic
from Variables import font_size
import paths

times_pi = []
ms = np.array([2**i for i in range(14, 25)])
parts = 20
time_array = np.zeros([len(ms),2])
circle = Hypersphere(dim=1)
z = CircleMean(space=circle)


times_pi = []
for m in ms:
    start = time.time()
    for i in range(m):
        rng = np.random.default_rng(i)
        minimum_uniform = rng.uniform(0, 1, size=1)
        N = parts
        n = rng.uniform(0.1, 0.59, size=1)
        particles = rng.uniform(minimum_uniform, minimum_uniform + n, size=N)
        k = pi_centrer(particles, np.ones_like(particles), 0, 1)
    end = time.time()
    times_pi.append(end - start)

time_array[:,0] = times_pi

times_frech = []
for m in ms:
    start = time.time()
    for i in range(m):
        rng = np.random.default_rng(i)
        minimum_uniform = rng.uniform(0, 1, size=1)
        N = parts
        n = rng.uniform(0.1, 0.59, size=1)
        particles = rng.uniform(minimum_uniform, minimum_uniform + n, size=N)
        x = geom_stats_intrinsic(particles,z)
    end = time.time()
    times_frech.append(end - start)

time_array[:,1] = times_frech

np.savetxt(paths.data / 'comp_time.txt',time_array)
np.savetxt(paths.data / 'ms.txt',ms)
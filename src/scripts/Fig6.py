import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Variables import font_size
import paths

matplotlib.rcParams.update({'font.size': font_size})

time_array = np.loadtxt(paths.data / 'comp_time.txt')
ms = np.loadtxt(paths.data / 'ms.txt')

fig, axis = plt.subplots(figsize = (4.03, 3.3))
axis.plot(ms / 1e6, time_array[:,0]/60, marker='.', label = 'Pseudo-centre of\nmass recentering', color = '#0173B2')
axis.plot(ms / 1e6, time_array[:,1]/60, marker='.', label = 'Intrinsic', color = 'k')
axis.set_xlabel('Number of calculations/$10^6$')
axis.set_ylabel('Time/minutes')
axis.set_xlim(0, None)
axis.set_ylim(0, None)
axis.set_xticks([0, 4, 8, 12, 16])
axis.set_xticklabels([0, 4, 8, 12, 16])
axis.set_yticks([0, 4, 8, 12])
axis.set_yticklabels([0, 4, 8, 12])
axis.legend()

plt.savefig(paths.figures / 'Fig6.pdf', bbox_inches='tight', pad_inches=0)
plt.close()
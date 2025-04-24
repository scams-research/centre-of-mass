import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Variables import font_size
import paths

matplotlib.rcParams.update({'font.size': font_size})

time_array = np.loadtxt('comp_time.txt')
ms = np.loadtxt('ms.txt')

fig, axis = plt.subplots(figsize = (4.03, 3.3))
axis.plot(ms, time_array[:,0]/60, label = 'Pseudo-recentering',color = '#0173B2')
axis.plot(ms, time_array[:,1]/60, label = 'Intrinsic', color = 'k')
axis.set_xlabel('Number of calculations')
axis.set_ylabel('Time (minutes)')
axis.legend()

plt.savefig(paths.figures / 'Fig6.pdf', bbox_inches='tight', pad_inches=0)
plt.close()
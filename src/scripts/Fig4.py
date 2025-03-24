import numpy as np
import matplotlib.gridspec as gridspec
from Functions import test_pi_centrer
from Variables import blue_set,cross_size, font_size, s_set, l_set
import matplotlib.pyplot as plt
import matplotlib
import paths

matplotlib.rcParams.update({'font.size': font_size})

# Top sub figure
coords = np.array([0.25,0.34,0.75])
x = np.arange(0, 1, 0.0001)
masses = np.array([1,1,1])

y = np.zeros(x.shape[0])
y[np.where(x==coords[0])] += 1
y[np.where(x==coords[1])] += 1
y[np.where(x==coords[2])] += 1

fft = np.fft.fft(y)
first = np.zeros_like(fft)
first[1:2] = fft[1:2]

ifft = np.fft.ifft(first).real


# Bottom sub figure
coords2 = np.array([0.25,0.5,0.75])
x2 = np.arange(0, 1, 0.0001)
masses2 = np.array([1,1,1])

y2 = np.zeros(x2.shape[0])
y2[np.where(x2==coords2[0])] += 1
y2[np.where(x2==coords2[1])] += 1
y2[np.where(x2==coords2[2])] += 1

fft2 = np.fft.fft(y2)
first2 = np.zeros_like(fft2)
first2[1:2] = fft2[1:2]

ifft2 = np.fft.ifft(first2).real


# Centre of mass calcs
True_com, Polar_com, _ = test_pi_centrer(masses,coords, box_min = 0, box_max=1)
True_com2, Polar_com2, _ = test_pi_centrer(masses2,coords2, box_min = 0, box_max=1)

height = 0.6
particles = 3

gs = gridspec.GridSpec(4, 1, width_ratios=[1], height_ratios=[1,2,1,2], wspace=0, hspace=0)
fig = plt.figure(figsize=(4.03, 3.3))
ax0 = fig.add_subplot(gs[0,0])
ax1 = fig.add_subplot(gs[1,0])
ax2 = fig.add_subplot(gs[2,0])
ax3 = fig.add_subplot(gs[3,0])


# Tick removal, horizontal line plotting
y1 = 0 
xmin = 0 
xmax = 0.99

for z in [ax0,ax2]:
    z.set_xlim(0,1)
    z.set_ylim(-1,1)
    z.hlines(y1, xmin, xmax, color = '#0173B2', alpha = 0.4, linewidth = blue_set)
    z.vlines(xmin, y1 - height / 2., y1 + height / 2.,color = '#0173B2', alpha = 0.4,linewidth = blue_set)
    z.vlines(xmax, y1 - height / 2., y1 + height / 2.,color = '#0173B2', alpha = 0.4,linewidth = blue_set)
    z.set_yticklabels([])
    z.set_yticks([])
    z.set_xticklabels([])
    z.set_xticks([])
    z.axis('off')

for l in [ax1,ax3]:
    l.set_xlim(0,1)
    l.set_ylim(0,1.05)
    l.spines['top'].set_visible(False)
    l.spines['right'].set_visible(False)

sine_wave2 = (ifft2 - np.min(ifft2))/np.ptp(ifft2)

ax3.plot(x2, sine_wave2, 'k--')

ax2.scatter(coords2,np.zeros(particles), s = s_set, color = ['black','#0173B2','black'], zorder=8)
ax2.plot(coords2,np.zeros(particles), color = 'grey', alpha = 0.8, linewidth = l_set)

ax2.axvline(coords2[1], ymax = 0.5, color = 'orange')
ax3.axvline(coords2[1], ymin = sine_wave2.max()/1.05, color = 'orange')

peak_x2 = x2[y2 == 1]
ax3.hist(peak_x2, bins=np.linspace(0, 1, 150), alpha=0.7, )

ax2.scatter(Polar_com2.round(3),0, s=s_set, color = 'orange', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=0.8,linewidths=cross_size)
ax2.scatter(True_com2,0, s=s_set, color = '#029E73', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=0.8,linewidths=cross_size)
ax3.set_xlabel('Particle coordinates')
ax3.set_ylabel('Mass density')

#-----------------------

ax1.spines['bottom'].set_visible(False)
ax1.set_xticks([])

# ---------------------------
sine_wave = (ifft - np.min(ifft))/np.ptp(ifft)

ax1.plot(x, sine_wave, 'k--')

ax0.scatter(coords,np.zeros(particles), s = s_set, color = ['black','#0173B2','black'], zorder=8)
ax0.plot(coords,np.zeros(particles), color = 'grey', alpha = 0.8, linewidth = l_set)

ax0.axvline(coords[1], ymax = 0.5, color = 'orange')
ax1.axvline(coords[1], ymin = sine_wave.max()/1.05, color = 'orange')

peak_x = x[y == 1]
ax1.hist(peak_x, bins=np.linspace(0, 1, 150), alpha=0.7)

ax0.scatter(Polar_com,0, s=s_set, color = 'orange', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=0.8,linewidths=cross_size)
ax0.scatter(True_com,0, s=s_set, color = '#029E73', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=0.8,linewidths=cross_size)
ax1.set_ylabel('Mass density')


plt.savefig(paths.figures / 'Fig4.pdf', bbox_inches='tight', pad_inches=0)
plt.close()
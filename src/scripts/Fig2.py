import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.transformations import wrap
from MDAnalysis.transformations.boxdimensions import set_dimensions
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.transformations.boxdimensions import set_dimensions
import matplotlib.gridspec as gridspec
import matplotlib as matplotlib
import matplotlib.patches as patches
from Functions import pi_centrer, xi_zeta_calc, xi_zeta_bar_calc, test_pi_centrer
from Variables import *
import paths

# Main

matplotlib.rcParams.update({'font.size': font_size})

MDA_store = np.zeros(100)
pi_store = np.zeros(100)
corrected_store = np.zeros(100)

# range x moves through
xrange = np.linspace(5.0001, 14.999, 100)
xrange[50] = 5.0001 + (
    14.999 - 5.001) / 2  # force part of the xrange to be the zero errror point
masses = np.array((1, 1, 1))

for i, x in enumerate(xrange):

    # Positions of particles in system
    positions = np.array((5.0001, x, 14.999))

    loop_md, loop_pi, loop_cor = test_pi_centrer(masses, positions)
    MDA_store[i] = loop_md
    pi_store[i] = loop_pi
    corrected_store[i] = loop_cor

gs = gridspec.GridSpec(2,
                       1,
                       width_ratios=[1],
                       height_ratios=[1, 5],
                       wspace=0,
                       hspace=0)
fig = plt.figure(figsize=(4.03, 3.3))
ax1 = fig.add_subplot(gs[0, 0])

# plotting options
xmin = 4
xmax = 16
s_set = 50

# Tick removal, horizontal line plotting
ax1.set_xlim(4, 16.1)
ax1.set_ylim(-1, 1)
ax1.hlines(y1, xmin, xmax, color='#0173B2', alpha=0.4, linewidth=blue_set)
ax1.set_yticklabels([])
ax1.set_yticks([])
ax1.set_xticklabels([])
ax1.set_xticks([])
ax1.axis('off')

# Selecting for specific molecule
particles = 3
coords = np.array([xrange[0], xrange[-1], xrange[20]])
masses = np.array([1, 1, 1])

ax1.scatter(coords,
            np.zeros(particles),
            s=s_set,
            color=['black', 'black', '#0173B2'],
            zorder=8)
ax1.plot(coords, np.zeros(particles), color='grey', alpha=0.8, linewidth=l_set)
ax1.scatter(MDA_store[20],
            0,
            s=s_set,
            color='#029E73',
            marker='x',
            label='Naive and Molecule Centre of Mass',
            zorder=10,
            alpha=0.8,
            linewidths=cross_size)
ax1.scatter(pi_store[20],
            0,
            s=s_set,
            color='#D55E00',
            marker='x',
            label='Naive and Molecule Centre of Mass',
            zorder=10,
            alpha=0.8,
            linewidths=cross_size)

ax0 = fig.add_subplot(gs[1, 0])

ax0.plot(xrange / 20,
         np.abs(pi_store - MDA_store) / 10,
         label='moving particle',
         color='#0173B2',
         solid_capstyle='round')
ax0.set_xlabel('Particle coordinates')
ax0.set_ylabel('Normalised centre of mass error')
ax0.set_ylim(0, 0.35)
ax0.set_xlim(0.2, 0.8)
ax0.set_yticks([0, 0.1, 0.2, 0.3])
ax0.set_xticks([0.25, 0.5, 0.75])

# Demonstration line
ax1.axvline(xrange[20], ymax=0.5, color='orange', linestyle='dashed')
ax0.axvline(xrange[20] / 20 - 0.0013,
            ymin=(np.abs(pi_store - MDA_store) / 10)[20] / 0.35,
            color='orange',
            linestyle='dashed')
ax0.axhline((np.abs(pi_store - MDA_store) / 10)[20],
            xmax=0.25,
            color='orange',
            linestyle='dashed',
            label='example error')
ax0.scatter(-20,
            20,
            s=s_set,
            color='orange',
            marker='x',
            label='B&B COM',
            zorder=10,
            alpha=0.8,
            linewidths=cross_size)
ax0.scatter(-20,
            -20,
            s=s_set,
            color='#029E73',
            marker='x',
            label='True COM',
            zorder=10,
            alpha=0.8,
            linewidths=cross_size)

ax0.spines['top'].set_visible(False)
plt.savefig(paths.figures / 'Fig2.pdf', bbox_inches='tight', pad_inches=0)
plt.close()

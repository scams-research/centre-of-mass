
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

# Function definitions 
def pi_centrer(coords,weights,box_min, box_width):
    """
    Calculates the centre of mass of a given set of coordinates and masses according to the bai and breen method.
    
    Parameters:
        coords (np.array): 1d array of particle coordinates
        weights (np.array): 1d array of particle masses
        box_min (float): Minimum x-coordinate of the box.
        box_width (float): Width of the box

    """
    frac_coords = (coords - box_min) / box_width
    theta = frac_coords * (2 * np.pi) 
    xi = np.cos(theta)
    zeta = np.sin(theta)
    xi_bar = np.average(xi,weights=weights)
    zeta_bar = np.average(zeta,weights=weights)
    theta_bar = np.arctan2(-zeta_bar, -xi_bar) + np.pi
    new_s_coords = (theta_bar) / (2 * np.pi)
    new_s_coords = new_s_coords*box_width + box_min

    return new_s_coords


def xi_zeta_calc(masses,coords, box_min = 0,box_width = 20):
    """
    Calculates the xi and zeta used within the bai and breen method.
    
    Parameters:
        coords (np.array): 1d array of particle coordinates
        weights (np.array): 1d array of particle masses
        box_min (float): Minimum x-coordinate of the box.
        box_width (float): Width of the box

    """

    frac_coords = (coords - box_min) / box_width 
    theta = frac_coords * (2 * np.pi) + np.pi
    xi = np.cos(theta)
    zeta = np.sin(theta)
    return xi, zeta, theta

def xi_zeta_bar_calc(masses,coords, box_min = 0, box_width = 20):
    """
    Calculates the xi bar, zeta bar and theta bar used within the bai and breen method.
    
    Parameters:
        coords (np.array): 1d array of particle coordinates
        weights (np.array): 1d array of particle masses
        box_min (float): Minimum x-coordinate of the box.
        box_width (float): Width of the box

    """
    frac_coords = (coords - box_min) / box_width
    theta = frac_coords * (2 * np.pi) + np.pi
    xi = np.cos(theta)
    zeta = np.sin(theta)    
    xi_bar = np.average(xi,weights=masses)
    zeta_bar = np.average(zeta,weights=masses)
    theta_bar = np.arctan2(-zeta_bar, -xi_bar) 
    return xi_bar, zeta_bar, theta_bar


def test_pi_centrer(masses,x_coords, box_min = 0, box_max = 20):
    """
    Calculates centre of mass of the given masses and coordinates by the Bai and Breen method, the method given in the paper and then using MDAnalysises unwrapping method that requires 
    bond information. 
    
    Parameters:
        coords (np.array): 1d array of particle coordinates
        weights (np.array): 1d array of particle masses
        box_min (float): Minimum x-coordinate of the box.
        box_max (float): Maximum x-coordinate of the box.

    """

    box_width = box_max - box_min
    particles = len(masses)

    yz = np.zeros((2, particles))
    x_coords = x_coords % 20
    coords = np.concatenate((x_coords[np.newaxis, :], yz), axis=0).T


    u = mda.Universe.empty(particles ,trajectory = True)
    u.add_TopologyAttr('masses',masses)
    u.transfer_to_memory()
    
    reader = MemoryReader(coords)
    u.trajectory = reader
    dim = np.array([box_max, box_max, box_max, 90, 90, 90])
    transform1 = mda.transformations.boxdimensions.set_dimensions(dim)
    transform2 = wrap(u.atoms)
    workflow  = [transform1,transform2]
    u.trajectory.add_transformations(*workflow)

    u.add_bonds([tuple(range(i, i+2)) for i in range(0, particles-1)])

    u.atoms.unwrap()
    MDA_com = u.atoms.center_of_mass()[0]
    com_pi = pi_centrer(coords[:,0], masses,box_min,box_width)
    corrected_com = (np.average(((coords[:,0] - (com_pi + 0.5*box_width) ) % box_max), weights = masses) + (com_pi + 0.5*box_width)) % box_max
    wrapped = ((coords[:,0] - (com_pi + 0.5*box_width) ) % box_max)

    return MDA_com, com_pi, corrected_com 


# Main

MDA_store = np.zeros(100)
pi_store = np.zeros(100)
corrected_store = np.zeros(100)

# range x moves through 
xrange = np.linspace(5.0001,14.999,100)
xrange[50] = 5.0001 + (14.999-5.001) / 2 # force part of the xrange to be the zero errror point
masses = np.array((1,1,1))

for i, x in enumerate(xrange):

    # Positions of particles in system
    positions = np.array((5.0001,x,14.999))

    loop_md,loop_pi,loop_cor = test_pi_centrer(masses, positions)
    MDA_store[i] = loop_md
    pi_store[i] = loop_pi
    corrected_store[i] = loop_cor


gs = gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[1,5], wspace=0, hspace=0)
fig = plt.figure(figsize=(4.03, 3.3))
ax1 = fig.add_subplot(gs[0,0])

# plotting options
xmin = 4
xmax = 16
y1 = 0
height = 1
blue_set = 0.5
l_set = 2
s_set = 50
marker_set = 15
font_size = 7
cross_size = 1.5

# Tick removal, horizontal line plotting
for x in [ax1]:
    x.set_xlim(4,16.1)
    x.set_ylim(-1,1)
    x.hlines(y1, xmin, xmax, color = '#0173B2', alpha = 0.4, linewidth = blue_set)
    x.set_yticklabels([])
    x.set_yticks([])
    x.set_xticklabels([])
    x.set_xticks([])
    x.axis('off')


# Selecting for specific molecule
particles = 3
coords = np.array([xrange[0],xrange[-1],xrange[20]])
masses =np.array([1,1,1])

ax1.scatter(coords,np.zeros(particles), s = s_set, color = ['black','black','#0173B2'], zorder=8)
ax1.plot(coords,np.zeros(particles), color = 'grey', alpha = 0.8, linewidth = l_set)
ax1.scatter(MDA_store[20],0, s=s_set, color = '#029E73', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=0.8,linewidths=cross_size)
ax1.scatter(pi_store[20],0, s=s_set, color = '#D55E00', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=0.8,linewidths=cross_size)


ax0 = fig.add_subplot(gs[1,0])



ax0.plot(xrange/20,np.abs(pi_store - MDA_store)/10, label = 'moving particle', color = '#0173B2',solid_capstyle='round')
ax0.set_xlabel('Particle coordinates')
ax0.set_ylabel('Centre of mass error / Span of molecule')
ax0.set_ylim(0,0.35)
ax0.set_xlim(0.2,0.8)



# Demonstration line
ax1.axvline(xrange[20], ymax = 0.5, color = 'orange', linestyle = 'dashed')
ax0.axvline(xrange[20]/20-0.0013, ymin = (np.abs(pi_store - MDA_store)/10)[20]/0.35, color = 'orange', linestyle = 'dashed')
ax0.axhline((np.abs(pi_store - MDA_store)/10)[20], xmax = 0.25, color = 'orange', linestyle = 'dashed', label = 'example error')
ax0.scatter(-20,20, s=s_set, color = 'orange', marker = 'x', label = 'B&B COM',zorder=10, alpha=0.8,linewidths=cross_size)
ax0.scatter(-20,-20, s=s_set, color = '#029E73', marker = 'x', label = 'True COM',zorder=10, alpha=0.8,linewidths=cross_size)



ax0.spines['top'].set_visible(False)
plt.show()
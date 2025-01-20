import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.transformations import wrap 
from MDAnalysis.transformations.boxdimensions import set_dimensions
import matplotlib.gridspec as gridspec
import matplotlib 
from matplotlib.patches import Arc


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
    return xi, zeta

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


def draw_dashed_box(ax, x_min, x_max, y_min, y_max, color='black', linestyle='--', linewidth = 0.5):
    """
    Draws a dashed rectangular box on the given axis.
    
    Parameters:
        ax (matplotlib.axes.Axes): The axis to draw the box on.
        x_min (float): Minimum x-coordinate of the box.
        x_max (float): Maximum x-coordinate of the box.
        y_min (float): Minimum y-coordinate of the box.
        y_max (float): Maximum y-coordinate of the box.
        color (str): Color of the dashed lines.
        linestyle (str): Style of the dashed lines.
    """
    # Draw the edges of the box
    ax.plot([x_min, x_max], [y_min, y_min], linestyle=linestyle, color=color, linewidth = linewidth)  # Bottom edge
    ax.plot([x_min, x_max], [y_max, y_max], linestyle=linestyle, color=color, linewidth = linewidth)  # Top edge
    ax.plot([x_min, x_min], [y_min, y_max], linestyle=linestyle, color=color, linewidth = linewidth)  # Left edge
    ax.plot([x_max, x_max], [y_min, y_max], linestyle=linestyle, color=color, linewidth = linewidth)  # Right edge


# Main

fig = plt.figure(figsize=(8.06, 8.06/4.6666667))

# box coordinates and line
xmin = 0
xmax = 20
y1 = 0
height = 1


gs = gridspec.GridSpec(6, 5, width_ratios=[1.9,0.18,1.05,0.18, 0.9], height_ratios=[1,1,1,1,1,1], wspace=0, hspace=0.2)

# Axes for fig 1b 
ax0 = fig.add_subplot(gs[0:3, 2])  # Top-left
ax1 = fig.add_subplot(gs[3:6, 2])  # Bottom-left

# Axes for fig 1c
ax2 = fig.add_subplot(gs[:, 4]) 

# Axes for fig 1a 
ax3 = fig.add_subplot(gs[0:2,0]) 
ax4 = fig.add_subplot(gs[2:4,0])
ax5 = fig.add_subplot(gs[4:6,0])


 # Alpha set point
blue_set = 0.5


# Remove tick labels and plot horizontal and plot blue lines
for x in [ax0,ax1]:
    x.set_xlim(-1,21)
    x.set_ylim(-1,1)
    x.hlines(y1, xmin, xmax, color = '#0173B2', alpha = 0.4, linewidth = blue_set)
    x.vlines(xmin, y1 - height / 2., y1 + height / 2.,color = '#0173B2', alpha = 0.4, linewidth = blue_set)
    x.vlines(xmax, y1 - height / 2., y1 + height / 2.,color = '#0173B2', alpha = 0.4, linewidth = blue_set)
    x.set_yticklabels([])
    x.set_yticks([])
    x.set_xticklabels([])
    x.set_xticks([])
    x.axis('off')

# Standardized plotting point settings
l_set = 2
s_set = 35
marker_set = 15
font_size = 7
cross_size = 1.5


# ------------------  Fig 1a plotting ---------------------
particles = 3
coords = np.array([3,6,8])
masses =np.array([12,12,16])

# draw lines and set twice the box
xmax_2 = xmax*2
y1 = 0
height = 1

# Remove tick labels and plot horizontal and plot blue lines
for ax in [ax3,ax4,ax5]:
    ax.set_xlim(-1,41)
    ax.set_ylim(-1,1)
    ax.hlines(y1, xmin, xmax_2, alpha = 0.4,linewidth = blue_set)
    ax.vlines(xmin, y1 - height / 2., y1 + height / 2.,color = '#0173B2', alpha = 0.4,linewidth = blue_set)
    ax.vlines(xmin+20, y1 - height / 3., y1 + height / 3.,color = '#0173B2', alpha = 0.4,linewidth = blue_set)
    ax.vlines(xmax_2, y1 - height / 2., y1 + height / 2.,color = '#0173B2', alpha = 0.4,linewidth = blue_set)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.axis('off')


for i,ax in enumerate([ax3,ax4,ax5]):
    ax.scatter(coords,np.zeros(particles), s = s_set, color = ['black','black','red'], zorder=8, lw = 0)
    ax.scatter(coords+20,np.zeros(particles), s = s_set, color = ['black','black','red'], zorder=8, alpha = 0.65, lw=0)
    ax.plot(coords+20,np.zeros(particles), color = 'grey', alpha = 0.35, linewidth = l_set,solid_capstyle='round')
    ax.plot(coords,np.zeros(particles), color = 'grey', alpha = 0.8, linewidth = l_set,solid_capstyle='round')

    for i in range(particles):
        ax.text(coords[i], 0.32, f'{i+1}', color='black', ha='center', fontsize=font_size)
        ax.text(coords[i] + 20, 0.32, f'{i+1}', color='black', ha='center', fontsize=font_size)

draw_dashed_box(ax3, 0, 20, -0.25, 0.25)
draw_dashed_box(ax4, 4.5, 24.5, -0.25, 0.25)
draw_dashed_box(ax5,7, 27, -0.25, 0.25)

com0 = np.average(coords, weights=masses)
com1 = np.average([coords[1],coords[2],coords[0]+20],weights = [masses[1],masses[2],masses[0]])
com2 = np.average([coords[2],coords[0]+20,coords[1]+20],weights = [masses[2],masses[0],masses[1]])

for ax,com in zip([ax3,ax4,ax5],[com0,com1,com2]):
    ax.scatter(com, 0, s=s_set, color='#029E73', marker='x', label='Naive Periodic Centre of Mass', zorder=10,linewidths=cross_size)
    ax.scatter(com+20, 0, s=s_set, color='#029E73', marker='x', zorder=10, alpha=0.7,linewidths=cross_size)


# ------------------  Fig 1b & 1c plotting ---------------------
particles = 3
coords = np.array([7,10,12])
masses =np.array([12,12,16])

ax1.set_xlabel('Box coordinate')

ax0.scatter(coords,np.zeros(particles), s = s_set, color = ['black','black','red'], zorder=8, lw = 0)
ax0.plot(coords,np.zeros(particles), color = 'grey', alpha = 0.8, linewidth = l_set)
ax0.scatter(np.average(coords, weights=masses),0, s=s_set, color = '#029E73', marker = 'x', label = 'Naive and Molecule Centre of Mass',zorder=10, alpha=1,linewidths=cross_size)


# Wrapping coordinates for particles across periodic boundary and calculation of COM 
coords = (coords-11) % 20
mda_com, com_pi, corrected_com = test_pi_centrer(masses,coords)

# Plotting of points and grey line across periodic boundary
ax1.scatter(coords[0:2],np.zeros(particles)[0:2],s = s_set, color = ['black','black'], zorder=8, lw = 0)
ax1.plot(coords[0:2],np.zeros(particles)[0:2], markersize=marker_set, color = 'grey', zorder=2,linewidth = l_set,alpha = 0.8)
ax1.plot([coords[1],19.8],np.zeros(particles)[0:2], markersize=marker_set, color = "grey", zorder=2,linewidth = l_set,alpha = 0.8)
ax1.plot([coords[2],0],np.zeros(particles)[0:2], markersize=marker_set, color = "grey", zorder=2,linewidth = l_set,alpha = 0.8, mew=0)
ax1.scatter(coords[2],0, s=s_set, color = "red", zorder=8, lw = 0)
ax1.scatter(np.average(coords, weights=masses),0, s=s_set, color = 'red', marker = 'x', label = 'Naive Centre of Mass', zorder=10,alpha=0.8,linewidths=cross_size)
ax1.scatter(mda_com,0, s=s_set, color = '#029E73', marker = 'x', label = 'Molecule Centre of Mass', zorder=10,alpha=0.9,linewidths=cross_size)


# Circular projection axis set up
ax2.set_xlim(-1.1,1.06)
ax2.set_ylim(-1.1,1.06)

ax2.set_yticklabels([])
ax2.set_yticks([])
ax2.set_xticklabels([])
ax2.set_xticks([])
ax2.axis('off')

# Custom marker orthogonal for circle
t = matplotlib.markers.MarkerStyle(marker='x')
t._transform = t.get_transform().rotate_deg(15)


# Circle plotting
circle = plt.Circle((0, 0), 1, fill=False,color = '#0173B2', alpha = 0.4, linewidth = blue_set)  # `fill=False` for an unfilled circle
ax2.add_artist(circle)
ax2.axhline(00,-0.1,0.1, color = '#0173B2', alpha = 0.6, linewidth = blue_set, zorder = 1)

# Conversion of particle coordinates to a circular projection
xi, zeta = xi_zeta_calc(masses,coords)

# Negative zeta just flips direction of projection
ax2.scatter(xi,-zeta, color=['black','black','red'], s = s_set, zorder=10, lw = 0)
ax2.scatter(0,0, color='#0173B2', s = 6, zorder=10, alpha = 0.6, lw = 0)

# B&B method centre of mass plotting
xi_bar, zeta_bar, theta_bar = xi_zeta_bar_calc(masses,coords)
xi_calc,zeta_calc = xi_zeta_calc(masses,com_pi)
ax2.scatter(xi_bar,-zeta_bar, color='orange', marker = t, s = s_set, zorder=10, label = 'Circular Projection Centre of Mass',linewidths=cross_size)
ax2.plot([0,xi_calc],[0,-zeta_calc], color = 'orange', alpha = 0.8, linewidth = 1.2, solid_capstyle='round')


arc = Arc((0, 0), 2, 2, theta1=-200, theta2=250, color='grey', alpha=0.7, linewidth =l_set, zorder=1)
ax2.add_patch(arc)
ax2.axis('off')


fig.text(.30, .05, '(a)', ha='center', size = font_size)
fig.text(.605, .05, '(b)', ha='center', size = font_size)
fig.text(.818, .05, '(c)', ha='center', size = font_size)

plt.show()
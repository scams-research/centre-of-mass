# Function definitions 
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.transformations import wrap 
from MDAnalysis.transformations.boxdimensions import set_dimensions
from geomstats.learning.frechet_mean import CircleMean
from geomstats.geometry.hypersphere import Hypersphere

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


def frechet_statistic(theta,theta0):
    return(np.minimum(np.abs(theta-theta0),1-np.abs(theta-theta0)))**2

def intrinsic_mean(coords,masses):
    """
    Assumes fractional coords i.e. 0 to 1

    Calculates analytical intrinsic mean via Frechet function minimisation

    """
    theta = coords #* 2 * np.pi
    k = np.arange(len(coords))
    ngon_means  = (np.average(theta,weights=masses) +  k/len(coords)) % 1 
    frechet_stat = np.array([np.sum(masses*frechet_statistic(means, theta)) for means in ngon_means])
    intrinsic_center_of_mass = ngon_means[np.where(frechet_stat==frechet_stat.min())]# / (2*np.pi)

    return intrinsic_center_of_mass[0]

def geom_stats_intrinsic(particles,circle_mean):
    """
    Calculates intrinsic mean via the geomstats implementation
    """
    theta  = (particles % 1) * (2 * np.pi)
    xi = np.cos(theta)
    zeta = np.sin(theta)
    data = np.array([xi, zeta]).T
    circle_mean.fit(data)
    x = circle_mean.estimate_
    return x
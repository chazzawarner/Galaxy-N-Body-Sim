import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
from latex2sympy2 import latex2sympy as l2s
from star_types import get_random_star
from astropy import units as u

# Method taken from https://bjerkeli.se/onewebmedia/thesis.pdf

# Define constants
G = 6.67430e-11 * u.m**3 / (u.kg * u.s**2)  # Gravitational constant

# Returns rotational velocity function at radius r
def v2_circ_from_potential(potential):
    # v2_circ = r * dPhi/dr
    G, M, a, r = sp.symbols("G M a r")
    v2_circ = sp.diff(potential, r) * r
    v2_circ = sp.lambdify((G, M, a, r), v2_circ, "numpy")
    
    return v2_circ # Returns a function that takes in G, M, a, and r and returns the circular velocity at that radius
    
# Define the equations of the Hernquist distribution
hernquist = {
    "potential": l2s(r"\Phi_r=-\frac{GM}{r+a}"), # Gravitational potential (Phi(r))
    "density": l2s(r"\rho_r=\frac{M}{2\pi}\frac{a}{r}\frac{1}{(r+a)^3}"), # Density distribution (rho(r))
    "mass_frac": l2s(r"M_f=\frac{r^2}{(r+a)^2}"), # Mass fraction (M_f = M(r)/M
}

miyamoto = {
    "potential": l2s(r"\Phi_{R,Z}=-\frac{GM}{\sqrt{R^2+(a+\sqrt{z^2+b^2})^2}}"), # Gravitational potential (Phi(R, z))
    "density": l2s(r"\rho_{R,Z}= \frac{b^2M}{4\pi} \frac{aR^2 + (a+3\sqrt{z^2+b^2}) (a+\sqrt{z^2+b^2})^2 }{ (R^2 + (a+\sqrt{z^2+b^2})^2)^{5/2} (z^2+b^2)^{3/2} }"), # Density distribution (rho(R, z))
}


def gen_points(num_bodies, density_func):
    # Use 


def main():
    # Test the elliptical galaxy generation function
    radius = 40e3 # parsecs
    num_bodies = 5
    total_mass = 1e9 # solar masses
    G_const = G.to_value(u.pc**3 / (u.M_sun * u.s**2))
    print(f"G: {G_const}")
    
    """# Plot density distribution
    radii = np.linspace(0.1, 10, 1000)
    density = hernquist["density"]
    M, a, r = sp.symbols("M a r")
    density = sp.lambdify((M, a, r), density, "numpy")
    density = density(100, 10, radii)
    
    #print(density)
    
    plt.plot(radii, density)
    plt.yscale("log")
    plt.show()"""
    

    R = np.linspace(0.1, 50, 1000)
    z = np.linspace(-20, 20, 1000)
    R_symbol, z_symbol, a, b, M = sp.symbols("R z a b M")
    density = miyamoto["density"]
    density = sp.lambdify((R_symbol, z_symbol, a, b, M), density, "numpy")

    # Create a 2D grid of R and z values
    R_grid, z_grid = np.meshgrid(R, z)

    # Evaluate the density function on the grid
    density = density(R_grid, z_grid, 1, 1, 1000)

    # Plot the Miyaomoto density distribution as a contour plot
    """plt.contour(R_grid, z_grid, density)
    plt.colorbar()
    plt.clabel()
    plt.xlabel("R")
    plt.ylabel("z")
    plt.title("Miyaomoto Density Distribution")
    plt.show()"""
    
    fig, ax = plt.subplots()
    levels = [0.0001, 0.001, 0.01, 0.1]
    CS = ax.contour(R_grid, z_grid, density, levels=levels)
    ax.clabel(CS, inline=1, fontsize=10)
    ax.set_title('Miyaomoto Density Distribution')
    plt.show()
    
    
    
    #positions, radii = elliptical_galaxy_gen(num_bodies, radius, total_mass, G=G_const, a=0.45)
    
    """# Plot a histogram of the radii
    plt.hist(radii, bins=100, density=True)
    plt.xscale("log")
    plt.xlabel("Radius")
    plt.ylabel("Density")
    plt.show()
    """
    
    """# Plot the positions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2])
    
    ## Plot a point at the origin
    ax.scatter(0, 0, 0, color="r")
    
    ## Set the aspect ratio of the plot to be equal
    ax.set_box_aspect([1,1,1])
    
    ## Limit the plot to the radius
    ax.set_xlim(-radius, radius)
    
    plt.show()"""
    
    # Remove points outside the radius*1.5
    positions = positions[np.linalg.norm(positions, axis=1) < radius*1.5]
    
    # Plot a 2x2 grid of views of the galaxy
    fig, axs = plt.subplots(2, 2)

    ## 3D plot
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], s=1)
    ax.set_title("3D view")

    ## Remove empty graph/axes for the 3D plot
    fig.delaxes(axs[0, 0])

    ## 2D Side view
    ax = axs[0, 1]
    ax.grid(True)
    ax.scatter(positions[:, 1], positions[:, 2], s=1)
    ax.set_title("2D Side view")
    ax.set_aspect('equal', 'box')

    ## 2D Top view
    ax = axs[1, 0]
    ax.grid(True)
    ax.scatter(positions[:, 0], positions[:, 1], s=1)
    ax.set_title("2D Top view")
    ax.set_aspect('equal', 'box')
    

    ## 2D Front view
    ax = axs[1, 1]
    ax.grid(True)
    ax.scatter(positions[:, 0], positions[:, 2], s=1)
    ax.set_title("2D Front view")
    ax.set_aspect('equal', 'box')

    plt.show()
    
main()
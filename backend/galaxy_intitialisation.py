import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
from latex2sympy2 import latex2sympy as l2s
from star_types import get_random_star
from astropy import units as u
from galpy.potential import HernquistPotential

# Method taken from https://bjerkeli.se/onewebmedia/thesis.pdf

# Define constants
G = 6.67430e-11 * u.m**3 / (u.kg * u.s**2)  # Gravitational constant


# Define the equations of the Hernquist distribution
hernquist = {
    "density_dist": l2s(r"\frac{M}{2\pi} \frac{a}{r} \frac{1}{(r+a)^3}"), # Density distribution (rho(r))
    "mass_frac": l2s(r"M_f=\frac{r^2}{(r+a)^2}"), # Mass fraction (M_f = M(r)/M)
    "grav_pot": l2s(r"-\frac{GM}{r+a}"), # Gravitational potential (Phi(r))
    "escape_vel": l2s(r"(\frac{2GM}{r+a})^{0.5}"), # Escape velocity i.e. max. velocity (v_e(r)
    "dist_func": {
        "function": l2s(r"\frac{M}{8\sqrt{2}\pi^3a^3{v_g}^3} \frac{1}{(1-q^2)^{5/2}} (3\arcsin{q} + q(1-q)^{1/2} (1-2q^2) (8q^4-8q^2-3) )"), # Distribution function (f(E))
        "q": l2s(r"\sqrt{\frac{-a * T}{G * M}}"), # q = sqrt(-aE/GM)
        "v_g": l2s(r"\sqrt{ \frac{G * M}{a} }"), # v_g = GM/a
    }
}

# Define elliptical galaxy generation function using the Hernquist distribution
#def elliptical_galaxy_gen(num_bodies, max_radius, total_mass, a=1, G=G.value):
def elliptical_galaxy_gen(num_bodies, a=1):
    # Generate positions
    ## Draw random mass fractions and solve for radius
    mass_fractions = np.random.rand(num_bodies)
    
    M_r, r, a_symbol = sp.symbols("M_r r a")
    mass_frac_eq = sp.Eq(hernquist["mass_frac"], M_r)
    radius_solutions = sp.solve(mass_frac_eq, r)
    radius_solutions_func = sp.lambdify((M_r, a_symbol), radius_solutions[0])
    #print(radius_solutions)
    
    radius = radius_solutions_func(mass_fractions, a)
    #print(f"Radius: {radius}")
    
    return radius
    
    """plt.hist(radius, bins=100, density=True)
    plt.xscale("log")
    plt.show()
    
    ## Randomly distribute on a sphere of corresponding radius
    ### Generate random positions
    
    positions = np.zeros((num_bodies, 3))
    positions[:, 0] = np.random.normal(0, 1, num_bodies)
    positions[:, 1] = np.random.normal(0, 1, num_bodies)
    positions[:, 2] = np.random.normal(0, 1, num_bodies)
    
    ### Normalise the positions
    norms = np.linalg.norm(positions, axis=1)
    positions = positions / norms.reshape(-1, 1)
    
    ### Scale the positions to the radius
    positions = positions * radius.reshape(-1, 1)
    
    
    # Generate masses
    ## Draw random masses
    masses = np.zeros(num_bodies)
    for i in range(num_bodies):
        masses[i] = get_random_star()

    ## Scale masses to the total mass
    masses = masses / np.sum(masses) * total_mass
    
    
    # Generate velocities for each position
    velocities = np.zeros((num_bodies, 3))
    
    # V^2_circ = GMr(r+a)^-2
    
    return positions, radius"""

def main():
    # Test the elliptical galaxy generation function
    radius = 40e3 # parsecs
    num_bodies = 1000
    total_mass = 1e9 # solar masses
    G_const = G.to_value(u.pc**3 / (u.M_sun * u.s**2))
    print(f"G: {G_const}")
    positions, radii = elliptical_galaxy_gen(num_bodies, radius, total_mass, G=G_const, a=0.45)
    
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
    
if __name__ == '__main__':
    main()
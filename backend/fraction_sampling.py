import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
from latex2sympy2 import latex2sympy as l2s
from galpy.potential import HernquistPotential
from scipy.optimize import fsolve
from scipy.optimize import minimize

# Define equations of the Hernquist distribution
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

# Sample Hernquist potential using mass fractions
def fraction_sampling(num_bodies, a=1, rejection_threshold=None):
    # Generate positions
    ## Draw random mass fractions and solve for radius
    mass_fractions = np.random.uniform(0, 1, num_bodies)
    
    radius_solution_func = lambda M_f, a: -a * ((np.sqrt(M_f) + M_f) / (np.sqrt(M_f) - 1))
    radius_solution_func = np.vectorize(radius_solution_func)
    
    radius = radius_solution_func(mass_fractions, a)
    #print(f"Radius: {radius}")
    
    if rejection_threshold is not None:
        # Reject samples that are above the rejection threshold
        radius = radius[radius < rejection_threshold]

    return radius

def main():
    num_bodies = 100000
    a = 0.45/10
    positions = fraction_sampling(num_bodies, a=a, rejection_threshold=10)
    
    
    radii = np.linspace(1e-5, np.max(positions), 10000)
    true_density = HernquistPotential(a=a, normalize=1).dens(radii, 0)
    area_under_density = np.trapz(true_density, radii)
    print(f"Area under density: {area_under_density}")
    normalised_density = true_density / area_under_density
    
    area_under_density = np.trapz(normalised_density, radii)
    print(f"Area under normalised density: {area_under_density}")
    
    
    fig, ax1 = plt.subplots()
    
    # Plot true density on first y-axis
    ax1.plot(radii, normalised_density, label="True density", color='blue')
    ax1.set_xlabel("Radius")
    ax1.set_ylabel("True Density", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_yscale("log")
    
    # Create second y-axis for histogram
    ax2 = ax1.twinx()
    
    # Plot sampled radii histogram on second y-axis
    ax2.hist(positions, bins=100, alpha=0.5, density=True, label="Sampled radii", edgecolor='black', color='orange')
    ax2.set_ylabel("Sampled Radii Density", color='orange')
    ax2.tick_params(axis='y', labelcolor='orange')
    ax2.set_yscale("log")
    
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':
    main()
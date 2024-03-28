import numpy as np
from matplotlib import pyplot as plt
from galpy.potential import HernquistPotential

"""
Using Hernquist potential, we can get the density function:
    rho(r) = M / (2 * pi) * a / r * 1 / (r + a)^3
Integrate over the volume and divide by the total mass to get mass fraction:
    M(r) / M = M_f = r^2 / (r + a)^2
Solve for r in terms of M_f to get the radius:
    r = -a * ((sqrt(M_f) + M_f) / (sqrt(M_f) - 1))
Can sample uniform mass fractions and solve for radius to get a sample of radii.
"""


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
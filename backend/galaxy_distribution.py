import numpy as np
import matplotlib.pyplot as plt
from star_types import get_random_star, get_average_star_mass

# Define Einasto Law for density profile of galaxy components
def einasto_law(a, rho_c, d_N, a_c, N):
    rho_a = rho_c * np.exp(-d_N * ((a / a_c)**N - 1))
    return rho_a

# Enclosed mass distribution using density profile
def enclosed_mass_distribution_density(radius_array, density_array):
        mass_distribution = np.zeros_like(radius_array)
        for i, radius in enumerate(radius_array):
            mass_distribution[i] = np.trapz(4 * np.pi * radius_array[:i+1]**2 * density_array[:i+1], radius_array[:i+1])
        return mass_distribution
    
"""# Enclosed mass distribution using bodies
def enclosed_mass_distribution_bodies(radius_array, bodies):
    mass_distribution = np.zeros_like(radius_array)
    for i, radius in enumerate(radius_array):
        mass_distribution[i] = np.sum([bodies["mass"][j] for j in range(len(bodies["mass"])) if bodies["distance"][j] <= radius])
    return mass_distribution"""

def enclosed_mass(radius, radius_array, mass_distribution):
    return np.interp(radius, radius_array, mass_distribution)

# Normalise mass distribution
def normalise_mass_distribution(mass_distribution, target_total_mass):
    normalised_distribution = mass_distribution * (target_total_mass / mass_distribution[-1])
    return normalised_distribution

def rotational_velocity(radius, mass_distribution, radius_array, g_const):
    enclosed_mass_at_radius = enclosed_mass(radius, radius_array, mass_distribution)
    velocity = np.sqrt(g_const * enclosed_mass_at_radius / radius)
    return velocity

# Inverse cdf of mass distribution
def inverse_mass_distribution(mass_distribution, max_distance, num_points=1000):
    radii = np.linspace(0, max_distance, num_points)
    
    cdf = np.cumsum(mass_distribution)
    cdf = cdf / cdf[-1]
    
    inverse_cdf_func = lambda p: np.interp(p, cdf, radii)
    
    return inverse_cdf_func

# Main function to generate galaxy (currently using Andromeda galaxy parameters)
def galaxy_generation(num_bodies=100, galaxy_diameter=46.56e3, g_const=4.3009e-3):
    num_bodies = num_bodies
    max_distance = galaxy_diameter/2 # Maximum distance from the centre of the galaxy in parsecs
    centre_vector = np.array([galaxy_diameter/2, galaxy_diameter/2])
    
    # Define mass distribution galaxy components
    bulge = {
        "a_c": 2.025e3,
        "q": 0.73,
        "N": 4,
        "d_N": 11.67,
        "rho_c": 2.2e-1,
        "M_comp": 4.9e10
    }
    
    disc = {
        "a_c": 11.35e3,
        "q": 0.1,
        "N": 1,
        "d_N": 2.67,
        "rho_c": 1.72e-2,
        "M_comp": 4.8e10
    }
    
    # Calculate the mass distribution of the galaxy
    def a(r, q, z=0):
        return np.sqrt(r**2 + (z/q)**2)
    
    radii = np.linspace(0, max_distance, 1000)
    
    bulge_density_profile = np.array([[r, einasto_law(a(r, bulge["q"]), bulge["rho_c"], bulge["d_N"], bulge["a_c"], bulge["N"])] for r in radii])
    disc_density_profile = np.array([[r, einasto_law(a(r, disc["q"]), disc["rho_c"], disc["d_N"], disc["a_c"], disc["N"])] for r in radii])
    
    bulge_mass_distribution = normalise_mass_distribution(enclosed_mass_distribution_density(radii, bulge_density_profile[:, 1]), bulge["M_comp"])
    disc_mass_distribution = normalise_mass_distribution(enclosed_mass_distribution_density(radii, disc_density_profile[:, 1]), disc["M_comp"])
    
    total_mass_distribution = bulge_mass_distribution + disc_mass_distribution
    
    
    print(f"Total mass of galaxy: {total_mass_distribution[-1]}")
    
    # Generate bodies
    inverse_cdf_func = inverse_mass_distribution(total_mass_distribution, max_distance)
    
    random_radii = np.random.rand(num_bodies)
    
    sampled_radii = inverse_cdf_func(random_radii)
    
    #return sampled_radii, total_mass_distribution, radii, bulge_mass_distribution, disc_mass_distribution

    # Generate random directions for bodies
    direction = np.random.rand(num_bodies, 2) * 2 - 1
    direction /= np.linalg.norm(direction, axis=1).reshape(-1, 1)
    
    # Generate positions for bodies
    positions = direction * sampled_radii.reshape(-1, 1) + centre_vector
    
    # Find the velocity of the bodies
    # Find the unit vector perpendicular to the direction of the body
    velocities_direction = np.array([-direction[:, 1], direction[:, 0]]).T
    velocities_direction /= np.linalg.norm(velocities_direction, axis=1).reshape(-1, 1)
    
    # Calculate the rotational velocity of the bodies
    velocities = velocities_direction * rotational_velocity(sampled_radii, total_mass_distribution, radii, g_const).reshape(-1, 1)
    
    # Generate the masses of the bodies, scaling so their total mass is that of the galaxy distribution
    masses = np.array([get_random_star() for _ in range(num_bodies)])
    mass_scale = total_mass_distribution[-1] / np.sum(masses)
    masses *= mass_scale
    
    
    # Return bodies
    bodies = {
        "mass": masses,
        "position": positions,
        "velocity": velocities
    }
    return bodies

def main():
    """sampled_radii, total_mass_distribution, radii, bulge_mass_distribution, disc_mass_distribution = galaxy_generation(num_bodies=100000)

    fig, ax1 = plt.subplots()

    # Plot distribution of stars on the first y-axis
    color = 'tab:blue'
    ax1.set_xlabel('Radii')
    ax1.set_ylabel('Count', color=color)
    ax1.hist(sampled_radii, bins=100, color=color, density=True, alpha=0.5, label="Sampled radii")
    ax1.tick_params(axis='y', labelcolor=color)

    # Create a second y-axis that shares the same x-axis
    ax2 = ax1.twinx()

    # Plot mass distribution on the second y-axis
    color = 'tab:red'
    ax2.set_ylabel('Mass Distribution', color=color)
    ax2.plot(radii, total_mass_distribution, color='tab:orange', label="Total mass distribution")
    ax2.plot(radii, bulge_mass_distribution, color='tab:green', label="Bulge mass distribution")
    ax2.plot(radii, disc_mass_distribution, color='tab:red', label="Disc mass distribution")
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.legend()

    fig.tight_layout()
    plt.show()"""
    
    andromeda_normal = galaxy_generation(num_bodies=10000, galaxy_diameter=46.56e3)
    a_normal_centre = np.array([46.56e3/2, 46.56e3/2])
    
    andromeda_expanded = galaxy_generation(num_bodies=100000, galaxy_diameter=46.56e3 * 2)
    a_expanded_centre = np.array([46.56e3, 46.56e3])
    
    andromeda_massive = galaxy_generation(num_bodies=100000, galaxy_diameter=46.56e3 * 10)
    a_massive_centre = np.array([46.56e3 * 5, 46.56e3 * 5])
    
    # Plot distribution of stars
    radii_normal = np.linalg.norm(andromeda_normal["position"] - a_normal_centre, axis=1)
    radii_normal /= 46.56e3
    plt.hist(radii_normal, bins=100, label="Usual diameter", alpha=0.5, density=True)
    
    andromeda_expanded_radii = np.linalg.norm(andromeda_expanded["position"] - a_expanded_centre, axis=1)
    andromeda_expanded_radii /= 46.56e3 * 2
    plt.hist(andromeda_expanded_radii, bins=100, label="Expanded diameter (x2)", alpha=0.5, density=True)
    
    andromeda_massive_radii = np.linalg.norm(andromeda_massive["position"] - a_massive_centre, axis=1)
    andromeda_massive_radii /= 46.56e3 * 10
    plt.hist(andromeda_massive_radii, bins=100, label="Massive diameter (x10)", alpha=0.5, density=True)
    
    plt.xlabel("Distance from centre (parsecs)")
    plt.ylabel("Count")
    
    plt.show()
    
    """# Plot galaxy
    plt.scatter(sampled_radii, np.zeros_like(sampled_radii), s=1)
    plt.show()"""
    
    galaxy_diameter = 46.56e3
    bodies = galaxy_generation(num_bodies=10000, galaxy_diameter=galaxy_diameter)
    
    """# Calculate the total mass of the bodies
    total_mass = np.sum(bodies["mass"])
    print(f"Total mass: {total_mass}")
    
    # Plot histogram of velocities of the bodies over distance from the centre
    centre_vec = np.array([galaxy_diameter/2, galaxy_diameter/2])
    distances = np.linalg.norm(bodies["position"] - centre_vec, axis=1)
    velocities = np.linalg.norm(bodies["velocity"], axis=1)
    
    plt.scatter(distances, velocities, s=1)
    
    plt.xlabel("Distance from centre (parsecs)")
    plt.ylabel("Velocity (parsecs/year)")
    
    plt.show()"""
    
    """# Plot the positions of the bodies
    plt.scatter(bodies["position"][:, 0], bodies["position"][:, 1], s=1)
    plt.show()"""
    
    
#main()
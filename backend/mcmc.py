import numpy as np
import matplotlib.pyplot as plt
from galpy.potential import MiyamotoNagaiPotential, HernquistPotential
from tqdm import tqdm

# Metropolis-Hastings MCMC algorithm for multivariate distributions i.e. Miyaoto-Nagai potential
def metropolis_hastings(target_density, dimensions, n_samples, burn_in=1000):
    # Initialise the chain
    n_samples += burn_in
    x0 = np.zeros(dimensions)
    xt = x0
    samples = np.zeros((n_samples, dimensions))
    
    # Run the chain
    for i in range(n_samples):
        # Generate a candidate sample
        xt_candidate = xt + np.random.normal(size=dimensions)
        
        # Accept or reject the candidate sample
        acceptance_probability = min(1, target_density(*xt_candidate) / target_density(*xt))
        
        # Accept the candidate sample with probability acceptance_probability
        if np.random.uniform() < acceptance_probability:
            xt = xt_candidate
            
        # Add the sample to the chain
        samples[i] = xt
            
    # Return the samples after burn-in
    return samples[burn_in:]


def main():
    """P = lambda x: 3 * np.exp(-x*x/2) + np.exp(-(x - 4)**2/2)
    P_norm = lambda x: P(x) / (3 * np.sqrt(2 * np.pi) + np.sqrt(2 * np.pi))
    samples_1k = metropolis_hastings(P, 1, 1000)
    samples_10k = metropolis_hastings(P, 1, 10000)
    samples_100k = metropolis_hastings(P, 1, 100000)

    x_vals = np.linspace(-10, 10, 1000)
    y_vals = P_norm(x_vals)
    plt.figure(1)
    plt.plot(x_vals, y_vals, 'r', label='P(x) normalised')
    
    plt.hist(samples_100k, bins=100, density=True, color='y', alpha=0.7, label='MCMC Samples (100k)')
    #plt.hist(samples_10k, bins=100, density=True, color='g', alpha=0.5, label='MCMC Samples (10k)')
    #plt.hist(samples_1k, bins=100, density=True, color='b', alpha=0.3, label='MCMC Samples (1k)')
    
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('Density')
    #plt.savefig('backend/mcmc_plots/p_x.png')
    plt.show()"""

    """# Plot density for Miyamoto-Nagai potential
    R = np.linspace(0.1, 10, 1000)
    z = np.linspace(-10, 10, 1000)
    a = 10
    b = 0.1
    R, z = np.meshgrid(R, z)
    potential = MiyamotoNagaiPotential(a=a, b=b, normalize=1)
    density = potential.dens(R, z, 0)
    plt.figure(2)
    plt.contour(R/a, z/a, density, 10)
    plt.title('Miyamoto-Nagai Density')
    plt.xlabel('R')
    plt.ylabel('z')
    plt.xlim(0, 5)
    plt.ylim(-1, 1)
    plt.savefig('backend/mcmc_plots/miyamoto_nagai_density.png')
    #plt.show()
    print('MCMC plots saved to backend/mcmc_plots')

    # Sample from Miyamoto-Nagai potential density using MCMC
    n_samples = 100000
    density = lambda R, z: potential.dens(R, z, 0)
    samples = metropolis_hastings(density, 2, n_samples)

    # Plot samples from Miyamoto-Nagai potential density
    plt.figure(3)
    xmin, xmax = 0, 5
    ymin, ymax = -1, 1
    plt.hexbin(samples[:, 0]/a, samples[:, 1]/a, gridsize=100, cmap='plasma', extent=[xmin, xmax, ymin, ymax], bins='log')
    plt.colorbar()
    plt.title('Miyamoto-Nagai Potential Samples')
    plt.xlabel('R/a')
    plt.ylabel('z/a')
    plt.xlim(0, 5)
    plt.ylim(-1, 1)
    plt.savefig('backend/mcmc_plots/miyamoto_nagai_samples_hexbin.png')
    plt.show()"""
    
    # Plot samples for Miyamoto-Nagai potential density but reduced to a sphere
    n_samples = 100000
    a = 0
    b = 10
    potential = MiyamotoNagaiPotential(a=a, b=b, normalize=1)
    density = lambda R, z: potential.dens(R, z, 0)
    samples = metropolis_hastings(density, 2, n_samples)
    plt.figure(4)
    plt.hexbin(samples[:, 0], samples[:, 1], gridsize=100, cmap='plasma', bins='log')
    plt.colorbar()
    plt.title('Miyamoto-Nagai Potential Samples (Reduced to Sphere)')
    plt.xlabel('R')
    plt.ylabel('z')
    
    plt.axis('equal')
    
    plt.show()


if __name__ == "__main__":
    main()


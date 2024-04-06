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
    # Plot samples for Miyamoto-Nagai potential density but reduced to a sphere (for testing)
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


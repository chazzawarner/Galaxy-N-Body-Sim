import numpy as np
import matplotlib.pyplot as plt
from galpy.potential import MiyamotoNagaiPotential, HernquistPotential
from tqdm import tqdm

# Metropolis-Hastings MCMC algorithm for multivariate distributions
def metropolis_hastings(target_density, dimensions, n_samples, burn_in=1000):
    n_samples += burn_in
    x0 = np.zeros(dimensions)
    xt = x0
    samples = np.zeros((n_samples, dimensions))
    for i in tqdm(range(n_samples)):
        xt_candidate = xt + np.random.normal(size=dimensions)
        acceptance_probability = min(1, target_density(*xt_candidate) / target_density(*xt))
        if np.random.uniform() < acceptance_probability:
            xt = xt_candidate
        samples[i] = xt
            
    return samples[burn_in:]


def main():
    P = lambda x: 3 * np.exp(-x*x/2) + np.exp(-(x - 4)**2/2)
    Z = 10.0261955464

    x_vals = np.linspace(-10, 10, 1000)
    y_vals = P(x_vals)
    plt.figure(1)
    plt.plot(x_vals, y_vals, 'r', label='P(x)')
    plt.legend(loc='upper right', shadow=True)
    plt.savefig('backend/mcmc_plots/p_x.png')

    # Plot density for Miyamoto-Nagai potential
    R = np.linspace(0.1, 10, 1000)
    z = np.linspace(-10, 10, 1000)
    a = 1.2
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
    print('MCMC plots saved to backend/mcmc_plots')

    # Sample from Miyamoto-Nagai potential density using MCMC
    n_samples = 1000000
    density = lambda R, z: potential.dens(R, z, 0)
    samples = metropolis_hastings(density, 2, n_samples)

    # Plot samples from Miyamoto-Nagai potential density
    plt.figure(3)
    xmin, xmax = 0, 5
    ymin, ymax = -1, 1
    plt.hexbin(samples[:, 0]/a, samples[:, 1]/a, gridsize=50, cmap='plasma', extent=[xmin, xmax, ymin, ymax], bins='log')
    plt.colorbar()
    plt.title('Miyamoto-Nagai Potential Samples')
    plt.xlabel('R/a')
    plt.ylabel('z/a')
    plt.xlim(0, 5)
    plt.ylim(-1, 1)
    plt.savefig('backend/mcmc_plots/miyamoto_nagai_samples_hexbin.png')


if __name__ == "__main__":
    main()


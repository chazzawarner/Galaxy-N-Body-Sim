import numpy as np
import matplotlib.pyplot as plt
from galpy.potential import HernquistPotential
import tqdm as tqdm

def mass_fraction_sampling(density, num_samples, total_mass=1):
    mass_fractions = np.linspace(0, 1, num_samples)
    
    
    
    
    return samples


def main():
    num_samples = 100
    hernquist = HernquistPotential(a=1, normalize=1)
    density = lambda r: hernquist.dens(r, 0)
    samples = mass_fraction_sampling(density, num_samples)
    #print(samples)
    plt.figure(1)
    plt.hist(samples, bins=100, density=True, color='r', alpha=0.7, label='Inverse Sampling')
    r = np.linspace(0, 2, 1000)
    plt.plot(r, density(r), 'r', label='Hernquist Density')
    plt.legend(loc='upper right', shadow=True)
    plt.title('Inverse Sampling of Hernquist Density')
    plt.xlabel('r')
    plt.ylabel('Density')
    plt.xlim(-0.1, 2)
    plt.savefig('backend/inverse_sampling_plots/hernquist_density.png')
    
    '''# Plot CDF of Hernquist density
    cdf, r = mass_fraction_sampling(density, num_samples)
    plt.figure(2)
    plt.plot(r, cdf, 'r', label='Hernquist CDF')
    plt.legend(loc='upper right', shadow=True)
    plt.title('Hernquist CDF')
    plt.xlabel('r')
    plt.ylabel('CDF')
    plt.savefig('backend/inverse_sampling_plots/hernquist_cdf.png')'''
    
if __name__ == '__main__':
    main()
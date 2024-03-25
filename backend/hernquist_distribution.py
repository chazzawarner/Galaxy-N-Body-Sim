import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
from latex2sympy2 import latex2sympy as l2s
from galpy.potential import HernquistPotential

mass_frac_l2s = l2s(r"M_f=\frac{r^2}{(r+a)^2}")

M_r, r, a_symbol = sp.symbols("M_r r a")
mass_frac_eq = sp.Eq(mass_frac_l2s, M_r)
radius_solutions = sp.solve(mass_frac_eq, r)
radius_solutions_func = sp.lambdify((M_r, a_symbol), radius_solutions[0])

def hernquist_mass_fraction_sampling(num_samples, a=1):
    mass_fractions = np.random.rand(num_samples)
    radius = radius_solutions_func(mass_fractions, a)
    return radius

def main():
    num_samples = 1000
    hernquist = HernquistPotential(a=1, normalize=1)
    density = lambda r: hernquist.dens(r, 0)
    samples = hernquist_mass_fraction_sampling(num_samples)
    
    plt.figure(1)
    plt.hist(samples, bins=100, density=True, color='r', alpha=0.7, label='Mass Fraction Sampling')
    r = np.linspace(0, 2, 1000)
    #plt.plot(r, density(r), 'r', label='Hernquist Density')
    plt.legend(loc='upper right', shadow=True)
    plt.title('Mass Fraction Sampling of Hernquist Density')
    plt.xlabel('r')
    plt.ylabel('Density')
    plt.xlim(-0.1, 2)
    #plt.yscale('log')
    plt.show()

    
    
    
if __name__ == '__main__':
    main()
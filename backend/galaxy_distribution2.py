import numpy as np
import matplotlib.pyplot as plt
import json
import astropy.units as u
from mcmc import metropolis_hastings
from galpy.potential import MiyamotoNagaiPotential, HernquistPotential, NFWPotential

# Galaxy generation function
class Galaxy:
    def __init__(self, name, components, num_bodies=1000000, g_const=1, pos_units=u.kpc, vel_units=u.km/u.s, mass_units=u.M_sun, density_units=u.M_sun/u.kpc**3, potential_units=u.km**2/u.s**2, time_units=u.Gyr, angle_units=u.rad):
        self.name = name
        self.components = components
        self.num_bodies = num_bodies
        self.g_const = g_const
        self.total_mass = self.get_total_mass()
        
        # Set component potentials to galpy potentials using astropy units
        for component in self.components:
            for potential in component['potentials']:
                if potential['type'] == 'Miyamoto-Nagai':
                    parameters = potential['parameters']
                    potential['galpy_potential'] = MiyamotoNagaiPotential(a=parameters['a']*pos_units, b=parameters['b']*pos_units, normalize=1)
                elif potential['type'] == 'Hernquist':
                    parameters = potential['parameters']
                    potential['galpy_potential'] = HernquistPotential(a=parameters['a']*pos_units)
                elif potential['type'] == 'NFW':
                    potential['parameters'].update({'normalise': 1})
                    potential['galpy_potential'] = NFWPotential(a=parameters['a']*pos_units)
        
        self.total_potential = self.get_total_potential()
        
        # Generate bodies from components
        for component in self.components:
            if not component['dark matter']:
                num_bodies_comp = round((self.get_component_mass(component)/self.total_mass) * self.num_bodies) # Number of bodies in component
                positions = self.generate_positions(component, num_bodies_comp)
                velocities = self.generate_velocities(positions, self.total_potential)
                component.update({'bodies': {'positions': positions, 'velocities': velocities}})
                    
    # Get total mass of galaxy from all components    
    def get_total_mass(self):
        total_mass = 0
        for component in self.components:
            if not component['dark matter']:
                for potential in component['potentials']:
                    total_mass += potential['parameters']['mass']
        return total_mass
    
    # Get mass of a single component
    def get_component_mass(self, component):
        component_mass = 0
        for potential in component['potentials']:
            component_mass += potential['parameters']['mass']
        return component_mass
    
    # Get total potential of galaxy from all components
    def get_total_potential(self):
        potentials = []
        for component in self.components:
            for potential in component['potentials']:
                potentials.append(potential['galpy_potential'])
        
        total_potential = potentials[0]
        for potential in potentials[1:]:
            total_potential += potential
        
        return total_potential
    
    # Generate positions of bodies from component
    def generate_positions(self, component, num_bodies):
        print(f'Generating {num_bodies} bodies for {component["name"]}')
        positions = np.empty((0, 3))
        for potential in component['potentials']:
            
            pot_mass = potential['parameters']['mass']
            pot_bodies = round((pot_mass/self.get_component_mass(component)) * num_bodies)
            
            if potential['type'] == 'Miyamoto-Nagai': # Miyamoto-Nagai (Axisymmetric potentials)
                density = lambda R, z: potential['galpy_potential'].dens(R, z, 0)
                samples = metropolis_hastings(density, 2, pot_bodies)
                theta = np.random.uniform(0, 2*np.pi, pot_bodies)
                
                R = samples[:, 0]
                z = samples[:, 1]
                
                x = R * np.cos(theta)
                y = R * np.sin(theta)
                
                positions = np.vstack((positions, np.column_stack((x, y, z))))
                
            else: # Hernquist or NFW (Spherical potentials) 
                density = lambda r: potential['galpy_potential'].dens(r + 1e-10, 0)
                samples = metropolis_hastings(density, 1, num_bodies)
                print(samples)
                direction = np.random.normal(size=(num_bodies, 3))
                direction /= np.linalg.norm(direction, axis=1).reshape(-1, 1)
                print(direction)
                
                positions = np.vstack((positions, samples * direction))
        
        return positions
    
    def generate_velocities(self, bodies, total_potential):
        pass
    
    def plot_rotational_vel(self):
        pass
    
    # Plot density of galaxy in a hexbin plot (as function of R and z)
    def plot_scatter_density(self):
        positions = np.empty((0, 3))
        for component in self.components:
            if 'bodies' in component:
                positions = np.vstack((positions, component['bodies']['positions']))
            
        R = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)
        
        plt.figure()
        plt.hexbin(R, positions[:, 2], gridsize=50, cmap='plasma', bins='log')
        plt.colorbar()
        plt.title(f'{self.name} Density')
        plt.xlabel('R (kpc)')
        plt.ylabel('z (kpc)')

        plt.savefig(f'backend/galaxy_plots/{self.name}_density_hexbin.png')
        
    def plot_component_scatter(self):
        plt.figure()
        for component in self.components:
            if 'bodies' in component:
                positions = component['bodies']['positions']
                R = np.sqrt(positions[:, 0]**2 + positions[:, 1]**2)
                print(component['name'])
                print(np.mean(R), np.mean(positions[:, 2]), np.std(R), np.std(positions[:, 2]))
                plt.scatter(R, positions[:, 2], label=component['name'], s=1)
        plt.legend()
        plt.title(f'{self.name} Component Scatter')
        plt.xlabel('R (kpc)')
        plt.ylabel('z (kpc)')
        #plt.xscale('log')
        plt.savefig(f'backend/galaxy_plots/{self.name}_component_scatter.png')
        
        
        
                    
    
    
    
    

def main():
    with open('backend/galaxies/generic.json', 'r') as f:
        galaxy_json = json.load(f)
    galaxy = Galaxy(galaxy_json['name'], galaxy_json['components'], num_bodies=1000000)
    print(galaxy.total_mass)
    galaxy.plot_scatter_density()
    galaxy.plot_component_scatter()
    
    
if __name__ == "__main__":
    main()
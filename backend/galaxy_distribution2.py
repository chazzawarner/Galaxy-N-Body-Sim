import numpy as np
import matplotlib.pyplot as plt
import json
import astropy.units as u
from mcmc import metropolis_hastings
from galpy.potential import MiyamotoNagaiPotential, HernquistPotential, NFWPotential, vcirc, tdyn, mass
from galaxy_intitialisation import elliptical_galaxy_gen
from tqdm import tqdm
from star_types import get_random_star

# Galaxy generation function
class Galaxy:
    def __init__(self, json_file, num_bodies=1000000, g_const=1, pos_units=u.kpc, vel_units=u.km/u.s, mass_units=u.M_sun, density_units=u.M_sun/u.kpc**3, potential_units=u.km**2/u.s**2, time_units=u.Gyr, angle_units=u.rad):
        # Load JSON file
        with open(json_file, 'r') as f:
            galaxy_json = json.load(f)
        
        self.file = json_file
        self.name = galaxy_json['name']
        self.components = galaxy_json['components']
        self.num_bodies = num_bodies
        self.g_const = g_const
        self.total_mass = self.get_total_mass()
        
        # Initialise galaxy
        self.initialise_galaxy(pos_units, mass_units)
        
    # Initialise galaxy from the provided JSON file
    def initialise_galaxy(self, pos_units, mass_units):
        # Set galpy potentials for each component
        self.set_galpy_potentials(pos_units, mass_units)
        
        # Calculate total potential of galaxy to use for rotational velocity
        self.total_potential = self.get_total_potential()
        self.rotational_velocities = []
        
        # Generate bodies from components
        for component in self.components:
            if not component['dark matter']:
                comp_mass = self.get_component_mass(component)
                num_bodies_comp = round((comp_mass/self.total_mass) * self.num_bodies) # Number of bodies in component
                
                positions = self.generate_positions(component, num_bodies_comp, comp_mass)
                velocities = self.generate_velocities(positions, self.total_potential)
                #masses = np.full(np.size(positions, axis=0), comp_mass/num_bodies_comp)
                masses = [get_random_star() for i in range(np.size(positions, axis=0))]
                component.update({'bodies': {'positions': positions, 'velocities': velocities, 'masses': masses}})
                
        # Update JSON file with galaxy data
        #self.export_galaxy(f'backend/galaxies/{self.name}.json')
    
    # Set galpy potentials for each component
    def set_galpy_potentials(self, pos_units, mass_units):
        # Set component potentials to galpy potentials using astropy units
        for component in self.components:
            for potential in component['potentials']:
                if potential['type'] == 'Miyamoto-Nagai':
                    parameters = potential['parameters']
                    potential['galpy_potential'] = MiyamotoNagaiPotential(amp=parameters['mass']*mass_units, a=parameters['a']*pos_units, b=parameters['b']*pos_units)
                elif potential['type'] == 'Hernquist':
                    parameters = potential['parameters']
                    potential['galpy_potential'] = HernquistPotential(amp=parameters['mass']*mass_units, a=parameters['a']*pos_units)
                elif potential['type'] == 'NFW':
                    potential['parameters'].update({'normalise': 1})
                    potential['galpy_potential'] = NFWPotential(amp=parameters['mass']*mass_units, a=parameters['a']*pos_units)
        
                    
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
        
        total_potential = potentials
        """for potential in potentials[1:]:
            total_potential.__add__(potential)"""
        
        print(f"Total potential: {total_potential}")
        
        return total_potential
    
    # Generate positions of bodies from component
    def generate_positions(self, component, num_bodies_comp, comp_mass):
        print(f'Generating {num_bodies_comp} bodies for component {component["name"]}')
        positions = np.empty((0, 3))
        for potential in component['potentials']:
            
            pot_mass = potential['parameters']['mass']
            print(f"Potential mass: {pot_mass}")
            pot_bodies = round((pot_mass/comp_mass) * num_bodies_comp)
            print(f"Potential bodies: {pot_bodies}")
            
            if potential['type'] == 'Miyamoto-Nagai': # Miyamoto-Nagai (Axisymmetric potentials)
                density = lambda R, z: potential['galpy_potential'].dens(R, z, 0)
                print(f"Num. bodies in potential {potential['type']}: {pot_bodies}")
                samples = metropolis_hastings(density, 2, pot_bodies)
                theta = np.random.uniform(0, 2*np.pi, pot_bodies)
                
                R = samples[:, 0]
                z = samples[:, 1]
                
                x = R * np.cos(theta)
                y = R * np.sin(theta)
                
                positions = np.vstack((positions, np.array([x, y, z]).T))
                print(positions)
                print(f"Generated {len(positions)} positions, num. bodies: {pot_bodies}")
                
                
            else: # Hernquist
                #density = lambda r: potential['galpy_potential'].dens(r + 1e-10, 0)
                #samples = metropolis_hastings(density, 1, num_bodies)
                print(f"Num. bodies in potential {potential['type']}: {pot_bodies}")
                samples = elliptical_galaxy_gen(pot_bodies, a=potential['parameters']['a'])
                
                """# Scale samples to the radius of "a" - a precautionary measure
                scale = potential['parameters']['a'] / np.percentile(samples, 99)
                print(f"Scale: {scale}")
                samples = samples * scale"""
                
                # Reject samples outside the maximum radius of the galaxy
                samples = samples[samples < self.get_galaxy_max_radius()]
                
                samples = samples.reshape(-1, 1)
                #print(samples)
                
                
                # Generate random directions for bodies
                direction = np.random.normal(size=(np.size(samples), 3))
                direction /= np.linalg.norm(direction, axis=1).reshape(-1, 1)
                #print(direction)
                
                positions = np.vstack((positions, samples * direction))
                #print(f"Positions: {positions}")
        
        return positions
    
    def get_galaxy_max_radius(self):
        max_radius = 0
        for component in self.components:
            for potential in component['potentials']:
                if potential['type'] == 'Miyamoto-Nagai':
                    a = potential['parameters']['a']
                    b = potential['parameters']['b']
                    if a > max_radius:
                        max_radius = a
                    if b > max_radius:
                        max_radius = b
                else:
                    if potential['parameters']['a'] > max_radius:
                        max_radius = potential['parameters']['a']
        return max_radius
    
    def generate_velocities(self, positions, total_potential):
        velocities = np.empty((0, 3))
        
        print('Calculating rotational velocities')
        for position in tqdm(positions):
            radius = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
            
            # Calculate the rotational speed of the body
            #v_circ = total_potential.vcirc(radius)
            v_circ = vcirc(total_potential, radius * u.kpc) #!!
            
            # Find the unit vector perpendicular to the direction of the body
            velocities_direction = np.array([-position[1], position[0], 0])
            velocities_direction /= np.linalg.norm(velocities_direction)
            
            # Calculate the rotational velocity of the body
            velocity = velocities_direction * v_circ
            velocities = np.vstack((velocities, velocity))
            self.rotational_velocities.append([v_circ, radius])
        
        return velocities
    
    def plot_rotational_vel(self):
        r = np.linspace(1e-10, self.get_galaxy_max_radius(), 1000)
        
        plt.figure()
        
        for component in self.components:
            total_component_potential = component['potentials'][0]['galpy_potential']
            for potential in component['potentials'][1:]:
                total_component_potential.__add__(potential['galpy_potential'])
            v_circ_component = total_component_potential.vcirc(r * u.kpc)
            #print(f"v_circ for {component['name']}: {v_circ_component}")
            
            plt.plot(r, v_circ_component, label=component['name'], linestyle='--')
        
        #v_circ_total = self.total_potential.vcirc(r * u.kpc)
        v_circ_total = vcirc(self.total_potential, r * u.kpc)
        plt.plot(r, v_circ_total, label='Total')
        plt.legend()
        plt.title(f'{self.name} Rotational Velocity')
        plt.xlabel('r (kpc)')
        plt.ylabel('v_circ (km/s)')
        plt.savefig(f'backend/galaxy_plots/{self.name}_rotational_velocity.png')
        plt.show()
            
            
                
        
    
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
        plt.show()
        
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
        #plt.savefig(f'backend/galaxy_plots/{self.name}_component_scatter.png')
        plt.show()
    
    # Plot galaxy as a 3D scatter plot and side on views (2x2 plot)
    def plot_full_galaxy(self):
        positions = np.empty((0, 3))
        for component in self.components:
            if 'bodies' in component:
                positions = np.vstack((positions, component['bodies']['positions']))
        
        fig, axs = plt.subplots(2, 2, figsize=(8, 8))  # Increase the figure size
        
        ## 3D plot
        ax = fig.add_subplot(2, 2, 1, projection='3d')
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], s=1)
        ax.set_title("3D view")
        
        ## Remove empty graph/axes for the 3D plot
        fig.delaxes(axs[0, 0])
        
        ## 2D Side view
        ax = axs[0, 1]
        #ax.scatter(positions[:, 1], positions[:, 2], s=1)
        hb = ax.hexbin(positions[:, 1], positions[:, 2], gridsize=50, cmap='plasma', bins='log')
        hb.set_facecolor('black')  # Set the hexbin plot background color to black
        ax.set_title("2D Side view")
        ax.set_aspect('equal', 'box')
        
        ## 2D Top view
        ax = axs[1, 0]
        #ax.scatter(positions[:, 0], positions[:, 1], s=1)
        hb = ax.hexbin(positions[:, 0], positions[:, 1], gridsize=50, cmap='plasma', bins='log')
        hb.set_facecolor('black')  # Set the hexbin plot background color to black
        ax.set_title("2D Top view")
        ax.set_aspect('equal', 'box')
        
        ## 2D Front view
        ax = axs[1, 1]
        #ax.scatter(positions[:, 0], positions[:, 2], s=1)
        ax.hexbin(positions[:, 0], positions[:, 2], gridsize=50, cmap='plasma', bins='log', facecolor='black')
        ax.set_title("2D Front view")
        ax.set_aspect('equal', 'box')
        
        plt.tight_layout()
        plt.show()
        
    # Calculate dynamical time of galaxy
    def get_t_dyn(self):
        t_dyn = tdyn(self.total_potential, self.get_galaxy_max_radius() * u.kpc)
        return t_dyn
    
    # Export galaxy 
    def export_galaxy(self, filename):
        components_copy = self.components.copy()
        for component in components_copy:
            for potential in component['potentials']:
                del potential['galpy_potential']
                
        with open(filename, 'w') as f:
            json.dump({'name': self.name, 'components': components_copy}, f, indent=4)
        
    # Return the galaxy bodies as a dictionary
    def get_galaxy(self):
        bodies = {
            "position": np.empty((0, 3)),
            "velocity": np.empty((0, 3)),
            "mass": np.empty(0)
        }
        
        for component in self.components:
            if 'bodies' in component:
                bodies['position'] = np.vstack((bodies['position'], component['bodies']['positions']))
                bodies['velocity'] = np.vstack((bodies['velocity'], component['bodies']['velocities']))
                bodies['mass'] = np.hstack((bodies['mass'], component['bodies']['masses']))
            
        return bodies
        
        
        
                    
    
    
    
    

def main():
    galaxy = Galaxy('backend/galaxies/basic_galaxy.json', num_bodies=10000)
    print(galaxy.total_mass)
    galaxy.plot_scatter_density()
    #galaxy.plot_component_scatter()
    galaxy.plot_rotational_vel()
    galaxy.plot_full_galaxy()
    """print(f"Dynamical time: {galaxy.get_t_dyn()} Gyr")
    #print(f"To return: {galaxy.get_galaxy()}")
    bodies = galaxy.get_galaxy()
    print(f"Total mass of bodies: {np.sum(bodies['mass'])}")
    print(f"Predicted total mass: {galaxy.total_mass}")
    print(f"Requested number of bodies: {galaxy.num_bodies}")
    print(np.size(bodies['mass']))
    print(np.size(bodies['position'], axis=0))
    print(np.shape(bodies['position']))
    print(np.size(bodies['velocity'])) 
    print(np.shape(bodies['velocity']))  """ 
    
if __name__ == "__main__":
    main()
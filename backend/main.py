import numpy as np
import csv
from nbody import NBody
from galaxy_distribution import Galaxy
from mergers import init_merger
import time
from datetime import datetime
import cProfile
import astropy.units as u
from astropy.constants import G

timesteps = 100 
timestep = (1 * u.Gyr).to(u.s).value / 1e5 
#print(f"Time step: {timestep}")

def main():
    print("Initialising simulation...")

    two_galaxies = True
    if two_galaxies:
        # Define the galaxy JSONs
        galaxy_jsons = ['backend/galaxies/andromeda.json', 'backend/galaxies/milkyway.json']
        
        # Calculate inital position and velocities
        distance_between_galaxies = 15e3 # parsecs
        approach_speed = 100 # km/s
        
        # Set positions
        galaxy_positions = [np.array([distance_between_galaxies/2, 1e3, 0]), np.array([-distance_between_galaxies/2, -1e3, 0])]
        
        # Calculate unit vectors for velocities
        galaxy_vec_12 = galaxy_positions[1] - galaxy_positions[0]
        galaxy_vec_12 = galaxy_vec_12 / np.linalg.norm(galaxy_vec_12)
        galaxy_vec_21 = galaxy_positions[0] - galaxy_positions[1]
        galaxy_vec_21 = galaxy_vec_21 / np.linalg.norm(galaxy_vec_21)
        
        # Set velocities
        galaxy_velocities = [galaxy_vec_12 * approach_speed/2, galaxy_vec_21 * approach_speed/2]
        
        # Get bodies
        bodies = init_merger(galaxy_jsons, galaxy_positions, galaxy_velocities, total_num_bodies=2000, check_csv=False)
        
        
    else:
        """galaxy = Galaxy('backend/galaxies/basic_galaxy.json', num_bodies=1000)
        bodies = galaxy.get_galaxy()"""
        
        bodies_csv = np.loadtxt('data/Jason the Galaxy.csv', delimiter=',', skiprows=1)
        bodies = {
            "mass": bodies_csv[:, 1],
            "position": bodies_csv[:, 2:5],
            "velocity": bodies_csv[:, 5:],
            "galaxy": np.full(len(bodies_csv), 0) # All from the same galaxy
        }
        
        # Convert postions to parsecs
        bodies["position"] = bodies["position"] * 1e3
    
    
    
    g_const = G.to_value(u.pc**3 / (u.M_sun * u.s**2))
    print(f"G: {g_const}")
    nbody = NBody(bodies)
    
    print("Total mass of bodies: ", np.sum(nbody.bodies["mass"]))
    
    # Run the simulation for timesteps iterations, saving time, mass, and position to data/output.csv
    print("Running simulation...")

    with open("data/output.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'mass', 'pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z', 'galaxy'])
    
        # Write the initial state
        for i in range(len(nbody.bodies["mass"])):
            writer.writerow([0, nbody.bodies["mass"][i], nbody.bodies["position"][i][0], nbody.bodies["position"][i][1], nbody.bodies["position"][i][2], nbody.bodies["velocity"][i][0], nbody.bodies["velocity"][i][1], nbody.bodies["velocity"][i][2], nbody.bodies["galaxy"][i]])
            
        # Write the subsequent states, updating the bodies at each step
        for i in range(timesteps-1):
            print(f'Step {i+1}/{timesteps}')
            nbody.update(timestep, remove_outliers=False)
            
            for j in range(len(nbody.bodies["mass"])):  # Write the time, mass, and position
                writer.writerow([i+1, nbody.bodies["mass"][j], nbody.bodies["position"][j][0], nbody.bodies["position"][j][1], nbody.bodies["position"][j][2], nbody.bodies["velocity"][j][0], nbody.bodies["velocity"][j][1], nbody.bodies["velocity"][j][2], nbody.bodies["galaxy"][j]])
        
    print("Simulation complete")
    
    
    
if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"Time taken in total: {end_time - start_time} seconds")
    #cProfile.run('main()')

import numpy as np
import csv
from nbody import NBody
from galaxy_distribution import Galaxy
import time
from datetime import datetime
import cProfile
import astropy.units as u
from astropy.constants import G

num_bodies = 100
timesteps = 100 
timestep = (1 * u.Gyr).to(u.s).value / 1e6
#print(f"Time step: {timestep}")
bounding_box = 46.56e3 # Size of the bounding box in parsecs (i.e. the diameter of the galaxy)
bh_mass = 3.5e9 # Mass of the black hole in solar masses
g_const = 4.3009e-3

def main():
    # Define 100 bodies with random mass, position, and velocity
    print("Creating bodies...")
    """bodies = []
    for _ in range(num_bodies):
        mass = get_random_star()
        centre_vec = np.array([bounding_box/2, bounding_box/2])
        direction = np.random.uniform(-1, 1, size=(2,))
        direction /= np.linalg.norm(direction)
        distance = np.random.uniform(0, 20e3) + 2500 #(np.random.exponential(scale=1) * bounding_box / 2) + 250 #0.1 * bounding_box
        position = centre_vec + direction * distance
        centre_to_pos = position - centre_vec
        c_to_p_perp = np.array([-centre_to_pos[1], centre_to_pos[0]])
        c_to_p_perp_norm = c_to_p_perp / np.linalg.norm(c_to_p_perp)
        distance_to_centre = np.linalg.norm(centre_to_pos)
        speed = np.sqrt( (g_const * bh_mass) / distance_to_centre ) #(1 - distance_to_centre / (bounding_box / 2)) * (np.linalg.norm(centre_to_pos) / (0.5 * bounding_box))
        print(f"Speed: {speed}")
        velocity = c_to_p_perp_norm * speed
        bodies.append(Body(mass, position, velocity))
    #bodies[0].mass = bh_mass  # Set the mass of the first body to be much larger than the others
    #bodies[0].position = np.array([bounding_box/2, bounding_box/2])  # Set the position of the first body to the origin
    #bodies[0].velocity = np.array([0.0, 0.0])  # Set the velocity of the first body to zero"""
    #bodies = galaxy_generation(num_bodies, bounding_box)
    
    print(f"Time step: {timestep} s")

    two_galaxies = True
    if two_galaxies:
        distance_between_galaxies = 5e3 # parsecs
        approach_velocity = 110 # km/s
        
        galaxy_1 = Galaxy('backend/galaxies/basic_galaxy.json', num_bodies=500)
        galaxy_1_bodies = galaxy_1.get_galaxy()
        galaxy_1_radius = np.percentile(np.linalg.norm(galaxy_1_bodies["position"], axis=1), 95)
        galaxy_1_bodies["position"] = galaxy_1_bodies["position"][:, :2] # In units of kpc
        galaxy_1_bodies["velocity"] = galaxy_1_bodies["velocity"][:, :2] # In units of km/s
        galaxy_1_bodies["position"] *= 1e3 # Into parsecs
        
        galaxy_2 = Galaxy('backend/galaxies/basic_galaxy.json', num_bodies=500)
        galaxy_2_bodies = galaxy_2.get_galaxy()
        galaxy_2_radius = np.percentile(np.linalg.norm(galaxy_2_bodies["position"], axis=1), 95)
        galaxy_2_bodies["position"] = galaxy_2_bodies["position"][:, :2] # In units of kpc
        galaxy_2_bodies["velocity"] = galaxy_2_bodies["velocity"][:, :2] # In units of km/s
        galaxy_2_bodies["position"] *= 1e3 # in parsecs instead of kpc
        
        total_distance = galaxy_1_radius + galaxy_2_radius + distance_between_galaxies
        galaxy_1_bodies["position"] += np.array([total_distance/2, 1e3])
        galaxy_2_bodies["position"] -= np.array([total_distance/2, 0])
        
        galaxy_1_bodies["velocity"] -= np.array([approach_velocity, 0])
        galaxy_2_bodies["velocity"] += np.array([approach_velocity, 0])
        
        bodies = {
            "mass": np.concatenate((galaxy_1_bodies["mass"], galaxy_2_bodies["mass"])),
            "position": np.concatenate((galaxy_1_bodies["position"], galaxy_2_bodies["position"])),
            "velocity": np.concatenate((galaxy_1_bodies["velocity"], galaxy_2_bodies["velocity"]))
        }
        
        
    else:
        galaxy = Galaxy('backend/galaxies/basic_galaxy.json', num_bodies=1000)
        bodies = galaxy.get_galaxy()
        
        # Only consider the x and y coordinates
        bodies["position"] = bodies["position"][:, :2] # In units of kpc
        bodies["velocity"] = bodies["velocity"][:, :2] # In units of km/s
        
        # Convert postions to parsecs
        bodies["position"] = bodies["position"] * 1e3
    
    
    
    g_const = G.to_value(u.pc**3 / (u.M_sun * u.s**2))
    print(f"G: {g_const}")
    nbody = NBody(bodies)
    
    print("Total mass of bodies: ", np.sum(nbody.bodies["mass"]))
    #print("Total mass of bodies (excluding black hole): ", np.sum(nbody.bodies["mass"][1:]))
    
    # Run the simulation for timesteps iterations, saving time, mass, and position to data/output.csv
    print("Running simulation...")
    
    now = datetime.now()
    date_time_str = str(int(now.timestamp()))
    filename = f'data/output_n{num_bodies}_ts{timesteps}_dt{timestep}_{date_time_str}.csv'

    #with open(filename, 'w', newline='') as f: # Sort out PyScript first
    with open("data/output.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'mass', 'pos_x', 'pos_y', 'vel_x', 'vel_y'])  # Write the header

        #writer.writerow([0, bodies[0].mass, bodies[0].position[0], bodies[0].position[1]])  # Write the initial state
        """for body in nbody.bodies:
                writer.writerow([0, body.mass, body.position[0], body.position[1]])
        
        for i in range(timesteps-1):
            print(f'Step {i+1}/{timesteps}')
            nbody.update(1)  # Update the bodies
            for body in nbody.bodies:
                writer.writerow([i+1, body.mass, body.position[0], body.position[1]])  # Write the time, mass, and position"""
                
        # Write the initial state
        for i in range(len(nbody.bodies["mass"])):
            writer.writerow([0, nbody.bodies["mass"][i], nbody.bodies["position"][i][0], nbody.bodies["position"][i][1], nbody.bodies["velocity"][i][0], nbody.bodies["velocity"][i][1]])
            
        # Write the subsequent states, updating the bodies at each step
        for i in range(timesteps-1):
            print(f'Step {i+1}/{timesteps}')
            nbody.update(timestep, remove_outliers=False)
            
            for j in range(len(nbody.bodies["mass"])):  # Write the time, mass, and position
                writer.writerow([i+1, nbody.bodies["mass"][j], nbody.bodies["position"][j][0], nbody.bodies["position"][j][1], nbody.bodies["velocity"][j][0], nbody.bodies["velocity"][j][1]])
        
    print("Simulation complete")
    
    
    
if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"Time taken in total: {end_time - start_time} seconds")
    #cProfile.run('main()')

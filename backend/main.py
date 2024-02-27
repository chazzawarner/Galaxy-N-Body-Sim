import numpy as np
import csv
from numpy import random
from nbody import NBody
from galaxy_distribution import galaxy_generation
import time
from datetime import datetime
import cProfile

num_bodies = 100
timesteps = 100
timestep = 1
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
    bodies = galaxy_generation(num_bodies, bounding_box)
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
            nbody.update(timestep)
            
            for j in range(len(nbody.bodies["mass"])):  # Write the time, mass, and position
                writer.writerow([i+1, nbody.bodies["mass"][j], nbody.bodies["position"][j][0], nbody.bodies["position"][j][1], nbody.bodies["velocity"][j][0], nbody.bodies["velocity"][j][1]])
        
    print("Simulation complete")
    
    
    
if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"Time taken in total: {end_time - start_time} seconds")
    #cProfile.run('main()')

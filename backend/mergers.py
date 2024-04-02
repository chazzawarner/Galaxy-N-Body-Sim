import numpy as np
import matplotlib.pyplot as plt
from galaxy_distribution import Galaxy
import json

# Merger initialisation function to place galaxies in initial positions
# Useful only when the galaxies aren't already created
def init_merger(galaxy_jsons, galaxy_positions, galaxy_velocities, total_num_bodies=1000, check_csv=True, csv_dir='backend/galaxies/'):
    # Initialise the galaxies
    galaxies = []
    """for galaxy in galaxy_jsons:
        num_bodies = np.round(total_num_bodies/len(galaxy_jsons))
        galaxies.append(Galaxy(galaxy, num_bodies=num_bodies))"""
        
    galaxies_bodies = []    
    
    # Check if the CSV file exists
    if check_csv:
        total_csv_bodies = 0
        for i in range(len(galaxy_jsons)):
            with open(galaxy_jsons[i]) as f:
                galaxy_json = json.load(f)
                print("Attempting to load " + csv_dir + galaxy_json["name"] + '.csv')
                try:
                    galaxy_csv = np.loadtxt(csv_dir + galaxy_json["name"] + '.csv', delimiter=',', skiprows=1)
                    total_csv_bodies += len(galaxy_csv)
                    
                    galaxy_bodies = {
                        "mass": galaxy_csv[:, 1],
                        "position": galaxy_csv[:, 2:5],
                        "velocity": galaxy_csv[:, 5:]
                    }
                    
                    # Convert postions to parsecs
                    galaxy_bodies["position"] = galaxy_bodies["position"] * 1e3
                    
                    galaxies_bodies.append(galaxy_bodies)
                    
                except:
                    print("CSV file not found")
                    break
                
        if not np.isclose(total_csv_bodies, total_num_bodies, atol=1):
            print(f"CSV files contain {total_csv_bodies} bodies, but {total_num_bodies} bodies are required")
            galaxies_bodies = []
            
    # If the CSV files don't exist, create the galaxies
    if len(galaxies_bodies) == 0:
        print("Creating galaxies...")
        
        galaxies_masses = []
        
        for i in range(len(galaxy_jsons)):
            with open(galaxy_jsons[i]) as f:
                galaxy_json = json.load(f)
                galaxy_mass = 0
                for component in galaxy_json["components"]:
                    if not component["dark matter"]:
                        for potential in component["potentials"]:
                            galaxy_mass += potential["parameters"]["mass"]
                galaxies_masses.append(galaxy_mass)
        
        for i in range(len(galaxy_jsons)):
            #num_bodies = np.round(total_num_bodies/len(galaxy_jsons))
            
            num_bodies = np.round(total_num_bodies * galaxies_masses[i] / np.sum(galaxies_masses))
            
            galaxies.append(Galaxy(galaxy_jsons[i], num_bodies=num_bodies))
        
        for galaxy in galaxies:
            temp = galaxy.get_galaxy()
            temp["position"] *= 1e3  # Convert positions to parsecs
            galaxies_bodies.append(temp)
    
    
    # Set the initial positions and velocities of the galaxies
    for i in range(len(galaxies_bodies)):
        galaxies_bodies[i]["position"] += galaxy_positions[i]
        galaxies_bodies[i]["velocity"] += galaxy_velocities[i]
        
    # Combine the galaxies
    bodies = {
        "mass": np.concatenate([galaxy["mass"] for galaxy in galaxies_bodies]),
        "position": np.concatenate([galaxy["position"] for galaxy in galaxies_bodies]),
        "velocity": np.concatenate([galaxy["velocity"] for galaxy in galaxies_bodies]),
        "galaxy": np.concatenate([np.full(len(galaxy["mass"]), i) for i, galaxy in enumerate(galaxies_bodies)]) # Assign each body to its galaxy
    }
    
    return bodies

def main():
    # Define the galaxy JSONs
    galaxy_jsons = ['backend/galaxies/paper_galaxy.json', 'backend/galaxies/paper_galaxy.json']
    
    # Calculate inital position and velocities
    distance_between_galaxies = 15e3 # parsecs
    approach_speed = 110 # km/s
    
    # Set positions
    galaxy_positions = [np.array([distance_between_galaxies/2, 2e3, 1e3]), np.array([-distance_between_galaxies/2, -2e3, -1e3])]
    
    # Calculate unit vectors for velocities
    galaxy_vec_12 = galaxy_positions[1] - galaxy_positions[0]
    galaxy_vec_12 = galaxy_vec_12 / np.linalg.norm(galaxy_vec_12)
    galaxy_vec_21 = galaxy_positions[0] - galaxy_positions[1]
    galaxy_vec_21 = galaxy_vec_21 / np.linalg.norm(galaxy_vec_21)
    
    # Set velocities
    galaxy_velocities = [galaxy_vec_12 * approach_speed/2, galaxy_vec_21 * approach_speed/2]
    
    # Get bodies
    bodies = init_merger(galaxy_jsons, galaxy_positions, galaxy_velocities)
    
    # Plot the initial positions of the galaxies
    fig, axs = plt.subplots(ncols=2, nrows=2)
    axs[0, 0].scatter(bodies["position"][:, 0], bodies["position"][:, 1], s=1)
    axs[0, 0].set_title("Top-down view (x-y)")
    axs[0, 0].set_aspect('equal')
    
    axs[0, 1].scatter(bodies["position"][:, 0], bodies["position"][:, 2], s=1)
    axs[0, 1].set_title("Side view (x-z)")
    axs[0, 1].set_aspect('equal')
    
    axs[1, 0].scatter(bodies["position"][:, 1], bodies["position"][:, 2], s=1)
    axs[1, 0].set_title("Front view (y-z)")
    axs[1, 0].set_aspect('equal')
    
    plt.show()

if __name__ == '__main__':
    main()
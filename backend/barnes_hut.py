import numpy as np
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool
from tqdm import tqdm

# Define a helper function for multiprocessing
def compute_force_helper(args):
    return compute_force(*args)



# Define the Node class
## Nodes can be either a region with >1 body (internal), a region with 1 body (external), or an empty region
class Node:
    def __init__(self, centre, size):
        self.centre = centre
        self.size = size
        self.centre_of_mass = 0
        self.mass = 0
        self.children = None
        
    def create_children(self):
        self.children = {
            "NW": Node(self.centre + np.array([-self.size/4, self.size/4]), self.size/2),
            "NE": Node(self.centre + np.array([self.size/4, self.size/4]), self.size/2),
            "SW": Node(self.centre + np.array([-self.size/4, -self.size/4]), self.size/2),
            "SE": Node(self.centre + np.array([self.size/4, -self.size/4]), self.size/2)
        }
    
    def in_quadrant(self, position):
        return (self.centre[0] - self.size/2 <= position[0] < self.centre[0] + self.size/2 and 
                self.centre[1] - self.size/2 <= position[1] < self.centre[1] + self.size/2)
    
    def insert(self, body):
        epsilon = 1e-3  # Adjust this value as needed
        body["position"] += epsilon * (np.random.rand(2) - 0.5)  # Add a small random offset to the position to avoid divide by zero errors
        
        # If the node is empty, insert the body
        if self.children is None and self.mass == 0:
            self.centre_of_mass = body["position"]
            self.mass = body["mass"]
        
        # If internal node, update the centre of mass and mass, recursively insert the body into the appropriate quadrant
        elif self.children is not None:
            # Update the centre of mass and mass
            self.mass += body["mass"]
            self.centre_of_mass = (self.centre_of_mass * self.mass + body["position"] * body["mass"]) / (self.mass + body["mass"])
            
            # Insert the body into the appropriate quadrant
            for child in self.children:
                if self.children[child].in_quadrant(body["position"]):
                    self.children[child].insert(body)
        
        # If external node, subdivide the region further by creating four children, then recursively insert both the body and the body already there, then update the centre of mass and mass
        elif self.children is None and self.mass > 0:
            self.create_children()
            old_body = {
                "mass": self.mass,
                "position": self.centre_of_mass
            }
            
            for child in self.children:
                if self.children[child].in_quadrant(old_body["position"]):
                    self.children[child].insert(old_body)
                if self.children[child].in_quadrant(body["position"]):
                    self.children[child].insert(body)
            
            
        else:
            print("Error: Node is neither empty, internal, nor external")
        
# Define quadtree function
def create_quadtree(bodies):
    pos_x_max = max(body[0] for body in bodies["position"])
    pos_x_min = min(body[0] for body in bodies["position"])
    pos_y_max = max(body[1] for body in bodies["position"])
    pos_y_min = min(body[1] for body in bodies["position"])
    
    bounding_box = max(pos_x_max - pos_x_min, pos_y_max - pos_y_min)
    centre = np.array([pos_x_min + bounding_box/2, pos_y_min + bounding_box/2])
    
    quadtree = Node(centre, bounding_box)
    for i in range(len(bodies["mass"])):
        body = {
            "mass": bodies["mass"][i],
            "position": bodies["position"][i]
        }
        
        quadtree.insert(body)
    
    return quadtree

def compute_force(node, body, theta, smoothing_length, g_const):
    force = np.array([0.0, 0.0])
    epsilon = 1e-5
    
    # Making sure the body is not interacting with itself
    if not np.array_equal(node.centre_of_mass, body["position"]):
    # If node is external (and not the body itself), calculate the force exerted by the body, adding it to the body's net force
        if node.children is None and node.mass > 0:
            r = node.centre_of_mass - body["position"]
            distance = np.linalg.norm(r) + epsilon
            force += (g_const * node.mass * body["mass"] / (smoothing_length**2 + distance**2)**(3/2)) * r
        
        # Calculate s/d, if < theta, use the node as a single body, calculate the force exerted by the body, adding it to the body's net force
        elif np.linalg.norm(node.centre - body["position"]) / node.size < theta:
            r = node.centre_of_mass - body["position"]
            distance = np.linalg.norm(r) + epsilon
            #force += (g_const * node.mass * body["mass"] / distance**3) * r
            force += (g_const * node.mass * body["mass"] / (smoothing_length**2 + distance**2)**(3/2)) * r
        
        # If s/d > theta and node is internal, run procedure recursively on each of the node's children
        elif node.children is not None:
            for child in node.children:
                force += compute_force(node.children[child], body, theta, smoothing_length, g_const)
                
        # If node is empty, do nothing
        else:
            pass

    return force



# Use the Barnes-Hut algorithm to calculate the forces on each body
def compute_all_forces(bodies, theta=0.5, g_const=4.3009e-3):
    
    
    
    
    time_start = time.time()
    quadtree = create_quadtree(bodies)
    time_end = time.time()
    print(f"Time taken to create quadtree: {time_end - time_start} seconds")
    
    forces = np.zeros((len(bodies["mass"]), 2))
    
    smoothing_length = 1.5 * (len(bodies["mass"]) ** (-0.44)) # Optimal smoothing length for Barnes-Hut algorithm using Hernquist
    
    time.start = time.time()
    
    """for i in range(len(bodies["mass"])):
        body = {
            "mass": bodies["mass"][i],
            "position": bodies["position"][i]
        }
        
        forces[i] = compute_force(quadtree, body, theta)
        #print(f"Forces on body {i}: {forces[i]}")"""
        
    print("Computing forces...")
    
    # Create a multiprocessing Pool
    with Pool() as p:
        # Create a list of arguments for each body
        args = [(quadtree, {"mass": bodies["mass"][i], "position": bodies["position"][i]}, theta, smoothing_length, g_const) for i in range(len(bodies["mass"]))]
        # Use the Pool's map function to compute the forces in parallel
        forces = list(tqdm(p.imap(compute_force_helper, args), total=len(args)))
        
    time.end = time.time()
    print(f"Time taken to compute forces: {time.end - time.start} seconds")
    
    
    return forces

def main():
    num_bodies = 100
    timesteps = 10
    timestep = 1
    bbox = [0, 100] # Bounding box, min/max x and y coordinates
    
    bodies = {
        "mass": [np.random.uniform(0.1, 10) for _ in range(num_bodies)],
        "position": np.random.uniform(0, 100, size=(num_bodies, 2)),
        "velocity": np.zeros((num_bodies, 2))
    }
    
    #quadtree = create_quadtree(bodies, np.array([bbox[0]/2, bbox[0]/2]), bbox[1] - bbox[0])
    
    # Plot the quadtree, drawing squares for each node
    def plot_qt(quadtree):
        if quadtree.children is not None:
            for child in quadtree.children:
                plot_qt(quadtree.children[child])
        
        plt.plot([quadtree.centre[0] - quadtree.size/2, quadtree.centre[0] + quadtree.size/2], [quadtree.centre[1] - quadtree.size/2, quadtree.centre[1] - quadtree.size/2], 'k-', linewidth=0.5)
        plt.plot([quadtree.centre[0] - quadtree.size/2, quadtree.centre[0] + quadtree.size/2], [quadtree.centre[1] + quadtree.size/2, quadtree.centre[1] + quadtree.size/2], 'k-', linewidth=0.5)
        plt.plot([quadtree.centre[0] - quadtree.size/2, quadtree.centre[0] - quadtree.size/2], [quadtree.centre[1] - quadtree.size/2, quadtree.centre[1] + quadtree.size/2], 'k-', linewidth=0.5)
        plt.plot([quadtree.centre[0] + quadtree.size/2, quadtree.centre[0] + quadtree.size/2], [quadtree.centre[1] - quadtree.size/2, quadtree.centre[1] + quadtree.size/2], 'k-', linewidth=0.5)
        
    #plot_qt(quadtree)
    
    

    
    forces = compute_all_forces(bodies)
    
    print(forces)
    
    updated_bodies = {
        "mass": np.copy(bodies["mass"]),
        "position": np.copy(bodies["position"]),
        "velocity": np.copy(bodies["velocity"])
    }
    
    updated_bodies["velocity"] += forces * timestep / updated_bodies["mass"].reshape(-1, 1)
    updated_bodies["position"] += updated_bodies["velocity"] * timestep
    
    print(bodies["position"] == updated_bodies["position"])
    
    plt.scatter(bodies["position"][:, 0], bodies["position"][:, 1], s=2)
    plt.scatter(updated_bodies["position"][:, 0], updated_bodies["position"][:, 1])
    plt.xlim(bbox)
    plt.ylim(bbox)
    plt.axis('equal')
    
    plt.show()
    

if __name__ == '__main__':
    main()
    
import numpy as np
import time
from barnes_hut import compute_all_forces as barnes_hut

g_const = 4.3009e-3  # Gravitational constant in units of parsecs * (km/s)^2 / solar_mass
    
# Define compute_force method with brute force algorithm
def brute_force(self):
    # Compute the forces between all pairs of bodies using the brute force algorithm
    # NOT IMPLEMENTED FOR NEW VECTORISED VERSION
    time_start = time.time()
    for i in range(self.n):
        for j in range(self.n):
            if i != j:
                r = self.bodies[j].position - self.bodies[i].position
                r_norm = np.linalg.norm(r)
                force = r * self.bodies[i].mass * self.bodies[j].mass * g_const / r_norm**3
                self.forces[i] += force
    time_end = time.time()
    print(f"Time taken to compute forces: {time_end - time_start} seconds (brute force method)")


    
# Defining the Body class
class Body:
    def __init__(self, mass, position, velocity):
        self.mass = mass
        self.position = position
        self.velocity = velocity

    def update(self, force, dt):
        # Update the velocity and position of the body based on the applied force and time step
        self.velocity += force * dt / self.mass
        self.position += self.velocity * dt
    
    def __str__(self):
        return f'Body(mass={self.mass}, position={self.position}, velocity={self.velocity})'
    
    def __repr__(self):
        return str(self)    
 
# Defining the NBody class
class NBody:
    def __init__(self, bodies):
        self.bodies = bodies # List of Body objects
        """{
            "mass": [body.mass for body in bodies],
            "position": [body.position for body in bodies],
            "velocity": [body.velocity for body in bodies]
            } #bodies """
        self.n = len(bodies)
        self.forces = np.zeros((self.n, 2))
    
    def compute_force(self, method="barnes_hut"):
        if method == "brute_force":
            forces = brute_force(self.bodies)
        elif method == "barnes_hut":
            forces = barnes_hut(self.bodies, g_const=g_const, theta=0.3)
        return forces
           
    # Remove outlying bodies if body is more than x standard deviations from the mean
    def remove_outliers(self, std_devs_lim=5):
        distances = np.linalg.norm(self.bodies["position"], axis=1)
        mean_position = np.mean(self.bodies["position"], axis=0)
        mean_distance = np.mean(distances)
        std_dev_distance = np.std(distances)
        
        bodies_to_remove = []
        
        for i in range(len(self.bodies["position"])):
            # SORT THIS OUT
            if np.linalg.norm(self.bodies["position"][i] - mean_position) > std_devs_lim * std_dev_distance:
                bodies_to_remove.append(i)
                
        """for i in bodies_to_remove:
            self.bodies["position"] = np.delete(self.bodies["position"], i, axis=0)
            self.bodies["velocity"] = np.delete(self.bodies["velocity"], i, axis=0)
            self.bodies["mass"] = np.delete(self.bodies["mass"], i)
            self.n -= 1
            print(f"Body {i} removed as outlier")"""

        if len(bodies_to_remove) > 0:
            self.bodies["position"] = np.delete(self.bodies["position"], bodies_to_remove, axis=0)
            self.bodies["velocity"] = np.delete(self.bodies["velocity"], bodies_to_remove, axis=0)
            self.bodies["mass"] = np.delete(self.bodies["mass"], bodies_to_remove)
            self.n -= len(bodies_to_remove)
            print(f"Bodies {bodies_to_remove} removed as outliers")
                
        
                    
    def update(self, dt, remove_outliers=True):
        # Update the bodies' positions and velocities based on the computed forces and time step
        if remove_outliers:
            self.remove_outliers()
        
        if self.forces is None:
            self.forces = np.zeros_like(self.bodies["position"])  # Initialise self.forces if it's None

        self.forces = self.compute_force(method="barnes_hut")
        #print(f"Forces: {self.forces}")
        
        epsilon = 1e-8
        
        self.bodies["velocity"] += self.forces * dt / (np.array(self.bodies["mass"]).reshape(-1, 1) + epsilon) # Update the velocity of the bodies
        #self.bodies["velocity"][0] = np.array([0.0, 0.0]) # Set BH velocity to zero
        self.bodies["position"] += self.bodies["velocity"] * dt # Update the position of the bodies
        
        print("Black hole force: ", self.forces[0])
        print("Black hole position: ", self.bodies["position"][0])
        print("Black hole mass: ", self.bodies["mass"][0])
        print("Black hole velocity: ", self.bodies["velocity"][0])
        
            
    def __str__(self):
        return f'NBody(bodies={self.bodies})'
    
    def __repr__(self):
        return str(self)
import numpy as np
import matplotlib.pyplot as plt

# Define the star types and their mass ranges (Harvard spectral classification)
star_types = {
    "O": {
        "proportion": 0.00003,
        "mass-range": [16, 100],
    },
    "B": {
        "proportion": 0.12,
        "mass-range": [2.1, 16],
    },
    "A":{ 
        "proportion": 0.61,
        "mass-range": [1.4, 2.1],
    },
    "F": {
        "proportion": 3.0,
        "mass-range": [1.04, 1.4],
    },
    "G": {
        "proportion": 7.6,
        "mass-range": [0.8, 1.04],
    },
    "K": {
        "proportion": 12,
        "mass-range": [0.45, 0.8],
    },
    "M": {
        "proportion": 76,
        "mass-range": [0.08, 0.45],
    },
}

total_proportion = sum(data["proportion"] for data in star_types.values())
#print(f"Total proportion: {total_proportion}")
proportion_scale = 100 / total_proportion

# Return a random star mass from the star types dictionary
def get_random_star():
    rand = np.random.rand() * 100
    #print(f"Random number: {rand}")
    cum_prop = 0
    for data in reversed(star_types.values()):
        if rand <= (data["proportion"] + cum_prop) * proportion_scale:
            #print(f"Star type: {data['mass-range']} ")
            mass = np.random.uniform(data["mass-range"][0], data["mass-range"][1])
            #print(f"Random star mass: {mass}")
            return mass
        
        cum_prop += data["proportion"]
        
# Get the average star mass of the star types dictionary
def get_average_star_mass():
    return sum(data["proportion"] * np.mean(data["mass-range"]) for data in star_types.values()) / 100
        
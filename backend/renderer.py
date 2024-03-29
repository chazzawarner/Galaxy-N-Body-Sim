import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# Render merger data to an animation
def render_merger(filename, output_filename="output.gif"):
    # Load the data, removing the header
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    
    # Split the data into time, mass, and position
    num_timesteps = np.max(data[:, 0])
    num_timesteps = int(num_timesteps)
    timesteps = np.empty(num_timesteps, dtype=object)
    print(f"Number of timesteps: {num_timesteps}")
    print(timesteps)
    
    for timestep in range(num_timesteps):
        timestep_data = data[data[:, 0] == timestep]
        masses = timestep_data[:, 1]
        positions = timestep_data[:, 2:]
        
        timesteps[timestep] = {
            "masses": masses,
            "positions": positions
        }
    
    # Set up the figure and axis
    fig, ax = plt.subplots()
    max_x = np.max(np.abs(timesteps[0]["positions"][:, 0]))
    max_y = np.max(np.abs(timesteps[0]["positions"][:, 1]))
    max_pos = np.max([max_x, max_y])
    ax.set_xlim(-max_pos, max_pos)
    ax.set_ylim(-max_pos, max_pos)
    
    ax.set_aspect('equal')
    
    # Animate the data
    def animate(timestep):
        ax.clear()
        ax.set_xlim(-max_pos, max_pos)
        ax.set_ylim(-max_pos, max_pos)
        ax.set_aspect('equal')
        
        masses = timesteps[timestep]["masses"]
        positions = timesteps[timestep]["positions"]
        
        ax.scatter(positions[:, 0], positions[:, 1], s=masses)
        ax.text(-max_pos, max_pos, f"Time: {timestep}")
        
        print(f"Animating timestep {timestep}")
        
        return ax
        
        
    print("Animating...")
    ani = animation.FuncAnimation(fig, animate, frames=num_timesteps, interval=100)
    
    print("Saving animation...")
    
    # Save animation as gif
    writer = animation.PillowWriter(fps=10, bitrate=1800)
    ani.save(output_filename, writer=writer)
        
    print("Animation saved to ", output_filename)
    
def main():
    render_merger("data/output.csv", output_filename="data/output.gif")

if __name__ == '__main__':
    main()
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# Render merger data to an animation
def render_merger(filename, output_filename="output.gif", trails=False):
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
    
    # Set up the figure and axis - 2x2 grid
    fig, axs = plt.subplots(nrows=2, ncols=2)
    max_x = np.max(np.abs(timesteps[0]["positions"][:, 0]))
    max_y = np.max(np.abs(timesteps[0]["positions"][:, 1]))
    max_z = np.max(np.abs(timesteps[0]["positions"][:, 2]))
    max_pos_top = max(max_x, max_y)
    max_pos_side = max(max_x, max_z)
    max_pos_front = max(max_y, max_z)
    
    # Plot a straight line between current and previous position for each body
    def plot_trails(plot, components, timestep):
        if timestep != 0:
            for i in range(len(timesteps[timestep]["positions"])):
                prev_pos = timesteps[timestep - 1]["positions"][i]
                pos = timesteps[timestep]["positions"][i]
                
                axs[plot[0], plot[1]].plot([prev_pos[components[0]], pos[components[0]]], [prev_pos[components[1]], pos[components[1]]], 'k-', linewidth=0.5, zorder=1, alpha=0.2)
        
    
    # Animate the data
    def animate(timestep):
        masses = timesteps[timestep]["masses"]
        positions = timesteps[timestep]["positions"]
        
        # Top-down view (x-y)
        axs[0,0].clear()
        axs[0,0].set_xlim(-max_pos_top, max_pos_top)
        axs[0,0].set_ylim(-max_pos_top, max_pos_top)
        axs[0,0].set_aspect('equal')
        
        if trails:
            plot_trails([0, 0], [0, 1], timestep)
        
        axs[0,0].scatter(positions[:, 0], positions[:, 1], s=masses)
        
        axs[0,0].set_title("Top-down view (x-y)", fontsize=12)
        
        # Side view (x-z)
        axs[0,1].clear()
        axs[0,1].set_xlim(-max_pos_side, max_pos_side)
        axs[0,1].set_ylim(-max_pos_side, max_pos_side)
        axs[0,1].set_aspect('equal')
        
        if trails:
            plot_trails([0, 1], [0, 2], timestep)
        
        axs[0,1].scatter(positions[:, 0], positions[:, 2], s=masses)
        axs[0,1].set_title("Side view (x-z)", fontsize=12)
        
        # Front view (y-z)
        axs[1,0].clear()
        axs[1,0].set_xlim(-max_pos_front, max_pos_front)
        axs[1,0].set_ylim(-max_pos_front, max_pos_front)
        axs[1,0].set_aspect('equal')
        
        if trails:
            plot_trails([1, 0], [1, 2], timestep)
        
        axs[1,0].scatter(positions[:, 1], positions[:, 2], s=masses)
        axs[1,0].set_title("Front view (y-z)", fontsize=12)
        
        # Put the timestep in the corner
        axs[1,1].clear()
        axs[1,1].axis('off')
        axs[1,1].text(0.5, 0.5, f"Time: {timestep}", fontsize=12, ha='center')
        
        print(f"Animating timestep {timestep}")
        
        return axs
        
        
    print("Animating...")
    ani = animation.FuncAnimation(fig, animate, frames=num_timesteps, interval=100)
    
    print("Saving animation...")
    
    # Save animation as gif
    writer = animation.PillowWriter(fps=10, bitrate=1800)
    ani.save(output_filename, writer=writer)
        
    print("Animation saved to ", output_filename)
    
def main():
    render_merger("data/output.csv", output_filename="data/output.gif", trails=True)

if __name__ == '__main__':
    main()
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib.offsetbox
from matplotlib.lines import Line2D
from star_types import star_types

# Class to render scalebar
class AnchoredHScaleBar(matplotlib.offsetbox.AnchoredOffsetbox):
    """ size: length of bar in data units
        extent : height of bar ends in axes units """
    def __init__(self, size=1, extent = 0.03, label="", loc=2, ax=None,
                 pad=0.4, borderpad=0.5, ppad = 0, sep=2, prop=None, 
                 frameon=True, linekw={}, **kwargs):
        if not ax:
            ax = plt.gca()
        trans = ax.get_xaxis_transform()
        size_bar = matplotlib.offsetbox.AuxTransformBox(trans)
        line = Line2D([0,size],[0,0], **linekw)
        vline1 = Line2D([0,0],[-extent/2.,extent/2.], **linekw)
        vline2 = Line2D([size,size],[-extent/2.,extent/2.], **linekw)
        size_bar.add_artist(line)
        size_bar.add_artist(vline1)
        size_bar.add_artist(vline2)
        txt = matplotlib.offsetbox.TextArea(label)
        self.vpac = matplotlib.offsetbox.VPacker(children=[size_bar,txt],  
                                 align="center", pad=ppad, sep=sep) 
        matplotlib.offsetbox.AnchoredOffsetbox.__init__(self, loc, pad=pad, 
                 borderpad=borderpad, child=self.vpac, prop=prop, frameon=frameon,
                 **kwargs)

# Render merger data to an animation
def render_merger(filename, output_filename="output.gif", trails=False, colours=["blue", "red"], depth_alpha=True):
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
        positions = timestep_data[:, 2:5]
        galaxy = timestep_data[:, 8]
        
        timesteps[timestep] = {
            "masses": masses,
            "positions": positions,
            "galaxy": galaxy
        }
    
    # Set up the figure and axis - 2x2 grid
    fig, axs = plt.subplots(nrows=2, ncols=2)
    max_x = np.max(np.abs(timesteps[0]["positions"][:, 0]))
    max_y = np.max(np.abs(timesteps[0]["positions"][:, 1]))
    max_z = np.max(np.abs(timesteps[0]["positions"][:, 2]))
    max_pos_top = max(max_x, max_y)
    max_pos_side = max(max_x, max_z)
    max_pos_front = max(max_y, max_z)
    
    fig.set_size_inches(8, 8)
    fig.set_tight_layout(True)
    
    # Plot a straight line between current and previous position for each body
    def plot_trails(plot, components, timestep):
        if timestep != 0:
            for i in range(len(timesteps[timestep]["positions"])):
                prev_pos = timesteps[timestep - 1]["positions"][i]
                pos = timesteps[timestep]["positions"][i]
                
                axs[plot[0], plot[1]].plot([prev_pos[components[0]], pos[components[0]]], [prev_pos[components[1]], pos[components[1]]], 'k-', linewidth=0.5, zorder=1, alpha=0.2)
    
    max_mass = 16 # solar masses, taken from star_types.py, minimum mass of an O-type star
    min_mass = 0.4

    # Animate the data
    def animate(timestep):
        masses = timesteps[timestep]["masses"]
        positions = timesteps[timestep]["positions"]
        galaxy = timesteps[timestep]["galaxy"]
        
        # Top-down view (x-y)
        axs[0,0].clear()
        axs[0,0].set_xlim(-max_pos_top, max_pos_top)
        axs[0,0].set_ylim(-max_pos_top, max_pos_top)
        axs[0,0].set_aspect('equal')
        
        if trails:
            plot_trails([0, 0], [0, 1], timestep)
            
        if depth_alpha:
            max_z = np.max(np.abs(positions[:, 2]))
            min_z = np.min(np.abs(positions[:, 2]))
            alpha = np.interp(positions[:, 2], [min_z, max_z], [0.1, 0.8])
        
        axs[0,0].scatter(positions[:, 0], positions[:, 1], s=masses, c=[colours[int(g)] for g in galaxy], zorder=2, alpha=alpha)
        
        axs[0,0].set_title("Top-down view (x-y)", fontsize=12)
        
        # Side view (x-z)
        axs[0,1].clear()
        axs[0,1].set_xlim(-max_pos_side, max_pos_side)
        axs[0,1].set_ylim(-max_pos_side, max_pos_side)
        axs[0,1].set_aspect('equal')
        
        if trails:
            plot_trails([0, 1], [0, 2], timestep)
            
        if depth_alpha:
            max_y = np.max(np.abs(positions[:, 1]))
            min_y = np.min(np.abs(positions[:, 1]))
            alpha = np.interp(positions[:, 1], [min_y, max_y], [0.1, 0.8])
        
        axs[0,1].scatter(positions[:, 0], positions[:, 2], s=masses, c=[colours[int(g)] for g in galaxy], zorder=2, alpha=alpha)
        axs[0,1].set_title("Side view (x-z)", fontsize=12)
        
        # Front view (y-z)
        axs[1,0].clear()
        axs[1,0].set_xlim(-max_pos_front, max_pos_front)
        axs[1,0].set_ylim(-max_pos_front, max_pos_front)
        axs[1,0].set_aspect('equal')
        
        if trails:
            plot_trails([1, 0], [1, 2], timestep)
            
        if depth_alpha:
            max_x = np.max(np.abs(positions[:, 0]))
            min_x = np.min(np.abs(positions[:, 0]))
            alpha = np.interp(positions[:, 0], [min_x, max_x], [0.1, 0.8])
        
        axs[1,0].scatter(positions[:, 1], positions[:, 2], s=masses, c=[colours[int(g)] for g in galaxy], zorder=2, alpha=alpha)
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
    
# Render specific frames from the merger data
def render_specific(filename, output_filename, frames, trails=False, colours=["blue", "red"], depth_alpha=True, vertical=True):
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
        positions = timestep_data[:, 2:5]
        galaxy = timestep_data[:, 8]
        
        timesteps[timestep] = {
            "masses": masses,
            "positions": positions,
            "galaxy": galaxy
        }
    
    # Find the maximum position in each dimension
    max_x = np.max(np.abs(timesteps[0]["positions"][:, 0]))
    max_y = np.max(np.abs(timesteps[0]["positions"][:, 1]))
    max_pos_top = max(max_x, max_y)
    
    # Define trails function
    def plot_trails(plot_index, components, frame):
        if frame != 0:
            for i in range(len(timesteps[frame]["positions"])):
                prev_pos = timesteps[frame - 1]["positions"][i]
                pos = timesteps[frame]["positions"][i]
                
                axs[plot_index].plot([prev_pos[components[0]], pos[components[0]]], [prev_pos[components[1]], pos[components[1]]], 'k-', linewidth=0.5, zorder=1, alpha=0.1)
    
    # Define depth alpha function
    def depth_alpha(positions, component):
        max_component = np.max(np.abs(positions[:, component]))
        min_component = np.min(np.abs(positions[:, component]))
        alpha = np.interp(positions[:, component], [min_component, max_component], [0.1, 0.8])
        
        return alpha
    
    # Initialise the figure and subplots
    num_frames = len(frames)
    plot_size = 2.5 # inches
    if vertical:
        fig, axs = plt.subplots(nrows=int(num_frames/2), ncols=2)
        fig.set_size_inches(plot_size * 2, plot_size * int(num_frames/2))
    else:
        fig, axs = plt.subplots(nrows=2, ncols=int(num_frames/2))
        fig.set_size_inches(plot_size * int(num_frames/2), plot_size * 2)
    axs = axs.flatten()
    
    for a in axs:
        a.set_aspect('equal')
        a.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)
        a.set_xlim(-max_pos_top, max_pos_top)
        a.set_ylim(-max_pos_top, max_pos_top)
        a.set_xticks(np.linspace(a.get_xlim()[0], a.get_xlim()[1], num=5))
        a.set_yticks(np.linspace(a.get_ylim()[0], a.get_ylim()[1], num=5))
        a.set_xticklabels([])
        a.set_yticklabels([])
    
    # Plot the frames
    for i, frame in enumerate(frames):
        # Get the data for the frame
        masses = timesteps[frame]["masses"]
        positions = timesteps[frame]["positions"]
        galaxy = timesteps[frame]["galaxy"]
        
        # If trails are enabled, plot them
        if trails:
            plot_trails(i, [0, 1], frame)
            
        # If depth alpha is enabled, calculate it
        if depth_alpha:
            alpha = depth_alpha(positions, 2)
        else:
            alpha = 0.8
        
        # Plot the frame from the top down (x-y)
        axs[i].scatter(positions[:, 0], positions[:, 1], s=masses, c=[colours[int(g)] for g in galaxy], zorder=2, alpha=alpha)
        axs[i].text(0.05, 0.9, f"t = {frame}", transform=axs[i].transAxes, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
        
    # Add scalebar
    barsize = 10e3 # parsecs
    ob = AnchoredHScaleBar(size=barsize, label=f"{int(barsize/1e3)} kpc", loc=4, frameon=False, pad=0.6, sep=5, linekw=dict(color="black"))
    axs[0].add_artist(ob)
        
 
    # Adjust the layout
    fig.tight_layout()
    fig.subplots_adjust(hspace=0, wspace=0)
    
    plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
    plt.show()
        

# Main function    
def main():
    render_merger("data/output.csv", output_filename="data/output.gif", trails=True)
    
    #render_specific("data/output.csv", "backend/galaxy_plots/???", [0, 1] , trails=True, depth_alpha=False, vertical=False)

if __name__ == '__main__':
    main()
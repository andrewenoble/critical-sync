import scipy as sp
import pylab as pl
from matplotlib import animation


# Define plotting and animation parameters.
L = 256 # Defines the size of the simulated 2D network of N = L x L nodes.

m_min = -0.5 # Defines a range of synchronization order parameter values 
m_max = 0.5  # for the density plot in each frame of the animation.

num_frames = 200
fps = 4
bitrate = 1600
blit = True


# Load the simulation output.
m_matrix = []

for i in xrange(num_frames):     
                                                        
    m_matrix.append(sp.loadtxt('../simulation_output/' 
        + 'm_' + str(i) + '.txt'))
    

# Create a base figure that will be updated to create each frame in the 
# animation sequence.
fig = pl.figure(figsize=(8., 8.))
ax = fig.add_axes([0, 0, 1, 1])

# The base figure will be the first frame (i = 0) in the animation sequence.    
i = 0

# Plot a density plot of the synchronization order parameter field.
print min(m_matrix[i]), max(m_matrix[i])
im = ax.imshow(m_matrix[i].reshape(L,L), cmap=pl.get_cmap('Greys'), vmin=m_min, 
    vmax=m_max, interpolation='nearest')
pl.axis('off')    


# Define a dummy function to initialize the animation.
def init(): return im,    

# Define a function called to update the base figure, creating the ith frame
# of the animation.
def animate(i):
    
    im.set_data(m_matrix[i].reshape(L,L))

    return im,

# Generate the animation.
anim = animation.FuncAnimation(fig, animate, init_func=init, 
    frames=num_frames, blit=blit)


# Save the animation.
anim.save('critical_sync.mp4', fps=fps, bitrate=bitrate,
    extra_args=['-vcodec', 'libx264'])

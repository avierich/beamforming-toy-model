import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


frames = 1000  # Number of frames in the animation
grid_size = (128, 128)  # Size of heatmap grid



OMEGA = 0.5

NUM_TRANSMITTERS = 4

THETA = np.pi*(180+30)/180.0
TX_OFFSET = 0

transmitters = [(x,TX_OFFSET) for x in range(int(grid_size[0]/4), int(3*grid_size[0]/4), int(grid_size[0]/(4*NUM_TRANSMITTERS)))]


def tx(tx_pos, pos, theta, t):
    # Compute radial distance from tx
    r = np.sqrt(np.power(tx_pos[0]-pos[0], 2) + np.power(tx_pos[1]-pos[1],2))
    beta_x = -1.0*OMEGA*(tx_pos[0]-grid_size[0]/4)*np.sin(theta)
    beta_x2 = -1.0*OMEGA*(tx_pos[0]-grid_size[0]/4)*np.sin(np.pi*(180-20)/180.0)
    return 0.5*np.exp(1.0j*(OMEGA*t - r*OMEGA + beta_x)) + 0.5*np.exp(1.0j*(OMEGA*t - r*OMEGA + beta_x2))

field = np.empty(grid_size)

def generate_data(theta, t):
    for index, output, in np.ndenumerate(field):
        field[index] = 0
        for transmitter in transmitters:
            field[index] += np.real(tx(transmitter,index,theta, t))
    return field

def generate_tx_points():
    points_x = []
    points_y = []
    for transmitter in transmitters:
        points_x.append(transmitter[0])
        points_y.append(transmitter[1])
    return points_x, points_y

# Initialize figure
fig, ax = plt.subplots()
heatmap_data = generate_data(THETA,0)
heatmap = ax.imshow(heatmap_data, vmin=-1.0*NUM_TRANSMITTERS, vmax=1.0*NUM_TRANSMITTERS, cmap="coolwarm", interpolation='sinc', animated=True)
tx_points = generate_tx_points()
scatter = ax.scatter(tx_points[1], tx_points[0], marker = 'v',edgecolors='black')

# Update function for animation
def update(frame):
    # heatmap.set_array(generate_data(THETA + 2*np.pi*frame/frames, frame))  # Update heatmap data]
    heatmap.set_array(generate_data(THETA, frame))  # Update heatmap data]

    # scatter.set_offsets(tx_points[1], tx_points[0])
    return [heatmap, scatter]

# Create animation
ani = animation.FuncAnimation(fig, update, frames=frames, interval=30, blit=True)

# Save as GIF (optional)
# ani.save("animated_heatmap.gif", writer="pillow")

plt.show()

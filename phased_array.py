import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


frames = 400  # Number of frames in the animation
RESOLUTION = 512
grid_size = (RESOLUTION, RESOLUTION)  # Size of heatmap grid




NOISE = 0.000001

NUM_TRANSMITTERS = 16

THETA = np.pi*(180+20)/180.0
THETA_B = np.pi*(180-10)/180.0

TX_OFFSET = 0

# transmitters = [(x,TX_OFFSET) for x in range(int(grid_size[0]/4)+int(grid_size[0]/32), int(3*grid_size[0]/4), int(grid_size[0]/(4*NUM_TRANSMITTERS)))]

OMEGA = 0.25

LAMBDA = 2*np.pi*1.0/OMEGA

tx_spacing = LAMBDA/2.0

transmitters = [(n*tx_spacing + RESOLUTION/2 - NUM_TRANSMITTERS*tx_spacing/2+tx_spacing/2, TX_OFFSET) for n in range(0,NUM_TRANSMITTERS)]


# Precompute transmitter to cell distances
tx_distances=[]

for tx in transmitters:
    distances = np.zeros(grid_size)
    for index, cell, in np.ndenumerate(distances):
        r = np.sqrt(np.power(tx[0]-index[0], 2) + np.power(tx[1]-index[1],2))
        distances[index] = r
    tx_distances.append(distances)



def tx(r, beta_x, t):
    # Compute radial distance from tx
    return np.exp(1.0j*(OMEGA*t - r*OMEGA + beta_x))/NUM_TRANSMITTERS


N_SAMPLES = 1000
def rx_probe(rx_pos, tx_theta):
    power = 0
    t = 0
    for t in range(N_SAMPLES):
        signal_point = 0
        for idx, transmitter in enumerate(transmitters):
            beta_x = -1.0*OMEGA*(transmitter[0]-grid_size[0]/4)*np.sin(tx_theta+np.pi)
            signal_point += tx(tx_distances[idx][rx_pos], beta_x, t)
        power += np.abs(signal_point**2)
    return power/N_SAMPLES
        



field = np.empty(grid_size)



## New generate data
def generate_data(theta, t):
    field = np.zeros(grid_size)
    for idx, transmitter in enumerate(transmitters):
        beta_x = -1.0*OMEGA*(transmitter[0]-grid_size[0]/4)*np.sin(theta+np.pi)
        field += np.real(np.exp(1.0j*(OMEGA*t - tx_distances[idx]*OMEGA + beta_x))/NUM_TRANSMITTERS)
    return field

def generate_tx_points():
    points_x = []
    points_y = []
    for transmitter in transmitters:
        points_x.append(transmitter[0])
        points_y.append(transmitter[1])
    return points_x, points_y

# Initialize figure
fig, (ax, ax2) = plt.subplots(1,2,figsize=(14, 6))
heatmap_data = generate_data(THETA,0)
heatmap = ax.imshow(heatmap_data, vmin=-1.0, vmax=1.0, cmap="coolwarm", interpolation='sinc', animated=True)
tx_points = generate_tx_points()
scatter = ax.scatter(tx_points[1], tx_points[0], marker = 'v',edgecolors='black')
rx_scatter = ax.scatter([100,100], [23,23], marker = 'v',edgecolors='black', color = 'orange')

annotationA = ax.annotate(
    'Rx A', xy=(1,0), xytext=(-1,0),
)
annotationB = ax.annotate(
    'Rx B', xy=(1,0), xytext=(-1,0),
)

START_ANGLE = np.pi*(-60)/180.0 
ANGLE_DELTA = np.pi*(120)/180.0

sum_rate_plot, = ax2.plot(np.linspace(START_ANGLE, START_ANGLE+ANGLE_DELTA, frames), np.zeros(frames))
ax2.set_ylim([0,30])
ax2.set_xlim([0,120])
ax2.set_title("Sum Rate")
ax2.set_ylabel("Rate")
ax2.set_xlabel("$\Delta\\theta$ (Degrees)")

fig.tight_layout()


angles = []
sum_rates = []
# Update function for animation
def update(frame):
    # heatmap.set_array(generate_data(THETA + 2*np.pi*frame/frames, frame))  # Update heatmap data]

    print("Frame: "+str(frame))

    theta_a = np.pi*(60)/180.0
    theta_b = START_ANGLE + ANGLE_DELTA*frame/frames

    # rx a position:
    rx_pos_a = (int(grid_size[0]/2-(RESOLUTION/2)*np.sin(theta_a)), int((RESOLUTION/2)*np.cos(theta_a))+TX_OFFSET)
    rx_pos_b = (int(grid_size[0]/2-(RESOLUTION/2)*np.sin(theta_b)), int((RESOLUTION/2)*np.cos(theta_b))+TX_OFFSET)
    print(rx_pos_a)

    heatmap.set_array(0.5*generate_data(theta_a, frame)+0.5*generate_data(theta_b, frame))  # Update heatmap data]



    # snir = rx_probe((23,100),THETA)/rx_probe((23,100),THETA_B)
    # print(np.log2(1+SNIR))


    # snir_b = 20*np.log10(rx_probe(rx_pos_b,theta_b)/rx_probe(rx_pos_b,theta_a))
    # snir_a = 20*np.log10(rx_probe(rx_pos_a,theta_a)/rx_probe(rx_pos_a,theta_b))



    rate_a = np.log2(1+rx_probe(rx_pos_b,theta_b)/(NOISE+rx_probe(rx_pos_b,theta_a)))
    rate_b = np.log2(1+rx_probe(rx_pos_b,theta_b)/(NOISE+rx_probe(rx_pos_b,theta_a)))

    print("Theta a: " + str(180*theta_a/np.pi) + " Theta b: " + str(180*theta_b/np.pi))
    # print("Power aa: " + str(rx_probe(rx_pos_a,theta_a)) + " Power ab: " + str(rx_probe(rx_pos_a,theta_b)))
    print("RX B Pos: " + str(rx_pos_b))
    print("Power bb: " + str(rx_probe(rx_pos_b,theta_b)) + " Power ba: " + str(rx_probe(rx_pos_b,theta_a)))
    # print(rx_pos_a)

    # print(snir_a)
    sum_rates.append(rate_a+rate_b)
    angles.append(180*(theta_a-theta_b)/np.pi)


    sum_rate_plot.set_data(angles, sum_rates)

    # print("Probe pos a: "+str(rx_probe(rx_pos_a,theta_b)))
    # print("Probe pos b: "+str(rx_probe(rx_pos_b,theta_b)))


    rx_scatter.set_offsets((np.c_[[rx_pos_a[1],rx_pos_b[1]], [rx_pos_a[0],rx_pos_b[0]]]))

    annotationA.set_position((rx_pos_a[1]-8,rx_pos_a[0]+12))
    annotationA.xy = (rx_pos_a[1],rx_pos_a[0])

    
    annotationB.set_position((rx_pos_b[1]-8,rx_pos_b[0]+12))
    annotationB.xy = (rx_pos_b[1],rx_pos_b[0])

    return [heatmap, scatter, rx_scatter, annotationA, annotationB, sum_rate_plot]

# Create animation
ani = animation.FuncAnimation(fig, update, frames=frames, interval=1, blit=False)

# Save as GIF (optional)
# ani.save("animated_heatmap_b.gif", writer="pillow")

plt.show()

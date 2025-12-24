import struct
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

data_debug = False

####################################################################################################
# BINARY FILE READING
####################################################################################################
struct_format = 'didddddddd'  # double int double...
struct_size = struct.calcsize(struct_format)

script_dir = os.path.dirname(os.path.abspath(__file__))
filename = os.path.join(script_dir, "global_data.bin")
file_size = os.path.getsize(filename)

num_records = file_size // struct_size
print(f"Data read from {num_records} records")

# read binary file - returns separate arrays for each data element
def readBodyPos(filename):
    ts = []
    idx = []
    x = []
    y = []
    vx = []
    vy = []
    ax = []
    ay = []
    fx = []
    fy = []

    with open(filename, "rb") as file:
        while True:
            data = file.read(struct_size)
            if len(data) < struct_size:
                break
            t, body_index, x_pos, y_pos, v_x, v_y, a_x, a_y, f_x, f_y = struct.unpack(struct_format, data)
            ts.append(t)
            idx.append(body_index)
            x.append(x_pos)
            y.append(y_pos)
            vx.append(v_x)
            vy.append(v_y)
            ax.append(a_x)
            ay.append(a_y)
            fx.append(f_x)
            fy.append(f_y)

    return ts, idx, x, y, vx, vy, ax, ay, fx, fy

####################################################################################################
# PLOTTING
####################################################################################################
# global variables
timestamps, body_indices, x, y, vx, vy, ax, ay, fx, fy = readBodyPos(filename)

if data_debug:
    for i in range(len(timestamps)):
        print(f"Body {body_indices[i]}: ({x[i]:.3f}, {y[i]:.3f}) at t={timestamps[i]:.3f}")

# get unique body indices
unique_body_indices = set(body_indices)

fig = plt.figure(figsize=(16, 8))
gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1, 1], height_ratios=[1, 1])

ax1 = fig.add_subplot(gs[:, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 1])

# plot for each body
for body_idx in sorted(unique_body_indices):
    # get x, y positions and timestamps for this body
    body_timestamps = []
    body_x_positions = []
    body_y_positions = []
    for i in range(len(body_indices)):
        if body_indices[i] == body_idx:
            body_timestamps.append(timestamps[i])
            body_x_positions.append(x[i])
            body_y_positions.append(y[i])

    # subplot 1
    ax1.plot(body_x_positions, body_y_positions, label=f'Body {body_idx}', marker='o', markersize=0.1)

    # subplot 2
    ax2.plot(body_timestamps, body_x_positions, label=f'Body {body_idx}')

    # subplot 3
    ax3.plot(body_timestamps, body_y_positions, label=f'Body {body_idx}')

# configure subplot 1
ax1.set_xlabel('X Position')
ax1.set_ylabel('Y Position')
ax1.set_title('Orbital Trajectories')
ax1.legend()
ax1.grid(True)
ax1.axis('equal')

# configure subplot 2
ax2.set_xlabel('Time')
ax2.set_ylabel('X Position')
ax2.set_title('X Position vs Time')
ax2.legend()
ax2.grid(True)

# configure subplot 3
ax3.set_xlabel('Time')
ax3.set_ylabel('Y Position')
ax3.set_title('Y Position vs Time')
ax3.legend()
ax3.grid(True)

plt.tight_layout()
plt.show()

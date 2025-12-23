import struct
import os
import matplotlib.pyplot as plt

####################################################################################################
# BINARY FILE READING
####################################################################################################
struct_format = 'didd'  # double int double double
struct_size = struct.calcsize(struct_format)

script_dir = os.path.dirname(os.path.abspath(__file__))
filename = os.path.join(script_dir, "body_pos_data.bin")
file_size = os.path.getsize(filename)

num_records = file_size // struct_size
print(f"File contains {num_records} records")

# read binary file - returns list of (timestamp, body_index, x, y) tuples
def readBodyPos(filename):
    records = []
    with open(filename, "rb") as file:
        while True:
            data = file.read(struct_size)
            if len(data) < struct_size:
                break
            records.append(struct.unpack(struct_format, data))
    return records

data_list = readBodyPos(filename)

for t, body_idx, x, y in data_list:
    print(f"Body {body_idx}: ({x:.3f}, {y:.3f}) at t={t:.3f}")

####################################################################################################
# PLOTTING
####################################################################################################
# get unique body indices
body_indices = set(record[1] for record in data_list)

# plot orbital trajectories
plt.figure(figsize=(10, 10))
for body_idx in sorted(body_indices):
    # get x, y positions for this body from data_list
    x_positions = [x for t, idx, x, y in data_list if idx == body_idx]
    y_positions = [y for t, idx, x, y in data_list if idx == body_idx]
    plt.plot(x_positions, y_positions, label=f'Body {body_idx}', marker='o', markersize=2)

plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Orbital Trajectories')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

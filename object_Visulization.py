import numpy as np
import matplotlib.pyplot as plt

# Load data from the .txt file
data_file = 'data.txt'  # Replace with your file path
points = []
normals = []

with open(data_file, 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) == 6:  # Assuming format "x y z nx ny nz"
            points.append([float(parts[0]), float(parts[1]), float(parts[2])])
            normals.append([float(parts[3]), float(parts[4]), float(parts[5])])

points = np.array(points)
normals = np.array(normals)

points = points[0:5000]
normals = normals[0:5000]

# Plotting
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot points
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b', marker='o', label='Points')

# Plot normals as vectors
scale_factor = 0.1  # Adjust the length of the normal vectors for visualization
for i in range(len(points)):
    ax.plot([points[i, 0], points[i, 0] + normals[i, 0]*scale_factor],
            [points[i, 1], points[i, 1] + normals[i, 1]*scale_factor],
            [points[i, 2], points[i, 2] + normals[i, 2]*scale_factor],
            color='r')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Points and Normal Vectors Visualization')

ax.legend(['Normals', 'Points'])

plt.show()


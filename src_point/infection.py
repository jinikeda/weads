import numpy as np
from scipy.spatial import KDTree

# Define the positions of the nodes (x, y coordinates)
# Define the range for x and y
x = np.arange(6)
y = np.arange(6)

# Create a meshgrid
xv, yv = np.meshgrid(x, y)

# Stack arrays in sequence horizontally (column wise)
node_positions = np.column_stack((xv.flatten(), yv.flatten()))

print(node_positions)

node_states = np.array([[0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0], [
                       0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
node_states = node_states.flatten()
# node_positions = np.array([
#     [0, 0], [1, 1], [2, 2], [3, 3], [4, 4],
#     [5, 5], [6, 6], [7, 7], [8, 8], [9, 9]
# ])

# Define the initial states of the nodes (0: not infected, 1: infected)
#node_states = np.array([1, 0, 0, 0, 0, 0, 1, 0, 0, 0])

# Define the infection distance threshold
infection_distance = 2.0

# Create a KDTree for efficient neighbor search
tree = KDTree(node_positions)

# Find all nodes within the infection distance of the infected node(s)
# Get the indices of infected nodes here node_state a tuple
infected_nodes = np.where(node_states == 1)[0]
for node in infected_nodes:
    neighbors = tree.query_ball_point(
        node_positions[node],
        infection_distance,
        p=2.0)  # p=2 (circle) for Euclidean distance
    for neighbor in neighbors:
        if neighbor != node:  # Exclude the node itself
            node_states[neighbor] = 1

print("Node positions:")
print(node_positions)
print("Node states after infection spread:")
print(node_states)

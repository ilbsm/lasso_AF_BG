import numpy as np

with open('align/_representatives.pim.pim', 'r') as f:
    lines = f.read().split('\n')[6:-1]
names = []
matrix_distances = []
for line in lines:
    line = line.split()
    name = line[1]
    names.append(name)
    line_distances = [float(x) for x in line[2:]]
    matrix_distances.append(line_distances)
matrix_distances = np.array(matrix_distances)
print(matrix_distances)
print(matrix_distances.shape)

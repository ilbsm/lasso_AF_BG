import os
from collections import defaultdict

folder = 'clusters'
lasso_types = defaultdict(list)
files = os.listdir(folder)
for infile in files:
    with open('{}/{}'.format(folder, infile), 'r') as f:
        cluster = infile.split('.')[0]
        lassos = set([])
        for line in f.readlines():
            lasso = line.strip().split()[1]
            lassos.add(lasso)
        for lasso in lassos:
            lasso_types[lasso].append(cluster)

with open('lassos_clusters.txt', 'w') as f:
    for lasso in sorted(lasso_types.keys()):
        f.write('{} {}\n'.format(lasso, ','.join(lasso_types[lasso])))

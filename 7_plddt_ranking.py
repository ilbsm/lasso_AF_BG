import numpy as np
import os

records = []
infiles = sorted(os.listdir('clusters2'))
for infile in infiles:
    if infile[0] == '.':
        continue
    with open('clusters2/{}'.format(infile), 'r') as f:
        name = '{:10}'.format(infile.split('.')[0].split('_')[1])
        c1 = []
        c2 = []
        mean_plddt_local = []
        min_plddt_local = []
        lassos = set([])
        for line in f.readlines():
            uniprot_id, lasso, bridge, *plddts = line.strip().split()
            if bridge != '-':
                plddts = [float(x) for x in plddts]
                c1.append(plddts[0])
                c2.append(plddts[1])
                mean_plddt_local.append(plddts[2])
                min_plddt_local.append(plddts[3])
                lassos.add(lasso)
        c1 = min(c1)
        c2 = min(c2)
        mean_plddt_local = round(np.mean(mean_plddt_local),1)
        min_plddt_local = min(min_plddt_local)
        lassos = ';'.join(sorted(lassos))
        score = int(round(c1+c2+min_plddt_local, 0))
        record = [score, name, c1, c2, mean_plddt_local, min_plddt_local, lassos]
        records.append(record)

records.sort(reverse=True)
with open('plddt_ranking.txt', 'w') as f:
    for record in records:
        f.write(' '.join([str(x) for x in record]) + '\n')

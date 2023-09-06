import numpy as np
import os

records70 = []
records80 = []
records90 = []
infiles = sorted(os.listdir('clusters3'))
for infile in infiles:
    if infile[0] == '.':
        continue
    with open('clusters3/{}'.format(infile), 'r') as f:
        name = '{:10}'.format(infile.split('.')[0].split('_')[1])
        c1 = []
        c2 = []
        mean_local = []
        min_local = []
        lassos = set([])
        uniprot_ids = []
        lines = f.readlines()
        for line in lines[1:]:
            line = line.strip().split()
            uniprot_id, lasso = line[:2]
            plddts = line[-8:-2]
            bridge = line[-9]
            if bridge != '-':
                plddts = [float(x) for x in plddts]
                c1.append(plddts[0])
                c2.append(plddts[1])
                mean_local.append(plddts[2])
                min_local.append(plddts[3])
                lassos.add(lasso)
                uniprot_ids.append(uniprot_id)
        c1 = min(c1), max(c1)
        c2 = min(c2), max(c2)
        mean_local = round(min(mean_local),1), round(max(mean_local),1)
        min_local = min(min_local), max(min_local)
        lassos = ';'.join(sorted(lassos))
        uniprot_ids = ';'.join(uniprot_ids)
        for_scoring = [c1[1], c2[1], mean_local[0]]
        score = int(round(sum(for_scoring), 0))
        c1_str = '-'.join(map(str, c1))
        c2_str = '-'.join(map(str, c2))
        mean_local_str = '-'.join(map(str, mean_local))
        min_local_str = '-'.join(map(str, min_local))
        record = [score, name, c1_str, c2_str, mean_local_str, min_local_str, lassos, uniprot_ids]
        if all([x>70 for x in for_scoring]):
            if all([x>80 for x in for_scoring]):
                if all([x>90 for x in for_scoring]):
                    records90.append(record)
                else:
                    records80.append(record)
            else:
                records70.append(record)

records70.sort(reverse=True)
records80.sort(reverse=True)
records90.sort(reverse=True)
with open('plddt_ranking.txt', 'w') as f:
    f.write('score protein_id c1 c2 avg min lassos\n')
    for records in [records90, records80, records70]:
        f.write('=====================================\n')
        for record in records:
            f.write(' '.join([str(x) for x in record]) + '\n')

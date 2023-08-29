from collections import defaultdict

clusters = defaultdict(list)
cluster_size = {}
with open('protein2clusters_mapping.txt', 'r') as f:
    for line in f.readlines():
        record_id, uniprot_id, uniref_id, representant, group_size = line.strip().split()
        if uniref_id != '-':
            clusters[uniref_id].append(uniprot_id)
            if not uniref_id in cluster_size:
                cluster_size[uniref_id] = int(group_size)

results = []
for cluster, proteins in clusters.items():
    size = int(cluster_size[cluster])
    proteins_num = len(proteins)
    percent = round(proteins_num/size,2)
    result = [cluster, proteins_num, size, percent, ','.join(proteins)]
    results.append(result)
results.sort(key=lambda x:x[1])
with open('clusters.txt', 'w') as f:
    for result in results:
        f.write(' '.join([str(x) for x in result]) + '\n')

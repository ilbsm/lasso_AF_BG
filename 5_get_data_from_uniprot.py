import numpy as np
import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from collections import defaultdict
from tqdm import tqdm


with open('lassos_clusters_most_interesting.txt', 'r') as f:
    clusters = []
    for line in f.readlines():
        c = line.strip().split()[1]
        clusters += c.split(',')
    clusters = sorted(set(clusters))
fasta_representatives = []
for cluster in tqdm(clusters, total = len(clusters)):
    api_url = 'https://rest.uniprot.org/uniprotkb/search?query=uniref_cluster_50:{}&size=500'.format(cluster)
    response = requests.get(api_url)
    res = response.json()['results']
    cluster_data = {}
    for cluster_elem in res:
        uniprot_id = cluster_elem['primaryAccession']
        fasta = cluster_elem['sequence']['value'].replace('\n','')
        domains = []
        sdomains = []
        for ref in cluster_elem['uniProtKBCrossReferences']:
            if ref['database'] == 'Pfam':
                domains.append(ref['id'])
            elif ref['database'] == 'SUPFAM':
                sdomains.append(ref['id'])
        cluster_data[uniprot_id] = defaultdict(lambda: '-')
        cluster_data[uniprot_id]['fasta'] = fasta
        if domains:
            cluster_data[uniprot_id]['domains'] = ','.join(domains)
        if sdomains:
            cluster_data[uniprot_id]['sdomains'] = ','.join(sdomains)
    fasta_representatives.append('>{}\n{}'.format(cluster, cluster_data[cluster.split('_')[1]]['fasta']))
    with open('clusters/{}.txt'.format(cluster), 'r') as f:
        for line in f.readlines():
            uniprot_id, lasso, bridge = line.strip().split()
            cluster_data[uniprot_id]['lasso'] = lasso
            cluster_data[uniprot_id]['bridge'] = bridge
    fasta_records = []
    for uniprot_id in cluster_data.keys():
        uniprot_id_dict = cluster_data[uniprot_id]
        fasta_record = '>{}\n{}'.format(uniprot_id, uniprot_id_dict['fasta'])
        fasta_records.append(fasta_record)
        with open('fasta/{}.txt'.format(uniprot_id), 'w') as g:
            g.write(fasta_record)
    with open('fasta_clusters/{}.txt'.format(cluster), 'w') as g:
        g.write('\n'.join(fasta_records))
    with open('clusters2/{}.txt'.format(cluster), 'w') as g:
        for uniprot_id in cluster_data.keys():
            r = cluster_data[uniprot_id]
            record = [uniprot_id, r['lasso'], r['bridge'], r['domains'], r['sdomains']]
            g.write('{:10} {:28} {:7} {} {}\n'.format(*record))
with open('representatives_fasta.txt', 'w') as g:
    g.write('\n'.join(fasta_representatives))















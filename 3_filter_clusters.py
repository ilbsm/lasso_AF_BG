import requests
from tqdm import tqdm

data = {}
with open('lassos_surfaces_no_dupl.tsv', 'r') as f:
    for line in f.readlines():
        record_id, lasso_type, *_ = line.strip().split()
        uniprot_id = record_id.split('-')[1]
        bridge = record_id.split('_')[1]
        data[uniprot_id] = '{} {}'.format(lasso_type, bridge)

with open('clusters.txt', 'r') as f:
    for line in tqdm(f.readlines(), total=5504):
        cluster, num, family, perc, members =  line.strip().split()
        if int(num) > 1 and float(perc) >= .5:
            api_url = 'https://rest.uniprot.org/uniprotkb/search?query=uniref_cluster_50:{}'.format(cluster)
            response = requests.get(api_url)
            res = response.json()['results']
            uniprot_ids = []
            if res:
                with open('clusters/{}.txt'.format(cluster), 'w') as g:
                    for r in res:
                        uniprot_id = r['primaryAccession']
                        lasso_data = '- -'
                        if uniprot_id in data:
                            lasso_data = data[uniprot_id]
                        g.write('{} {}\n'.format(uniprot_id, lasso_data))

            

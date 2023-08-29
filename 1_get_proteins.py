import requests
from tqdm import tqdm
from collections import defaultdict

last_one = None
try:
    with open('protein2clusters_mapping.txt', 'r') as g:
        for line in g.readlines():
            last_one = line.split()[0]
except:
    pass

with open('lassos_surfaces_no_dupl.tsv', 'r') as f:
    for line in tqdm(f.readlines(), total=38602):
        record_id = line.split()[0]
        if last_one:
            if last_one == record_id:
                last_one = None
        else:
            uniprot_id = record_id.split('-')[1]
            api_url = 'https://rest.uniprot.org/uniref/search?query=uniprot_id:{}'.format(uniprot_id)
            response = requests.get(api_url)
            res = response.json()['results']
            uniref_id = '-'
            member_count = '0'
            representant = '-'
            for r in res:
                if r['entryType']=='UniRef50':
                    uniref_id = r['id']
                    member_count = str(r['memberCount'])
                    representant = r['representativeMember']['memberId']
            to_write = [record_id, uniprot_id, uniref_id, representant, member_count]
            with open('protein2clusters_mapping.txt', 'a+') as g:
                g.write('{}\n'.format(' '.join(to_write)))

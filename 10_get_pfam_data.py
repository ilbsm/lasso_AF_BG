import requests
import os
from tqdm import tqdm
from collections import defaultdict

domains_dict = defaultdict(list)
sdomains_dict = defaultdict(list)
infiles = sorted(os.listdir('clusters3'))
for infile in tqdm(infiles, total=len(infiles)):
    if infile[0] == '.': continue
    with open('clusters3/{}'.format(infile), 'r') as f:
        name = infile.split('.')[0].split('_')[1]
        lines = f.read().split('\n')
        lines[0] += ' Pfam   SUPFAM'
        for i, line in enumerate(lines[1:], start=1):
            if line:
                uniprot_id = line.split()[0]
                api_url = 'https://rest.uniprot.org/uniprotkb/search?query=accession:{}'.format(uniprot_id)
                response = requests.get(api_url)
                res = response.json()['results']
                domains = []
                sdomains = []
                for reference in res[0]['uniProtKBCrossReferences']:
                    if reference['database'] == 'Pfam':
                        domains.append(reference['id'])
                    elif reference['database'] == 'SUPFAM':
                        sdomains.append(reference['id'])
                domains.sort()
                sdomains.sort()
                if domains:
                    for domain in domains:
                        domains_dict[domain].append(name)
                    domains = ','.join(domains)
                else:
                    domains = '-'
                if sdomains:
                    for sdomain in sdomains:
                        sdomains_dict[sdomain].append(name)
                    sdomains = ','.join(sdomains)
                else:
                    sdomains = '-'
                lines[i] = '{} {} {}'.format(line, domains, sdomains)
    with open('clusters4/{}'.format(infile), 'w') as g:
        g.write('\n'.join(lines))

with open('domains.txt', 'w') as g:
    for k,v in domains_dict.items():
        g.write('{:10} {}\n'.format(k, ','.join(v)))
with open('sdomains.txt', 'w') as g:
    for k,v in sdomains_dict.items():
        g.write('{:10} {}\n'.format(k, ','.join(v)))




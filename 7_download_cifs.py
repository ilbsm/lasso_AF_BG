import os
import requests
import wget
from tqdm import tqdm
from collections import defaultdict
from urllib.error import HTTPError, URLError

clusters = sorted(os.listdir('clusters2'))
cifs_to_download = []
for cluster in clusters:
    if cluster[0] == '.':
        continue
    with open('clusters2/{}'.format(cluster), 'r') as f:
        for line in f.readlines():
            uniprot_id = line.split()[0]
            cifs_to_download.append(uniprot_id)

i = 0
with open('failed_downloads.txt', 'a+') as f:
    while i<len(cifs_to_download)-1:
        for i, uniprot_id in tqdm(enumerate(cifs_to_download), total=len(cifs_to_download)):
            filename = 'AF-{}-F1-model_v4.cif'.format(uniprot_id)
            if not os.path.isfile('cifs/{}'.format(filename)):
                url = 'https://alphafold.ebi.ac.uk/files/{}'.format(filename)
                try:
                    wget.download(url, out='cifs')
                except HTTPError:
                    f.write('{}\n'.format(uniprot_id))
                except URLError:
                    print('again')
                    break

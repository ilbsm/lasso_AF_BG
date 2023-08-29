import os
import requests
import wget
from tqdm import tqdm
from collections import defaultdict
from urllib.error import HTTPError, URLError

clusters = []
with open('lassos_clusters_most_interesting.txt', 'r') as f:
    for line in f.readlines():
        clusters += line.strip().split()[1].split(',')
clusters = sorted(set(clusters))

cifs_to_download = []
for cluster in clusters:
    with open('clusters/{}.txt'.format(cluster), 'r') as f:
        for line in f.readlines():
            cifs_to_download.append(line.strip().split()[0])

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

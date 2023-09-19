import os

TEMPLATES_PATH = "/home/bagren/dev/lasso_AF_BG/cifs_from_pdb"
MODEL_PATH = "/home/bagren/dev/lasso_AF_BG/cifs_from_pdb"
to_align = []

with open('proteins_from_lassoprot.txt', 'r') as f:
    for line in f.readlines():
        pdb, chain, *_ = line.split(';')
        to_align.append(pdb + '.cif')

for i, model in enumerate(to_align[:-1]):
    for template in to_align[i+1:]:
        os.system(f"pymol -c pymol_superimpose.pml -- {TEMPLATES_PATH.rstrip('/')}/{template} {MODEL_PATH.rstrip('/')}/{model}")


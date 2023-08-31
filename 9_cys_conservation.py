import os
from tqdm import tqdm

def find_indexes(protein):
    fasta = protein['fasta']
    cys_ndxes = protein['cys']
    list_ndxes = []
    resnum = 0
    for ndx, aa in enumerate(fasta):
        if aa != '-':
            resnum += 1
        if resnum in cys_ndxes:
            list_ndxes.append(ndx)
            if len(list_ndxes) == 2:
                return list_ndxes

infiles = sorted(os.listdir('clusters2'))

for infile in tqdm(infiles, total = len(infiles)):
    if infile[0] == '.': continue
    name = infile.split('.')[0]
    uniprot_ids = []
    uniprot_ids_AF = []
    proteins = {}
    with open('clusters2/{}'.format(infile), 'r') as f:
        for line in f.readlines():
            line = line.strip().split()
            uniprot_id, lasso, bridge, lasso_range, *plddts = line
            uniprot_ids.append(uniprot_id)
            if bridge != '-':
                uniprot_ids_AF.append(uniprot_id)
                bridge = tuple([int(x) for x in bridge.split('-')])
                proteins[uniprot_id] = {'cys':bridge, 'record':line, 'fasta':''}
            else:
                proteins[uniprot_id] = {'record':line, 'fasta':''}
    with open('align/{}.aln-clustal_num.clustal_num'.format(name), 'r') as f:
        for line in f.readlines():
            line = line.strip().split()
            if len(line) == 3:
                uniprot_id = line[0].split('_')[0]
                if uniprot_id in uniprot_ids:
                    proteins[uniprot_id]['fasta'] += line[1]
    c1 = set([])
    c2 = set([])
    for uniprot_id in uniprot_ids_AF:
        c1_ndx_list, c2_ndx_list = find_indexes(proteins[uniprot_id])
        c1.add(c1_ndx_list)
        c2.add(c2_ndx_list)
        proteins[uniprot_id]['cys_aligned'] = '{:d}-{:d}'.format(c1_ndx_list, c2_ndx_list)
    cysteines = sorted(c1 | c2)
    for uniprot_id in uniprot_ids:
        conserved = []
        print(uniprot_id)
        for cys in cysteines:
            print(cys, proteins[uniprot_id]['fasta'])
            if proteins[uniprot_id]['fasta'][cys] == 'C':
                conserved.append('+')
            else:
                conserved.append('-')
        proteins[uniprot_id]['conserved'] = ' '.join(conserved)
    with open('clusters3/{}'.format(infile), 'w') as g:
        cysteines_str = ' '.join([str(x) for x in cysteines])
        g.write('uniprot_id  lasso_type                   {} bridge_aligned bridge lass_range plddts(c1 c2 mean_local min_local mean_global min_global)\n'.format(cysteines_str))
        for uniprot_id in uniprot_ids:
            to_write1 = proteins[uniprot_id]['record'][:40]
            to_write2 = proteins[uniprot_id]['conserved']
            to_write3 = proteins[uniprot_id]['cys_aligned']
            to_write4 = proteins[uniprot_id]['record'][40:]
            g.write('{} {} {} {}\n'.format(to_write1, to_write2, to_write3, to_write4))
    









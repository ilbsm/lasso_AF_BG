import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm

def read_cif(cif_file):
    cifdict = MMCIF2Dict(cif_file)
    ndx_list = cifdict['_atom_site.label_seq_id']
    atom_list = cifdict['_atom_site.label_atom_id']
    plddt_list = cifdict['_atom_site.B_iso_or_equiv']
    fasta = cifdict['_entity_poly.pdbx_seq_one_letter_code']
    fasta2 = cifdict['_entity_poly.pdbx_seq_one_letter_code_can']
    if fasta != fasta2:
        raise
    fasta = fasta[0].replace('\n','')
    ndx2plddt = {}
    lists = ndx_list, atom_list, plddt_list
    for ndx, atom, plddt in zip(*lists):
        if atom == 'CA':
            ndx2plddt[int(ndx)] = float(plddt)
    return ndx2plddt, fasta

def get_plddt(ndx2plddt, bridge, lasso, edges):
    plddt_global = np.mean(list(ndx2plddt.values()))
    plddt_global_min = min(list(ndx2plddt.values()))
    if bridge != '-':
        cys1, cys2 = [int(x) for x in bridge.split('-')]
        beg, end = edges
        plddt_cys1 = ndx2plddt[cys1]
        plddt_cys2 = ndx2plddt[cys2]
        plddt_short = [v for k,v in ndx2plddt.items() if beg <= k <= end]
        plddt_local = np.mean(plddt_short)
        plddt_local_min = min(plddt_short)
    else:
        plddt_local = -1
        plddt_local_min = -1
        plddt_cys1 = -1
        plddt_cys2 = -1
    return plddt_cys1, plddt_cys2, plddt_local, plddt_local_min, plddt_global, plddt_global_min

def find_lasso_edges():
    uniid2edges = {}
    with open('lassos_surfaces_no_dupl.tsv', 'r') as f:
        for line in f.readlines():
            record_id, lasso, _, _, Ncross, Ccross, *_ = line.strip().split()
            uniprot_id, bridge = record_id.split('_')
            uniprot_id = uniprot_id.split('-')[1]
            ndxes = set([int(x) for x in bridge.split('-')])
            if Ncross != '""':
                ndxes |= set([int(x.strip('*+-')) for x in Ncross.split(',')])
            if Ccross != '""':
                ndxes |= set([int(x.strip('*+-')) for x in Ccross.split(',')])
            uniid2edges[uniprot_id] = min(ndxes), max(ndxes)
    return uniid2edges


def prepare():
    clusters = set([])
    with open('lassos_clusters_most_interesting.txt', 'r') as f:
        for line in f.readlines():
            cluster = line.strip().split()[1]
            clusters |= set(cluster.split(','))
    clusters = sorted(clusters)
    with open('failed_downloads.txt', 'r') as f:
        absent = set(f.read().split('\n'))
    uniid2edges = find_lasso_edges()
    return clusters, absent, uniid2edges

if __name__ == '__main__':
    clusters, absent, uniid2edges = prepare()
    for cluster in tqdm(clusters, total = len(clusters)):
        with open('clusters/{}.txt'.format(cluster), 'r') as f:
            with open('clusters2/{}.txt'.format(cluster), 'w') as g:
                with open('fasta/{}.txt'.format(cluster), 'w') as gg:
                    to_write = []
                    to_write_no_struct = []
                    to_write_else = []
                    for line in f.readlines():
                        uniprot_id, lasso, bridge = line.strip().split()
                        if uniprot_id in absent:
                            lasso = 'no_struct'
                            plddts_str = '   '.join(['-1']*6)
                            edges = bridge
                        else:
                            ndx2plddt, fasta = read_cif('cifs/AF-{}-F1-model_v4.cif'.format(uniprot_id))
                            if uniprot_id in uniid2edges:
                                edges = uniid2edges[uniprot_id]
                            else:
                                edges = bridge
                            plddts = get_plddt(ndx2plddt, bridge, lasso, edges)
                            plddts = [str(round(x,1)) for x in plddts]
                            plddts_str = '{:4} {:4} {:4} {:4} {:4} {:4}'.format(*plddts)
                            gg.write('>{}\n'.format(uniprot_id))
                            gg.write(fasta + '\n')
                        edges = '-'.join([str(x) for x in edges])
                        record = [uniprot_id, lasso, bridge, edges, plddts_str]
                        if lasso[0] == 'L':
                            to_write.append(record)
                        elif lasso[0] == 'n':
                            to_write_no_struct.append(record)
                        else:
                            to_write_else.append(record)
                    to_write += to_write_else
                    to_write += to_write_no_struct
                    for record in to_write:
                        g.write('{:10} {:28} {:7} {:7} {}\n'.format(*record))















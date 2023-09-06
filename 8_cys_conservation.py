import os
import statistics
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm
from collections import defaultdict

def read_cif(cif_file):
    if os.path.isfile(cif_file):
        cifdict = MMCIF2Dict(cif_file)
        ndx_list = cifdict['_atom_site.label_seq_id']
        atom_list = cifdict['_atom_site.label_atom_id']
        plddt_list = cifdict['_atom_site.B_iso_or_equiv']
        ndx2plddt = {}
        lists = ndx_list, atom_list, plddt_list
        for ndx, atom, plddt in zip(*lists):
            if atom == 'CA':
                ndx2plddt[int(ndx)] = float(plddt)
    else:
        ndx2plddt = defaultdict(lambda: '-')
    return ndx2plddt

def get_plddt(ndx2plddt, bridge, edges):
    if ndx2plddt.values():
        plddt_global = statistics.mean(list(ndx2plddt.values()))
        plddt_global_min = min(list(ndx2plddt.values()))
        if bridge != '-':
            #cys1, cys2 = [int(x) for x in bridge.split('-')]
            cys1, cys2 = bridge
            beg, end = edges
            plddt_cys1 = ndx2plddt[cys1]
            plddt_cys2 = ndx2plddt[cys2]
            plddt_short = [v for k,v in ndx2plddt.items() if beg <= k <= end]
            plddt_local = statistics.mean(plddt_short)
            plddt_local_min = min(plddt_short)
        else:
            plddt_local = '-'
            plddt_local_min = '-'
            plddt_cys1 = '-'
            plddt_cys2 = '-'
        return plddt_cys1, plddt_cys2, plddt_local, plddt_local_min, plddt_global, plddt_global_min
    else:
        return '-','-','-','-','-','-'

def find_lasso_edges():
    uniid2edges = {}
    uniid2bridges = {}
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
            uniid2bridges[uniprot_id] = tuple([int(x) for x in bridge.split('-')])
    return uniid2edges, uniid2bridges

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

def read_fasta(name, uniprot_ids):
    fastas = defaultdict(str)
    with open('align/{}.aln-clustal_num.clustal_num'.format(name), 'r') as f:
        for line in f.readlines():
            line = line.strip().split()
            if len(line) == 3:
                uniprot_id = line[0].split('_')[0]
                if uniprot_id in uniprot_ids:
                    fastas[uniprot_id] += line[1]
    return fastas

def load_cluster(infile):
    proteins = {}
    uniprot_ids = []
    with open('clusters2/{}'.format(infile), 'r') as f:
        for line in f.readlines():
            line = line.strip().split()
            uniprot_id, lasso, bridge, domains, sdomains = line
            uniprot_ids.append(uniprot_id)
            proteins[uniprot_id] = defaultdict(lambda: '-')
#            if bridge != '-':
#                    uniprot_ids_AF.append(uniprot_id)
#                bridge = tuple([int(x) for x in bridge.split('-')])
#                proteins[uniprot_id]['cys'] = bridge
            proteins[uniprot_id]['record'] = line
    return proteins, uniprot_ids

def format_bridge(bridge):
    if bridge == '-':
        return '   -   '
    elif type(bridge) == str:
        return bridge
    else:
        return '{:3d}-{:3d}'.format(*bridge)

if __name__ == '__main__':
    uniid2edges, uniid2bridges = find_lasso_edges()
    infiles = sorted(os.listdir('clusters2'))
    for infile in tqdm(infiles, total = len(infiles)):
        if infile[0] == '.': continue
        name = infile.split('.')[0]
        proteins, uniprot_ids = load_cluster(infile)
        fastas = read_fasta(name, uniprot_ids)
        for uniprot_id, fasta in fastas.items():
            proteins[uniprot_id]['fasta'] = fasta
        for uniprot_id in uniprot_ids:
            bridge = proteins[uniprot_id]['cys']
            cif_file = 'cifs/AF-{}-F1-model_v4.cif'.format(uniprot_id)
            ndx2plddt = read_cif(cif_file)
            if uniprot_id in uniid2edges: # it means that lasso have been found
                proteins[uniprot_id]['range'] = uniid2edges[uniprot_id]
                proteins[uniprot_id]['cys'] = uniid2bridges[uniprot_id]
            edge = proteins[uniprot_id]['range']
            bridge = proteins[uniprot_id]['cys']
            plddts = get_plddt(ndx2plddt, bridge, edge)
            plddts = [str(round(x,1)) if type(x)==float else x for x in plddts]
            plddts_str = '{:4} {:4} {:4} {:4} {:4} {:4}'.format(*plddts)
            proteins[uniprot_id]['plddts'] = plddts_str
        c1 = set([])
        c2 = set([])
        for uniprot_id in uniprot_ids:
            if proteins[uniprot_id]['cys'] != '-' and proteins[uniprot_id]['fasta'] != '-':
                c1_ndx_list, c2_ndx_list = find_indexes(proteins[uniprot_id])
                c1.add(c1_ndx_list)
                c2.add(c2_ndx_list)
                proteins[uniprot_id]['cys_aligned'] = '{:3d}-{:3d}'.format(c1_ndx_list, c2_ndx_list)
        cysteines = sorted(c1 | c2)
        for uniprot_id in uniprot_ids:
            conserved = []
            for cys in cysteines:
                if proteins[uniprot_id]['fasta'][cys] == 'C':
                    conserved.append('+')
                else:
                    conserved.append('-')
            proteins[uniprot_id]['conserved'] = '   '.join(conserved) + ' '
        # sorting output
        a = []
        b = []
        c = []
        d = []
        for uniprot_id in uniprot_ids:
            if proteins[uniprot_id]['range'] != '-':
                a.append(uniprot_id)
            elif proteins[uniprot_id]['record'][1] != '-':
                b.append(uniprot_id)
            elif proteins[uniprot_id]['plddts'][-4] != '-':
                c.append(uniprot_id)
            else:
                d.append(uniprot_id)
        uniprot_ids = a+b+c+d
        with open('clusters3/{}'.format(infile), 'w') as g:
            cysteines_str = ' '.join([str(x) for x in cysteines])
            g.write('uniprot_id lasso_type                  {} bridge_aligned bridge lasso_range plddts(c1 c2 mean_local min_local mean_global min_global) Pfam SUPFAM\n'.format(cysteines_str))
            for uniprot_id in uniprot_ids:
                prot = proteins[uniprot_id]
                to_write1 = '{:10} {:28}'.format(*prot['record'][:2])
                to_write2 = prot['conserved']
                to_write3 = '{} {} {}'.format(format_bridge(prot['cys_aligned']),
                                              format_bridge(prot['cys']),
                                              format_bridge(prot['range']))
                to_write4 = proteins[uniprot_id]['plddts']
                to_write5 = '{} {}'.format(*proteins[uniprot_id]['record'][3:])
                to_write = [to_write1, to_write2, to_write3, to_write4, to_write5]
                g.write('{} {} {} {} {}\n'.format(*to_write))
        







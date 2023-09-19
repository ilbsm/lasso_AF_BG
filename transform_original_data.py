import os
import statistics
from copy import deepcopy
from collections import defaultdict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def read_cif(uniprot_id):
    cif_file = 'cifs/AF-{}-F1-model_v4.cif'.format(uniprot_id)
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
    cys1, cys2 = [int(x) for x in bridge.split('-')]
    beg, end = [int(x) for x in edges.split('-')]
    plddt_cys1 = ndx2plddt[cys1]
    plddt_cys2 = ndx2plddt[cys2]
    plddt_short = [v for k,v in ndx2plddt.items() if beg <= k <= end]
    plddt_local = statistics.mean(plddt_short)
    plddt_local_min = min(plddt_short)
    return plddt_cys1, plddt_cys2, plddt_local, plddt_local_min#, plddt_global, plddt_global_min

if __name__ == '__main__':
    with open('lassos_surfaces_no_dupl.tsv', 'r') as f:
        records = []
        for line in f.readlines():
            record_id, lasso, _,_, cross_N, cross_C, *_ = line.strip().split()
            uniprot_id, bridge = record_id.split('_')
            bridge = bridge.split('-')
            uniprot_id = uniprot_id.split('-')[1]
            resids = deepcopy(bridge)
            if cross_N != '""':
                resids += cross_N.replace('+','').replace('-','').replace('*','').split(',')
            if cross_C != '""':
                resids += cross_C.replace('+','').replace('-','').replace('*','').split(',')
            resids = [int(x) for x in resids]
            lasso_range = min(resids), max(resids)
            #bridge_fmt = '{:>3}-{:<4}'.format(*bridge)
            #lasso_range_fmt = '{:>3}-{:<4}'.format(*lasso_range)
            #ndx2plddt = read_cif(uniprot_id)
            #plddts = get_plddt(ndx2plddt, bridge_fmt, lasso_range_fmt)
            #plddts = [str(round(x,1)) if type(x)==float else x for x in plddts] 
            #plddts_str = '{:4} {:4} {:4} {:4}'.format(*plddts)  
            #record = '{:10} {} {} {} {}'.format(uniprot_id, bridge_fmt, lasso_range_fmt, plddts_str, lasso)
            record = '{:10} {:>3}-{:<4} {:>3}-{:<4} {}'.format(uniprot_id, *bridge, *lasso_range, lasso)
            records.append(record)
    with open('prots_gln_gt2_more.txt', 'w') as g:
        g.write('#uniprot_id bridge lasso_range plddt(c1 c2 avg min) lasso_type\n')
        g.write('\n'.join(records))
     

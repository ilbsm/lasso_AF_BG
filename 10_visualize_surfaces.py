import os

with open('plddt_ranking.txt', 'r') as f:
    for line_str in f.readlines():
        line = line_str.strip().split()
        if line[0].isnumeric():
            cluster = line[1]
            with open('clusters3/UniRef50_{}.txt'.format(cluster), 'r') as ff:
                record = ff.read().split('\n')[1].split()
                print(record)
                uniprot_id = record[0]
                bridge = record[-10]
            print('Plot lasso? (write anything to abort)')
            x = input()
            if x: break
            pdb_file = 'surfaces/{}/AF-{}-F1.pdb'.format(uniprot_id[-2:], uniprot_id)
            tcl_file = 'surfaces/{}/AF-{}-F1_{}.tcl'.format(uniprot_id[-2:], uniprot_id, bridge)
            print(line_str)
            os.system('vmd {} -e {}'.format(pdb_file, tcl_file))   

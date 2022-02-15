'''
Generate counts for sim

Author: Jorge Mestre Tomas
Date: 15/02/2022
Last update: 15/02/2022
'''

import argparse
import random


def main():
    parser = argparse.ArgumentParser(prog='sqanti3_sim.py', description="SQANTI-SIM: a simulator of controlled novelty and degradation of transcripts sequence by long-reads")
    parser.add_argument('--cat', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('--deleted', default = False,  help = '\t\tFile with deleted trans', required=True)
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPath for output files')


    args = parser.parse_args()

    n_trans = 10000
    counts = str(5)
    tpm = str(5)

    deleted = set()
    with open(args.deleted, 'r') as del_in:
        skip = del_in.readline()
        for line in del_in:
            line = line.split()
            deleted.add(line[0])
    del_in.close()


    tot_trans = deleted
    with open(args.cat, 'r') as cat_in:
        skip = cat_in.readline()
        for line in cat_in:
            trans = line.split()[0]
            if trans not in deleted:
                tot_trans.add(trans)
            if len(tot_trans) == n_trans:
                break

    cat_in.close()
    
    f_out = open(args.output, 'w')
    f_out.write('target_id\test_counts\ttpm\n')
    for trans in tot_trans:
        f_out.write(trans + '\t' + counts + '\t' + tpm + '\n')
    f_out.close()

if __name__ == '__main__':
    main()
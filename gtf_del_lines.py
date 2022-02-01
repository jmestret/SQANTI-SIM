#!/usr/bin/env python3
'''
gtf_del_lines.py
Modify original GTF deleting transcripts to simulate reads

Author: Jorge Mestre Tomas
Date: 19/01/2020
Last update: 01/02/2021 by Jorge Mestre
'''

__author__ = 'jormart2@alumni.uv.es'
__version__ = '0.0'

import random
from time import time


def target_trans(f_name, counts):
    # TODO: use sets intead of lists

    trans4SC = {
    }

    target_trans = []
    target_genes = []
    ref_trans = []
    ref_genes = []

    with open(f_name, 'r') as gtf:
        header = gtf.readline()
        for line in gtf:
            line_split = line.split()
            SC = line_split[2]

            if SC in list(trans4SC.keys()):
                trans4SC[SC].append(tuple(line_split))
            else:
                trans4SC[SC] = [tuple(line_split)]
    
    gtf.close()

    for SC in list(counts.keys()):
        if counts[SC] > 0:
            SCtrans = trans4SC[SC]
            random.Random(1234).shuffle(SCtrans)
            for trans in SCtrans:
                trans_id = trans[0]
                gene_id = trans[1]
                SC = trans[2]

                if SC in ['FSM', 'ISM'] and counts[SC] > 0 and trans_id not in ref_trans and gene_id not in ref_genes:
                    ref_t = trans[4]
                    if ref_t not in target_trans:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_trans.append(ref_t)
                        counts[SC] -= 1

                elif SC in ['NIC', 'NNC'] and counts[SC] > 0 and gene_id not in ref_genes:
                    ref_g = trans[3]
                    if ref_g not in target_genes:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_genes.append(ref_g)
                        counts[SC] -= 1
                
                if counts[SC] <= 0:
                    break

    return target_trans, ref_genes, ref_trans

def getGeneID(line):
    line_split = line.split()
    gene_id = line_split[line_split.index('gene_id') + 1]
    gene_id = gene_id.replace(';', '').replace('"', '')
    
    return gene_id

def getTransID(line):
    try:
        line_split = line.split()
        trans_id = line_split[line_split.index('transcript_id') + 1]
        trans_id = trans_id.replace(';', "").replace('"', "")
    except:
        trans_id = None

    return trans_id

def at_least_one_trans(gene_info):
    for i in gene_info:
        line_split = i.split()
        if line_split[2] == 'exon':
            return True
    return False

def modifyGTF(f_name_in, f_name_out, target_trans, ref_trans, ref_genes):
    f_out = open(f_name_out, 'w')

    gene_info = []
    prev_gene = str()

    with open(f_name_in, 'r') as gtf_in:
        for line in gtf_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                gene_id = getGeneID(line)
                if not prev_gene:
                    prev_gene = gene_id
                
                if gene_id == prev_gene:
                    gene_info.append(line)

                elif prev_gene in ref_genes:
                    for i in gene_info:
                        f_out.write(i)
                    prev_gene = gene_id
                    gene_info =  [line]

                else:
                    tmp = []
                    for i in gene_info:
                        trans_id = getTransID(i)
                        if trans_id not in target_trans:
                            tmp.append(i)
                    
                    if at_least_one_trans(tmp):
                        for i in tmp:
                            f_out.write(i)
                            
                    prev_gene = gene_id
                    gene_info =  [line]
    gtf_in.close()

    if prev_gene in ref_genes:
        for i in gene_info:
            f_out.write(i)
            prev_gene = gene_id
            gene_info =  [line]

    else:
        tmp = []
        for i in gene_info:
            trans_id = getTransID(i)
            if trans_id not in target_trans:
                tmp.append(i)
                    
        if at_least_one_trans(tmp):
            for i in tmp:
                f_out.write(i)
    
    f_out.close()

    return


def main():
    gtf2SC_name = '/home/jorge/Desktop/prueba_gtf2SC.txt'
    gtf_name = '/home/jorge/Desktop/simulacion/getSC/chr3.gencode.v38.annotation.gtf'
    gtf_out = '/home/jorge/Desktop/simulacion/modifgtf/chr3.gencode.v38.annotation.ISM.del.gtf'
    counts = {
        'FSM' : 0, 'ISM' : 0, 'NIC' : 100, 'NNC' : 0
    }
    target, ref_genes, ref_trans = target_trans(gtf2SC_name, counts)

    modifyGTF(gtf_name, gtf_out, target, ref_genes, ref_trans)

if __name__ == '__main__':
    t_ini = time()
    main()
    t_fin = time()
    print('[Execution time %s seconds]' %(t_fin-t_ini))
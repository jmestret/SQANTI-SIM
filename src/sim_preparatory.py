#!/usr/bin/env python3
'''
sim_preparatory.py
Generate counts for sim

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 03/03/2022
'''

import os
import sys
import subprocess
import random
import numpy
import pysam
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt


def target_trans(f_name: str, f_name_out: str, counts: dict, seed: int)-> tuple:
    '''
    Choose those transcripts that will be deleted from the original GTF
    to generate the modified file to use as the reference annotation

    Args:
        f_name (str) name of the file with the GTF classification
        counts (dict) dictinary with the number of transcripts of each SC to be
                      deleted
    '''

    trans_by_SC = defaultdict(lambda: [])
    trans_by_gene = defaultdict(lambda: [])

    target_trans = set()
    target_genes = set()
    ref_trans = set()
    ref_genes = set()

    # Build a list for each SC with all transcripts that were classified there
    with open(f_name, 'r') as cat:
        col_names = cat.readline()
        for line in cat:
            line_split = line.split()
            gene = line_split[1]
            SC = line_split[2]

            trans_by_SC[SC].append(tuple(line_split))
            trans_by_gene[gene].append(tuple(line_split))
    
    cat.close()

    # Select randomly the transcripts of each SC that are going to be deleted
    # It's important to make sure you don't delete its reference trans or gene
    for SC in counts:
        if counts[SC] > 0:
            SCtrans = trans_by_SC[SC]
            random.Random(seed).shuffle(SCtrans)
            for trans in SCtrans:
                trans_id = trans[0]
                gene_id = trans[1]
                SC = trans[2]
                ref_g = trans[3]
                ref_t = trans[4]

                if SC in ['full-splice_match', 'incomplete-splice_match'] and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and ref_t not in target_trans:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        ref_trans.add(ref_t)
                        counts[SC] -= 1

                elif SC in ['novel_in_catalog', 'novel_not_in_catalog', 'genic_intron'] and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and gene_id not in target_genes and ref_g not in target_genes:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        ref_genes.add(ref_g)
                        counts[SC] -= 1

                elif SC in ['fusion', 'antisense', 'genic'] and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and gene_id not in target_genes:
                        ref_g = trans[3].split('_')
                        for i in ref_g:
                            if i in target_genes:
                                break
                        else:
                            target_trans.add(trans_id)
                            target_genes.add(gene_id)
                            for i in ref_g:
                                ref_genes.add(i)
                            counts[SC] -= 1
                
                elif SC == 'intergenic' and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and gene_id not in target_genes:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        counts[SC] -= 1
                
                if counts[SC] <= 0:
                    break

    final_target = target_trans
    for gene in trans_by_gene:
        for trans in trans_by_gene[gene]:
            if trans[0] in target_trans:
                trans_by_gene[gene].remove(trans)
                if len(trans_by_gene[gene]) == 0:
                    final_target.add(gene)

    f_out = open(f_name_out, 'w')
    with open(f_name, 'r') as cat:
        col_names = cat.readline()
        col_names = col_names.split()
        col_names.append('sim_type')
        f_out.write('\t'.join(col_names) + '\n')
        for line in cat:
            line = line.split()
            trans_id = line[0]
            if trans_id in target_trans:
                line.append('novel')
            else:
                line.append('known')
            f_out.write('\t'.join(line) + '\n')
    
    cat.close()
    f_out.close()
    
    return final_target
    #return target_trans, ref_genes, ref_trans


def getGeneID(line: str)-> str:
    '''
    Returns the gene_id of a GTF line

    Args:
        line (str) line readed from GTF file

    Returns:
        gene_id (str) gene_id from that feature
    '''

    line_split = line.split()
    gene_id = line_split[line_split.index('gene_id') + 1]
    gene_id = gene_id.replace(';', '').replace('"', '')
    
    return gene_id


def getTransID(line: str)-> str:
    '''
    Returns the transcript_id of a GTF line

    Args:
        line (str) line readed from GTF file

    Returns:
        trans_id (str) transcript_id from that feature
    '''

    try:
        line_split = line.split()
        trans_id = line_split[line_split.index('transcript_id') + 1]
        trans_id = trans_id.replace(';', "").replace('"', "")
    except:
        trans_id = None

    return trans_id


def modifyGTF(f_name_in: str, f_name_out: str, target: list):
    '''
    Modify the original GTF deleting target transcripts to simulate specific
    SQANTI3 structural categorires

    Args:
        f_name_in (str) file name of the reference annotation GTF
        f_name_out (str) file name of the modified GTF generated
        target_trans (list) list of transcripts that will be deleted
        ref_genes (list) list of genes that can't be deleted
    '''
    
    f_out = open(f_name_out, 'w')

    with open(f_name_in, 'r') as gtf_in:
        for line in gtf_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                gene_id = getGeneID(line)
                trans_id = getTransID(line)
                if gene_id in target or trans_id in target:
                    pass
                else:
                    f_out.write(line)
    gtf_in.close()
    f_out.close()

    return


def simulate_gtf(args):
    print('***Writting modified GTF\n')
    counts = defaultdict(lambda: 0, {
        'full-splice_match': 0,
        'incomplete-splice_match': args.ISM,
        'novel_in_catalog': args.NIC,
        'novel_not_in_catalog':args.NNC,
        'fusion' : args.Fusion,
        'antisense': args.Antisense,
        'genic_intron': args.GI,
        'genic' :args.GG,
        'intergenic':args.Intergenic
    })

    gtf_modif = os.path.join(args.dir, (args.output + '_modified.gtf'))
    tmp = os.path.join(os.path.dirname(os.path.abspath(args.trans_index)),'tmp_preparatory.tsv')

    target = target_trans(args.trans_index, tmp, counts, args.seed)
    modifyGTF(args.gtf, gtf_modif, target)
    os.remove(args.trans_index)
    os.rename(tmp, args.trans_index)

    return counts


def summary_table_del(counts_ini: dict, counts_end: dict):
    for sc in counts_end:
        counts_ini[sc] -= counts_end[sc]
    
    print('\033[94m_' * 79 + '\033[0m')
    print('\033[92mS Q A N T I - S I M\033[0m \U0001F4CA')
    print()
    print('GTF modification summary Table \U0001F50E')
    print('\033[94m_' * 79 + '\033[0m')
    for k, v in counts_ini.items():
        print('\033[92m|\033[0m ' + k + ': ' + str(v))


def create_expr_file_fixed_count(f_idx: str, n_trans: int, read_count: int, expr_out: str, index_out:str):

    novel_trans = []
    known_trans = []

    with open(f_idx, 'r') as f_in:
        skip = f_in.readline()
        for line in f_in:
            line = line.split()
            sim_type = line[11]
            if sim_type == 'novel':
                novel_trans.append(line[0])
            else:
                known_trans.append(line[0])
    f_in.close()

    random.shuffle(known_trans)
    known_trans = known_trans[:(n_trans-len(novel_trans))]

    tot_trans = len(novel_trans) + len(known_trans)
    if tot_trans != n_trans:
        print('Warning: A higher number than annotated transcripts was requested to simulates, only %s transcript will be simulated' %(tot_trans))

    coverage = read_count//n_trans
    tpm = (1000000.0 * coverage) / (coverage * n_trans) # Not taking into account transcript length

    if f_idx == index_out:
        tmp = os.path.join(os.path.dirname(os.path.abspath(f_idx)),'tmp_preparatory.tsv')
    else:
        tmp = index_out
    f_out = open(tmp, 'w')
    with open(f_idx, 'r') as idx:
        col_names = idx.readline()
        col_names = col_names.split()
        col_names.extend(['requested_counts', 'requested_tpm'])
        f_out.write('\t'.join(col_names) + '\n')
        for line in idx:
            line = line.split()
            trans_id = line[0]
            if trans_id in novel_trans or trans_id in known_trans:
                line.extend([str(coverage), str(tpm)])
            else:
                line.extend(['0', '0'])
            f_out.write('\t'.join(line) + '\n')
    idx.close()
    f_out.close()

    if f_idx == index_out:
        os.remove(f_idx)
        os.rename(tmp, f_idx)

    sns.set(style="whitegrid")
    sns.histplot([coverage]*len(known_trans), color="skyblue", label="Known", kde=True)
    sns.histplot([coverage]*len(novel_trans), color="red", label="Novel", kde=True)
    plt.legend()
    plt.savefig(''.join(expr_out.split('.')[:-1]) + '.png')


def create_expr_file_nbinom(f_idx: str, n_trans, nbn_known, nbp_known, nbn_novel, nbp_novel, output: str, index_out: str):
    novel_trans = []
    known_trans = []

    with open(f_idx, 'r') as f_in:
        skip = f_in.readline()
        for line in f_in:
            line = line.split()
            sim_type = line[11]
            if sim_type == 'novel':
                novel_trans.append(line[0])
            else:
                known_trans.append(line[0])
    f_in.close()

    random.shuffle(known_trans)
    known_trans = known_trans[:(n_trans-len(novel_trans))]

    nb_known = numpy.random.negative_binomial(nbn_known,nbp_known,len(known_trans)).tolist()
    nb_known = [1 if n == 0 else n for n in nb_known] # minimum one count per transcript
    nb_novel = numpy.random.negative_binomial(nbn_novel,nbp_novel,len(novel_trans)).tolist()
    nb_novel = [1 if n == 0 else n for n in nb_novel] # minimum one count per transcript
    n_reads =sum(nb_known) + sum(nb_novel)

    if f_idx == index_out:
        tmp = os.path.join(os.path.dirname(os.path.abspath(f_idx)),'tmp_preparatory.tsv')
    else:
        tmp = index_out
    f_out = open(tmp, 'w')
    i_known = 0
    i_novel = 0
    with open(f_idx, 'r') as idx:
        col_names = idx.readline()
        col_names = col_names.split()
        col_names.extend(['requested_counts', 'requested_tpm'])
        f_out.write('\t'.join(col_names) + '\n')
        for line in idx:
            line = line.split()
            trans_id = line[0]
            if trans_id in novel_trans:
                coverage = nb_novel[i_novel]
                i_novel += 1
            elif trans_id in known_trans:
                coverage = nb_known[i_known]
                i_known += 1
            else:
                coverage = 0

            tpm = round(((1000000.0 * coverage) / n_reads), 2)
            line.extend([str(coverage), str(tpm)])
            f_out.write('\t'.join(line) + '\n')
    idx.close()
    f_out.close()

    if f_idx == index_out:
        os.remove(f_idx)
        os.rename(tmp, f_idx)

    sns.set(style="whitegrid")
    sns.histplot(nb_known, color="skyblue", label="Known", kde=True)
    sns.histplot(nb_novel, color="red", label="Novel", kde=True)
    plt.legend()
    plt.savefig(''.join(output.split('.')[:-1]) + '.png')


def create_expr_file_sample(f_idx: str, ref_trans,reads, output: str, index_out:str, tech):
    sam_file = ''.join(output.split('.')[:-1]) + '_' + tech + '.sam'

    if tech == 'pb':
        cmd = ['minimap2', ref_trans, reads, '-x', 'map-pb',
               '-a', '--secondary=no', '-o', sam_file]
    elif tech == 'ont':
        cmd=['minimap2', ref_trans, reads, '-x', 'map-ont',
             '-a', '--secondary=no', '-o', sam_file]

    if subprocess.check_call(cmd, shell=True)!=0:
        print('ERROR running minimap2: {0}'.format(cmd), file=sys.stderr)
        sys.exit(1)
    
    trans_counts = defaultdict(lambda: 0)

    with pysam.AlignmentFile(sam_file, 'r') as sam_file_in:
        for align in sam_file_in:
            trans_id = align.reference_name

            if align.reference_id == -1 or align.is_supplementary or align.is_secondary:
                continue
            trans_counts[trans_id] += 1
    os.remove(sam_file)

    expr_distr = list(trans_counts.values())
    expr_distr.sort()
    n_trans = len(expr_distr)

    novel_trans = []
    known_trans = []

    with open(f_idx, 'r') as f_in:
        skip = f_in.readline()
        for line in f_in:
            line = line.split()
            sim_type = line[11]
            if sim_type == 'novel':
                novel_trans.append(line[0])
            else:
                known_trans.append(line[0])
    f_in.close()

    random.shuffle(known_trans)
    known_trans = known_trans[:(n_trans-len(novel_trans))]

    if (len(novel_trans) + len(known_expr)) < n_trans:
        n_trans = len(novel_trans) + len(known_expr)
        expr_distr = expr_distr[-n_trans,]

    lim_novel = (len(novel_trans)//3)*2
    lim_known = (len(known_trans)//3)*2

    novel_expr = expr_distr[:lim_novel]
    known_expr =  expr_distr[-lim_known:]

    expr_distr = expr_distr[lim_novel:lim_known]
    random.shuffle(expr_distr)

    novel_expr.extend(expr_distr[:len(novel_trans)-lim_novel])
    known_expr.extend(expr_distr[-(len(known_trans)-lim_known):])
    n_reads =sum(novel_expr) + sum(known_expr)


    if f_idx == index_out:
        tmp = os.path.join(os.path.dirname(os.path.abspath(f_idx)),'tmp_preparatory.tsv')
    else:
        tmp = index_out
    f_out = open(tmp, 'w')
    i_known = 0
    i_novel = 0
    with open(f_idx, 'r') as idx:
        col_names = idx.readline()
        col_names = col_names.split()
        col_names.extend(['requested_counts', 'requested_tpm'])
        f_out.write('\t'.join(col_names) + '\n')
        for line in idx:
            line = line.split()
            trans_id = line[0]
            if trans_id in novel_trans:
                coverage = novel_expr[i_novel]
                i_novel += 1
            elif trans_id in known_trans:
                coverage = known_expr[i_known]
                i_known += 1
            else:
                coverage = 0

            tpm = round(((1000000.0 * coverage) / n_reads), 2)
            line.extend([str(coverage), str(tpm)])
            f_out.write('\t'.join(line) + '\n')
    idx.close()
    f_out.close()

    if f_idx == index_out:
        os.remove(f_idx)
        os.rename(tmp, f_idx)

    sns.set(style="whitegrid")
    fig = sns.kdeplot(novel_expr, shade=True, color="r")
    fig = sns.kdeplot(known_expr, shade=True, color="b")
    sns.histplot(known_expr, color="skyblue", label="Known", kde=True)
    sns.histplot(novel_expr, color="red", label="Novel", kde=True)
    plt.legend()
    plt.savefig(''.join(output.split('.')[:-1]) + '.png')
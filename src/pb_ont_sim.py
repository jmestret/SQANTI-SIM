#!/usr/bin/env python3
'''
pb_ont_sim.py
Generate counts for sim

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 20/02/2022
'''

import logging
import os
import sys
import subprocess
import random
import numpy
import pysam
from collections import defaultdict


#------------------------------------
# MODIFY GTF FUNCTIONS

def target_trans(f_name: str, f_name_out: str, counts: dict)-> tuple:
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
    with open(f_name, 'r') as gtf:
        header = gtf.readline()
        for line in gtf:
            line_split = line.split()
            gene = line_split[1]
            SC = line_split[2]

            trans_by_SC[SC].append(tuple(line_split))
            trans_by_gene[gene].append(tuple(line_split))
    
    gtf.close()

    f_out = open(f_name_out, 'w')
    f_out.write('TransID\tGeneID\tSC\tRefGene\tRefTrans\tTSS\tTTS\tDonors\tAcceptors\n')

    # Select randomly the transcripts of each SC that are going to be deleted
    # It's important to make sure you don't delete its reference trans or gene
    for SC in counts:
        if counts[SC] > 0:
            SCtrans = trans_by_SC[SC]
            random.Random(1234).shuffle(SCtrans)
            for trans in SCtrans:
                trans_id = trans[0]
                gene_id = trans[1]
                SC = trans[2]
                ref_g = trans[3]
                ref_t = trans[4]

                if SC in ['full-splice_match', 'incomplete-splice_match', 'genic_intron'] and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and gene_id not in target_genes and ref_t not in target_trans:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        ref_trans.add(ref_t)
                        counts[SC] -= 1
                        f_out.write('\t'.join(trans))
                        f_out.write('\n')

                elif SC in ['novel_in_catalog', 'novel_not_in_catalog'] and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and gene_id not in target_genes and ref_g not in target_genes:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        ref_genes.add(ref_g)
                        counts[SC] -= 1
                        f_out.write('\t'.join(trans))
                        f_out.write('\n')

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
                            f_out.write('\t'.join(trans))
                            f_out.write('\n')
                
                elif SC == 'intergenic' and counts[SC] > 0:
                    if trans_id not in ref_trans and gene_id not in ref_genes and gene_id not in target_genes:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        counts[SC] -= 1
                        f_out.write('\t'.join(trans))
                        f_out.write('\n')
                
                if counts[SC] <= 0:
                    break

    final_target = target_trans
    for gene in trans_by_gene:
        for trans in trans_by_gene[gene]:
            if trans[0] in target_trans:
                trans_by_gene[gene].remove(trans)
                if len(trans_by_gene[gene]) == 0:
                    final_target.add(gene)

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
    del_trans = os.path.join(args.dir, (args.output + '_deleted.txt'))

    target = target_trans(args.cat, del_trans, counts)
    modifyGTF(args.gtf, gtf_modif, target)

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


def create_expr_file_fixed_count(f_cat: str, f_del: str, n_trans: int, read_count: int, output: str):

    deleted = set()
    with open(f_del, 'r') as del_in:
        skip = del_in.readline()
        for line in del_in:
            line = line.split()
            deleted.add(line[0])
    del_in.close()


    tot_trans = deleted
    with open(f_cat, 'r') as cat_in:
        skip = cat_in.readline()
        for line in cat_in:
            trans = line.split()[0]
            if trans not in deleted:
                tot_trans.add(trans)
            if len(tot_trans) == n_trans:
                break

    cat_in.close()
    if len(tot_trans) != n_trans:
        print('Warning: A higher number than annotated transcripts was requested to simulates, only %s transcript will be simulated' %(len(tot_trans)))
    

    coverage = read_count//n_trans
    tpm = (1000000.0 * coverage) / (coverage * n_trans) # Not taking into account transcript length
    f_out = open(output, 'w')
    f_out.write('target_id\test_counts\ttpm\n')
    for trans in tot_trans:
        f_out.write(trans + '\t' + str(coverage) + '\t' + str(tpm) + '\n')
    f_out.close()

def create_expr_file_nbinom(f_cat: str, f_del: str, n_trans, nbn_known, nbp_known, nbn_novel, nbp_novel, output: str):
    deleted = set()
    with open(f_del, 'r') as del_in:
        skip = del_in.readline()
        for line in del_in:
            line = line.split()
            deleted.add(line[0])
    del_in.close()


    known_trans = set()
    with open(f_cat, 'r') as cat_in:
        skip = cat_in.readline()
        for line in cat_in:
            trans = line.split()[0]
            if trans not in deleted:
                known_trans.add(trans)
            if (len(known_trans) + len(deleted)) == n_trans:
                break

    cat_in.close()

    nb_known = numpy.random.negative_binomial(nbn_known,nbp_known,len(known_trans)).tolist()
    nb_novel = numpy.random.negative_binomial(nbn_novel,nbp_novel,len(deleted)).tolist()

    f_out = open(output, 'w')
    f_out.write('target_id\test_counts\ttpm\n')
    i_known = 0
    i_novel = 0
    tot_trans = set.union(deleted,known_trans)
    for trans in tot_trans:
        if trans in deleted:
            coverage = nb_novel[i_novel]
            i_novel += 1
        else:
            coverage = nb_known[i_known]
            i_known += 1
        if coverage == 0: # minimum one count per transcript
            coverage +=1
        tpm = (1000000.0 * coverage) / (coverage * n_trans)  
        f_out.write(trans + '\t' + str(coverage) + '\t' + str(tpm) + '\n')
    f_out.close()


def create_expr_file_sample(f_cat: str, f_del: str, ref_trans,reads, output: str, tech):
    sam_file = ''.join(output.split('.')[:-1]) + '_' + tech + '.sam'

    if tech == 'pb':
        res = subprocess.run(['minimap2', ref_trans, reads, '-x', 'map-pb',
                              '-a', '--secondary=no', '-o', sam_file
        ])
    elif tech == 'ont':
        res = subprocess.run(['minimap2', ref_trans, reads, '-x', 'map-ont',
                              '-a', '--secondary=no', '-o', sam_file
        ])
    
    if res.returncode != 0:
        logging.error('minimap2 failed with code %d' %(res.returncode))
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

    deleted = set()
    with open(f_del, 'r') as del_in:
        skip = del_in.readline()
        for line in del_in:
            line = line.split()
            deleted.add(line[0])
    del_in.close()


    known_trans = set()
    with open(f_cat, 'r') as cat_in:
        skip = cat_in.readline()
        for line in cat_in:
            trans = line.split()[0]
            if trans not in deleted:
                known_trans.add(trans)
            if (len(known_trans) + len(deleted)) == n_trans:
                break
    cat_in.close()

    print(len(deleted), len(known_trans), n_trans)

    lim_novel = (len(deleted)//3)*2
    lim_known = (len(known_trans)//3)*2

    novel_expr = expr_distr[:lim_novel]
    known_expr =  expr_distr[-lim_known:]

    print(lim_novel, len(novel_expr))
    print(lim_known, len(known_expr))

    expr_distr = expr_distr[lim_novel:lim_known]
    random.shuffle(expr_distr)

    novel_expr.extend(expr_distr[:len(deleted)-lim_novel])
    known_expr.extend(expr_distr[-(len(known_trans)-lim_known):])

    f_out = open(output, 'w')
    f_out.write('target_id\test_counts\ttpm\n')

    print(len(deleted), len(novel_expr))
    print(len(known_trans), len(known_expr))

    for i, trans in enumerate(deleted):
        coverage = novel_expr[i] 
        tpm = (1000000.0 * coverage) / (coverage * n_trans)
        f_out.write(trans + '\t' + str(coverage) + '\t' + str(tpm) + '\n')
    
    for i, trans in enumerate(known_trans):
        coverage = known_expr[i]
        tpm = (1000000.0 * coverage) / (coverage * n_trans)
        f_out.write(trans + '\t' + str(coverage) + '\t' + str(tpm) + '\n')

    f_out.close()

    

    

def pb_simulation(args):
    logging.info('***Simulating PacBio reads with IsoSeqSim')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    isoseqsim = os.path.join(src_dir, 'IsoSeqSim/bin/isoseqsim')
    util_dir = os.path.join(src_dir, 'IsoSeqSim/utilities/')
    res = subprocess.run([isoseqsim, '-g', str(args.genome),
                         '-a', str(args.gtf), '--expr', str(args.expr),
                         '--c5', os.path.join(util_dir, '5_end_completeness.PacBio-Sequel.tab'),
                         '--c3', os.path.join(util_dir, '3_end_completeness.PacBio-Sequel.tab'),
                         '-o', os.path.join(args.output, 'PacBio_simulated'),
                         '-t', os.path.join(args.output, 'PacBio_simulated.tsv'),
                         '--es 0.01731', '--ed 0.01090', '--ei 0.02204',
                         '-n', args.read_count,
                         '-m normal', '--cpu', str(args.cores)                    
    ])

    if res.returncode != 0:
        logging.error('***ERROR running IsoSeqSim, contact developers for support')
        return
    
    logging.info('***IsoSeqSim simulation done')
    return


def ont_simulation(args):
    if args.read_type == 'dRNA':
        model_name = 'human_NA12878_dRNA_Bham1_guppy'
        r_type = 'dRNA'
        uracil = True
    elif args.read_type == 'cDNA':
        model_name = 'human_NA12878_cDNA_Bham1_guppy'
        r_type = 'cDNA_1D2'
        uracil = False
    else:
        logging.error('***ERROR not valid read_type value %s' %(args.read_type))
        return
    
    src_dir = os.path.dirname(os.path.realpath(__file__))
    nanosim = os.path.join(src_dir, 'NanoSim')
    models = os.path.join(src_dir, 'NanoSim/pre-trained_models/')
    model_dir = models + model_name + '/'
    if not os.path.exists(model_dir):
        logging.info('***Untar NanoSim model')
        cwd = os.getcwd()
        os.chdir(models)
        res = subprocess.run(['/bin/tar -xzf', model_name + '.tar.gz'])
        os.chdir(cwd)
        if res.returncode != 0:
            logging.error('Unpacking NanoSim pre-trained model failed')

    logging.info('***Simulating ONT reads with NanoSim')
    cmd = [nanosim, 'transcriptome', '-rt', str(args.transcriptome),
           '-rg', str(args.genome), '-e', str(args.expr),
           '-c', str(model_dir + 'training'),
           '-o', os.path.join(args.output, 'ONT_simulated'),
           '-n', str(args.read_count), '-r', r_type,
           '-b guppy', '-t', str(args.cores), '--fastq'
    ]

    if uracil:
        cmd.append('--uracil')

    res = subprocess.run(cmd)

    if res.returncode != 0:
        logging.error('***ERROR running NanoSim, contact developers for support')
        return

    logging.info('***Renaming and counting ONT reads')
    ref_trans = set()
    ref_dict = defaultdict(lambda: str())
    with open(args.gtf, 'r') as f_in:
        for line in f_in:
            if not line.startswith('#'):
                line_split = line.split()
                feature = line_split[2]
                if feature == 'exon':
                    trans_id = line_split[line_split.index('transcript_id') + 1]
                    trans_id = trans_id.replace(';', '').replace('"', '')
                    short_id = trans_id.split('.')[0]
                    ref_trans.add(short_id)  # TODO: dont loose the whole transcript id
                    ref_dict[short_id] = trans_id
    f_in.close()

    fastqs = [os.path.join(args.output, "ONT_simulated_aligned_reads.fastq"),
              os.path.join(args.output, "ONT_simulated_unaligned_reads.fastq")]

    n_read = 0
    pair_id = []
    id_counts = defaultdict(lambda: 0)
    f_name = os.path.join(args.output, 'ONT_simulated.fastq')
    f_out = open(f_name, 'w')

    for f in fastqs:
        f_in = open(f, 'r')
        for line in f_in:
            if line.startswith('@'):
                line = line.lstrip('@')
                trans_id = line.split('_')[0]

                if trans_id not in ref_trans:
                    logging.warning(trans_id, 'was not found in the annotation')
                else:
                    trans_id = ref_dict[trans_id]
                    
                id_counts[trans_id] += 1
                read_id = 'ONT_simulated_read_{}'.format(n_read)
                n_read += 1
                pair_id.append((read_id, trans_id))

                f_out.write('@{}\n'.format(read_id))
            else:
                f_out.write(line)
    f_in.close()
    f_out.close()

    logging.info('***Saving counts and read-to-isoform files')
    f_name = os.path.join(args.output, 'ONT_simulated.read_to_isoform.tsv')
    f_out = open(f_name, 'w')

    for pair in pair_id:
        f_out.write(str(pair[0]) + '\t' + str(pair[1]) + '\n')
    f_out.close()

    f_name = os.path.join(args.output, 'ONT_simulated.isoform_counts.tsv')
    f_out = open(f_name, 'w')

    for k, v in id_counts:
        f_out.write(str(k) + '\t' + str(v) + '\n')
    f_out.close()

    logging.info('***NanoSim simulation done')
    return
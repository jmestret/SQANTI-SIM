#!/usr/bin/env python3
'''
sqanti3_sim.py

Wrapper of long-read RNA-seq simulators (NanoSim and IsoSeqSim) to simulate
controlled novelty and degradation of transcripts based on SQANTI3 structural
categories

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 19/01/2022
'''

__version__ = '0.0'

import argparse
import os
import sys
import random
import numpy
from collections import defaultdict
from src import classif_gtf
from src import pb_ont_sim
from src import sqanti3_stats

def classif(input):
    parser = argparse.ArgumentParser(prog='sqanti_sim.py classif', description='sqanti_sim.py classif parse options')
    parser.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format', required=True)
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('--min_ref_len', type=int, default=0, help='\t\tMinimum reference transcript length (default: 0 bp as in largasp challenge 1 evaluation)')
    parser.add_argument('-k', '--cores', default='1', type=int, help = '\t\tNumber of cores to run in parallel')

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print('classif mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    # Classify GTF transcripts in SQANTI3 structural categories
    trans_info = classif_gtf.classify_gtf(args)

    # Print summary table
    print('***Summary table from categorization\n')
    classif_gtf.summary_table_cat(trans_info)

    return

def preparatory(input):
    parser = argparse.ArgumentParser(prog='sqanti_sim.py preparatory', description='sqanti_sim.py preparatory parse options')

    subparsers = parser.add_subparsers(dest='subparser')

    # Arguments for each simulation mode
    parser_equal = subparsers.add_parser('equal', help = '\t\tSimulate with equal coverage for all reads')
    parser_equal.add_argument('-nt', '--n_trans', default = 10000, type=int,  help = '\t\tNumber of different transcripts to simulate')
    parser_equal.add_argument('--read_count', default = 50000, type=int,  help = '\t\tNumber of reads to simulate')

    parser_custom = subparsers.add_parser('custom', help = '\t\tSimulate with diferent negative binomial distributions for novel and known transcripts')
    parser_custom.add_argument('-nt', '--n_trans', default = 10000, type=int,  help = '\t\tNumber of different transcripts to simulate')
    parser_custom.add_argument('--nbn_known',default=50,type = float, help="Average read count per transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)")
    parser_custom.add_argument('--nbp_known',default=0.5,type = float, help="The parameter 'p' of the Negative Binomial distribution")
    parser_custom.add_argument('--nbn_novel',default=5,type = float, help="Average read count per transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)")
    parser_custom.add_argument('--nbp_novel',default=0.5,type = float, help="The parameter 'p' of the Negative Binomial distribution")
    

    parser_sample = subparsers.add_parser('sample', help = '\t\tSimulate using a real sample')
    parser_sample.add_argument('--rt', help='reference transcripts in FASTA format (for custom simulation and/or ont reads)', type=str)
    parser_sample.add_argument('-i', '--reads', default = False,  help = '\t\tInput reads for quantification')
    #parser_sample.add_argument('--pb', action='store_true', help = '\t\tIf used the program will simulate PacBio reads with IsoSeqSim')
    #parser_sample.add_argument('--ont', action='store_true', help = '\t\tIf used the program will simulate ONT reads with NanoSim')
    

    # General argments for all modes
    parser.add_argument('--cat', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format', required=True)
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('-k', '--cores', default='1', type=int, help = '\t\tNumber of cores to run in parallel')
    parser.add_argument('--ISM', default='1000', type=int, help = '\t\tNumber of incomplete-splice-matches to delete')
    parser.add_argument('--NIC', default='1000', type=int, help = '\t\tNumber of novel-in-catalog to delete')
    parser.add_argument('--NNC', default='1000', type=int, help = '\t\tNumber of novel-not-in-catalog to delete')
    parser.add_argument('--Fusion', default='0', type=int, help = '\t\tNumber of Fusion to delete')
    parser.add_argument('--Antisense', default='0', type=int, help = '\t\tNumber of Antisense to delete')
    parser.add_argument('--GG', default='0', type=int, help = '\t\tNumber of Genic-genomic to delete')
    parser.add_argument('--GI', default='0', type=int, help = '\t\tNumber of Genic-intron to delete')
    parser.add_argument('--Intergenic', default='0', type=int, help = '\t\tNumber of Intergenic to delete')    
    parser.add_argument("--seed", "-s", help="randomizer seed [123]", default=123, type=int)

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print('preparatory mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    random.seed(args.seed)
    numpy.random.seed(args.seed)
    
    # Modify GTF
    counts_end = pb_ont_sim.simulate_gtf(args)
    print('***Summary table from GTF modification\n')
    counts_ini = defaultdict(lambda: 0, {
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
    pb_ont_sim.summary_table_del(counts_ini, counts_end)

    expression_out = os.path.join(args.dir, (args.output + '_expression.tsv'))
    deleted_out = os.path.join(args.dir, (args.output + '_deleted.txt'))

    if args.subparser == 'equal':
        pb_ont_sim.create_expr_file_fixed_count(args.cat, deleted_out, args.n_trans,
                                            args.read_count, expression_out
        )
    elif args.subparser == 'custom':
        pb_ont_sim.create_expr_file_nbinom(args.cat, deleted_out, args.n_trans,
                                            args.nbn_known, args.nbp_known, 
                                            args.nbn_novel, args.nbp_novel,
                                            expression_out
        )
    
    elif args.subparser == 'sample':
        '''
        if args.pb:
            pb_ont_sim.create_expr_file_sample(args.cat, deleted_out, 
                                            args.rt, args.reads,
                                            expression_out, 'pb'
            )
        if args.ont:
            pb_ont_sim.create_expr_file_sample(args.cat, deleted_out, 
                                            args.rt, args.reads,
                                            expression_out, 'ont'
            )
        '''
        pb_ont_sim.create_expr_file_sample(args.cat, deleted_out, 
                                           args.rt, args.reads,
                                           expression_out, 'ont'
        )

    
    else:
        print('Not valid sim mode', file=sys.stderr)


def sim(input):
    parser = argparse.ArgumentParser(prog='sqanti_sim.py sim', description='sqanti_sim.py sim parse options')

    # General argments for all modes
    parser.add_argument('--genome', default = False,  help = '\t\tReference genome FASTA')
    parser.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format', required=True)
    parser.add_argument('--rt', help='reference transcripts in FASTA format (for custom simulation and/or ont reads)', type=str)
    parser.add_argument('--expr', default = False,  help = '\t\tExpression file', required=True)
    parser.add_argument('--read_type', default = 'dRNA', type=str,  help = '\t\tRead type for NanoSim simulation')
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('-k', '--cores', default='1', type=int, help = '\t\tNumber of cores to run in parallel')  
    group = parser.add_argument_group('Sequencing technology', 'Choose PacBio or Nanopore reads')
    group.add_argument('--pb', action='store_true', help = '\t\tIf used the program will simulate PacBio reads with IsoSeqSim')
    group.add_argument('--ont', action='store_true', help = '\t\tIf used the program will simulate ONT reads with NanoSim')
    parser.add_argument("--seed", "-s", help="randomizer seed [123]", default=123, type=int)
    parser.add_argument('--read_count', default = None, type=int,  help = '\t\tNumber of reads to simulate')

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print('sim mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    random.seed(args.seed)
    numpy.random.seed(args.seed)
    
    args.output = os.path.join(args.dir, args.output)

    if args.pb:
        pb_ont_sim.pb_simulation(args)
    if args.ont:
        pb_ont_sim.ont_simulation(args)


def eval(input):
    parser = argparse.ArgumentParser(prog='sqanti_sim.py eval', description='sqanti_sim.py eval parse options')
    parser.add_argument('--isoforms', default = False,  help = '\t\tGTF with trancriptome reconstructed with your pipeline', required=True)
    parser.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format')
    parser.add_argument('--genome', default = False,  help = '\t\tReference genome FASTA')
    parser.add_argument('--deleted', default = False,  help = '\t\tFile with deleted trans')
    parser.add_argument('--cat', default = False,  help = '\t\tFile with deleted trans')
    parser.add_argument('--expr', default = False,  help = '\t\tFile with deleted trans')
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('-k', '--cores', default='1', type=int, help = '\t\tNumber of cores to run in parallel')
    parser.add_argument("--min_ref_len", type=int, default=0, help="\t\tMinimum reference transcript length (default: 0 bp as in largasp challenge 1 evaluation)")
    
    args, unknown = parser.parse_known_args(input)

    if unknown:
        print('sim mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    sqanti3_stats.sqanti3_stats(args)

#####################################
#                                   #
#               MAIN                #
#                                   #
#####################################

print(
    '''                                                                      
      _____  ____            _   _ _______ _____      _____ _____ __  __  
     / ____|/ __ \     /\   | \ | |__   __|_   _|    / ____|_   _|  \/  | 
    | (___ | |  | |   /  \  |  \| |  | |    | |_____| (___   | | | \  / | 
     \___ \| |  | |  / /\ \ | . ` |  | |    | |______\___ \  | | | |\/| | 
     ____) | |__| | / ____ \| |\  |  | |   _| |_     ____) |_| |_| |  | | 
    |_____/ \___\_\/_/    \_\_| \_|  |_|  |_____|   |_____/|_____|_|  |_| 
                                                                          
              A SIMULATOR OF CONTROLLED NOVELTY AND DEGRADATION           
                    OF TRANSCRIPTS SEQUENCED BY LONG-READS                
    '''
    )

if len(sys.argv) < 2:
    print('usage: python sqanti_sim.py <mode> --help\n', file=sys.sys.stderr)
    print('modes: classif, preparatory, sim, eval\n', file=sys.sys.stderr)
    sys.exit(1)

else:
    mode = sys.argv[1].lower()
    input = sys.argv[2:]

if mode == 'classif':
    res = classif(input)

if mode == 'preparatory':
    res = preparatory(input)

if mode == 'sim':
    res = sim(input)

if mode == 'eval':
    res = eval(input)


if mode == '--version':
	sys.stderr.write('SQANTI-SIM v0.0.0\n')
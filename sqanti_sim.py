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
from collections import defaultdict
from src import classif_modif_gtf
from src import expr_file_fixed_count
from src import pb_ont_sim

def classif():
    parser = argparse.ArgumentParser(prog='sqanti_sim.py classif', description='sqanti_sim.py classif parse options')
    parser.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format', required=True)
    parser.add_argument('--cat', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument("--min_ref_len", type=int, default=0, help="\t\tMinimum reference transcript length (default: 0 bp as in largasp challenge 1 evaluation)")
    parser.add_argument('--ISM', default='0', type=int, help = '\t\tNumber of incomplete-splice-matches to delete')
    parser.add_argument('--NIC', default='0', type=int, help = '\t\tNumber of novel-in-catalog to delete')
    parser.add_argument('--NNC', default='0', type=int, help = '\t\tNumber of novel-not-in-catalog to delete')
    parser.add_argument('--Fusion', default='0', type=int, help = '\t\tNumber of Fusion to delete')
    parser.add_argument('--Antisense', default='0', type=int, help = '\t\tNumber of Antisense to delete')
    parser.add_argument('--GG', default='0', type=int, help = '\t\tNumber of Genic-genomic to delete')
    parser.add_argument('--GI', default='0', type=int, help = '\t\tNumber of Genic-intron to delete')
    parser.add_argument('--Intergenic', default='0', type=int, help = '\t\tNumber of Intergenic to delete')
    parser.add_argument('--read_only', action='store_true', help = '\t\tIf used the program will only categorize the GTF file but skipping writing a new modified GTF')
    parser.add_argument('-k', '--cores', default='1', type=int, help = '\t\tNumber of cores to run in parallel')

    args, unknown = parser.parse_known_args()

    if unknown:
        print('classif mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    if not args.cat:
        # Classify GTF transcripts in SQANTI3 structural categories
        trans_info = classif_modif_gtf.classify_gtf(args)
        cat_out = os.path.join(args.dir, (args.output + '_categories.txt'))
        args.cat = cat_out

    if not args.read_only:
        # Modify reference GTF
        counts_end = classif_modif_gtf.simulate_gtf(args)

    if not args.cat:
        # Print summary table
        print('***Summary table from categorization\n')
        classif_modif_gtf.summary_table_cat(trans_info)

    if not args.read_only:
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
        classif_modif_gtf.summary_table_del(counts_ini, counts_end)

    return


def expr():
    parser = argparse.ArgumentParser(prog='sqanti_sim.py expr', description='sqanti_sim.py expr parse options')
    parser.add_argument('--cat', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('--deleted', default = False,  help = '\t\tFile with deleted trans', required=True)
    parser.add_argument('-nt', '--n_trans', default = 10000, type=int,  help = '\t\tNumber of different transcripts to simulate', required=True)
    parser.add_argument('-c', '--coverage', default = 5, type=int,  help = '\t\tNumber of counts to simulate for each transcript (coverage)', required=True)
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('--read_only', action='store_true', help = '\t\tIf used the program will only categorize the GTF file but skipping writing a new modified GTF')

    args, unknown = parser.parse_known_args()

    if unknown:
        print('expr mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    f_out = os.path.join(args.dir, (args.output + '_expression.tsv'))

    expr_file_fixed_count.create_expr_file(args.cat, args.deleted, args.n_trans,
                                           args.coverage, f_out
    )

    return


def sim():
    parser = argparse.ArgumentParser(prog='sqanti_sim.py sim', description='sqanti_sim.py sim parse options')
    parser.add_argument('--genome', default = False,  help = '\t\tReference genome FASTA')
    parser.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format', required=True)
    parser.add_argument('--expr', default = False,  help = '\t\tExpression file as formated in README', required=True)
    parser.add_argument('--read_count', default = 50000, type=int,  help = '\t\tNumber of reads to simulate', required=True)
    parser.add_argument('--read_type', default = 'dRNA', type=str,  help = '\t\tRead type for NanoSim simulation', required=True)
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('-k', '--cores', default='1', type=int, help = '\t\tNumber of cores to run in parallel')
    group = parser.add_argument_group('Sequencing technology', 'Choose PacBio or Nanopore reads')
    group.add_argument('--pb', action='store_true', help = '\t\tIf used the program will simulate PacBio reads with IsoSeqSim')
    group.add_argument('--ont', action='store_true', help = '\t\tIf used the program will simulate ONT reads with NanoSim')
    
    args, unknown = parser.parse_known_args()

    if unknown:
        print('sim mode unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    args.output = cat_out = os.path.join(args.dir, args.output)

    if args.pb:
        pb_ont_sim.pb_simulation(args)
    if args.ont:
        pb_ont_sim.ont_simulation(args)


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
    print('modes: classif, expr, sim, stats\n', file=sys.sys.stderr)
    sys.exit(1)

else:
    mode = sys.argv[1].lower()

if mode == 'classif':
    try:
        res = classif()
    except:
        sys.exit(1)

if mode == 'expr':
    try:
        res = expr()
    except:
        sys.exit(1)

if mode == 'sim':
    try:
        res = sim()
    except:
        sys.exit(1)



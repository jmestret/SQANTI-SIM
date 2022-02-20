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
import multiprocessing
import os
import sys
from collections import defaultdict
from tqdm import tqdm
from src import classif_modif_gtf


def classif():
    def initializer():
            global min_ref_len
            min_ref_len = args.min_ref_len

    parser = argparse.ArgumentParser(prog='sqanti_sim.py classif', description='sqanti_sim.py classif parse option')
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
        print('classif unrecognized arguments: {}\n'.format(' '.join(unknown)), file=sys.stderr)

    cat_out = os.path.join(args.dir, (args.output + '_categories.txt'))
    gtf_modif = os.path.join(args.dir, (args.output + '_modified.gtf'))
    del_trans = os.path.join(args.dir, (args.output + '_deleted.txt'))

    if not args.cat:
        # parsing transcripts from GTF
        print('***Parsing transcripts from GTF reference annotation file\n')
        trans_by_chr = classif_modif_gtf.gtf_parser(args.gtf)

        # classify transcripts
        print('***Classifying transcripts according to its SQANTI3 structural category\n')
        trans_info = defaultdict(lambda: [])

        if args.cores <= 1:
            initializer()
            for chrom in trans_by_chr:
                print(chrom)
                for record in tqdm(range(len(trans_by_chr[chrom]))):
                    trans_by_region = trans_by_chr[chrom][record]
                    tmp = classif_modif_gtf.transcript_classification(trans_by_region)
                    for k in tmp:
                        trans_info[k].extend(tmp[k])

        else: # multiprocessing
            all_regions = []
            for chrom in trans_by_chr:
                all_regions.extend(trans_by_chr[chrom])
            
            pool = multiprocessing.Pool(args.cores, initializer, ())
            tmp = pool.map(classif_modif_gtf.transcript_classification, all_regions)
            for x in tmp:
                for k, v in x.items():
                    trans_info[k].extend(v)

        # Write category file
        print("***Writting structural category file\n")
        classif_modif_gtf.write_category_file(trans_info, cat_out)
        cat_in = cat_out

    if not args.read_only:
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
        target = classif_modif_gtf.target_trans(cat_in, del_trans, counts)
        classif_modif_gtf.modifyGTF(args.gtf, gtf_modif, target)

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
        classif_modif_gtf.summary_table_del(counts_ini, counts)

    return


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
    print('modes: classif, sim, stats\n', file=sys.sys.stderr)
    sys.exit(1)

else:
    mode = sys.argv[1].lower()

if mode == 'classif':
    try:
        res = classif()
    except:
        sys.exit(1)


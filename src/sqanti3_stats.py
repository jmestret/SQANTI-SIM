#!/usr/bin/env python3
'''
sqanti3_stats.py
Generate counts for sim

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 20/02/2022
'''

import argparse
from email.policy import default
import subprocess
import os
import sys
from collections import defaultdict
import pandas
from SQANTI3.utilities.short_reads import get_TSS_bed, get_ratio_TSS
from SQANTI3.sqanti3_qc import CAGEPeak


try:
    from STAR import STARJunctionReader
except:
    print("Unable to import STARJunctionReader! Please make sure cDNA_Cupcake/sequence/ is in $PYTHONPATH.", file=sys.stderr)
    sys.exit(-1)

try:
    from bx.intervals import IntervalTree
except ImportError:
    print('Unable to import bx-python! Please make sure bx-python is installed.', file=sys.stderr)
    sys.exit(-1)


def sqanti3_stats(args):
    def write_whithin_cage(row):
        return within_cage_dict[row['transcript_id']]

    def write_dist_cage(row):
        return dist_cage_dict[row['transcript_id']]

    def write_ratio_TSS(row):
        if row['transcript_id'] in ratio_TSS_dict:
            return ratio_TSS_dict[row['transcript_id']]['max_ratio_TSS']
        else:
            return 1

    print('***Running SQANTI3')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    sqanti3 = os.path.join(src_dir, 'SQANTI3/sqanti3_qc.py')

    min_ref_len = 0
    cmd =[sqanti3, args.isoforms, args.gtf, args.genome,
                          '-o', args.output, '-d', args.dir, '--cpus', str(args.cores),
                          '--min_ref_len', str(min_ref_len),
                          '--force_id_ignore']
    
    if args.cage_peak:
        cmd.append('--cage_peak')
        cmd.append(args.cage_peak)

    if args.short_reads:
        cmd.append('--short_reads')
        cmd.append(args.short_reads)
        
    cmd = ' '.join(cmd)
    if subprocess.check_call(cmd, shell=True)!=0:
        print('ERROR running SQANTI3: {0}'.format(cmd), file=sys.stderr)
        #sys.exit(1)

    trans_index = pandas.read_csv(args.trans_index, sep='\t', header=0)
    if args.cage_peak:
        print('***Parsing CAGE Peak data')
        cage_peak_data = CAGEPeak(args.cage_peak)

        within_cage_dict = defaultdict(lambda: False)
        dist_cage_dict = defaultdict(lambda: False)
        with open(args.trans_index, 'r') as index_file:
            header_names = index_file.readline()
            header_names = header_names.split()
            id_pos = header_names.index('transcript_id')
            chrom_pos = header_names.index('chrom')
            strand_pos = header_names.index('strand')
            start_pos = header_names.index('TSS_genomic_coord') # start and end coordinates already swapped for negative strand
            for line in index_file:
                line = line.split()
                within_cage, dist_cage = cage_peak_data.find(line[chrom_pos], line[strand_pos], line[start_pos])
                within_cage_dict[line[id_pos]] = within_cage
                dist_cage_dict[line[id_pos]] = dist_cage
        index_file.close()

        trans_index['dist_to_cage_peak'] = trans_index.apply(write_dist_cage, axis=1)
        trans_index['within_cage_peak'] = trans_index.apply(write_whithin_cage, axis=1)
    
    if args.short_reads:
        # Short Read Coverage
        # TODO: get min_cov and min_cov_pos
        
        # Short reads ratio TSS
        star_out = os.path.join(args.dir, '/STAR_mapping/')
        star_index = os.path.join(args.dir, '/STAR_index/')
        chr_order = star_index + "/chrNameLength.txt"
        inside_bed, outside_bed = get_TSS_bed(args.gtf, chr_order)
        bams=[]
        for filename in os.listdir(star_out):
            if filename.endswith('.bam'):
                bams.append(star_out + '/' + filename)
        ratio_TSS_dict = get_ratio_TSS(inside_bed, outside_bed, bams, chr_order)
        trans_index['ratio_TSS'] = trans_index.apply(write_ratio_TSS, axis=1)


    trans_index.to_csv(args.trans_index, sep='\t', na_rep='NA', header=True, index=False)

    print('***Generating SQANTI-SIM report')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    classification_file = os.path.join(args.dir, (args.output + '_classification.txt'))
    junctions_file = os.path.join(args.dir, (args.output + '_junctions.txt'))

    cmd=['Rscript', os.path.join(src_dir,'SQANTI_SIM_report.R'),
         classification_file, junctions_file, args.trans_index, src_dir]

    cmd = ' '.join(cmd)
    if subprocess.check_call(cmd, shell=True)!=0:
        print('ERROR running SQANTI-SIM report generation: {0}'.format(cmd), file=sys.stderr)
        sys.exit(1)
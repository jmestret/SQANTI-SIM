#!/usr/bin/env python3
'''
sqanti3_stats.py
Generate counts for sim

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 20/02/2022
'''

import argparse
import subprocess
import os
import sys

def sqanti3_stats(args):
    print('***Running SQANTI3')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    sqanti3 = os.path.join(src_dir, 'SQANTI3/sqanti3_qc.py')

    cmd =[sqanti3, args.isoforms, args.gtf, args.genome,
                          '-o', args.output, '-d', args.dir, '--cpus', str(args.cores),
                          '--min_ref_len', str(args.min_ref_len),
                          '--force_id_ignore']
    
    
    if subprocess.check_call(cmd, shell=True)!=0:
        print('ERROR running SQANTI3: {0}'.format(cmd), file=sys.stderr)
        sys.exit(1)

    print('***Generating SQANTI-SIM report')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    classification_file = os.path.join(args.dir, (args.output + '_classification.txt'))
    junctions_file = os.path.join(args.dir, (args.output + '_junctions.txt'))

    cmd=['Rscript', os.path.join(src_dir,'SQANTI_SIM_report.R'),
         classification_file, junctions_file, args.trans_index, src_dir]

    res = subprocess.run()

    if subprocess.check_call(cmd, shell=True)!=0:
        print('ERROR running SQANTI-SIM report generation: {0}'.format(cmd), file=sys.stderr)
        sys.exit(1)
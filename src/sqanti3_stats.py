#!/usr/bin/env python3
'''
sqanti3_stats.py
Generate counts for sim

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 20/02/2022
'''

import argparse
import subprocess
import logging
import os

def sqanti3_stats(args):
    logging.info('***Running SQANTI3')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    sqanti3 = os.path.join(src_dir, 'SQANTI3/sqanti3_qc.py')

    cmd =[sqanti3, args.isoforms, args.gtf, args.genome,
                          '-o', args.output, '-d', args.dir, '--cpus', str(args.cores),
                          '--force_id_ignore']
    if args.min_ref_len != 0:
        cmd.append('--min_ref_len')
        cmd.append(str(args.min_ref_len))
    
    res = subprocess.run(cmd)

    if res.returncode != 0:
        logging.error('***ERROR running SQANTI3, contact developers for support')
        return
    
    logging.info('***Generating SQANTI-SIM report')

    src_dir = os.path.dirname(os.path.realpath(__file__))
    classification_file = os.path.join(args.dir, (args.output + '_classification.txt'))
    junctions_file = os.path.join(args.dir, (args.output + '_junctions.txt'))

    res = subprocess.run([os.path.join(src_dir,'SQANTI_SIM_report.R'),
                          classification_file, junctions_file, args.deleted,
                          args.cat, args.expr, src_dir
    ])
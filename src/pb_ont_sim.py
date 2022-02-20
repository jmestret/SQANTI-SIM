'''
Generate counts for sim

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 19/01/2020
'''

from asyncio.log import logger
import logging
import os
from pyexpat import model
import subprocess
from weakref import ref


def pb_simulation(args):
    logging.info('***Simulating PacBio reads with IsoSeqSim')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    isoseqsim = os.path.join(src_dir, 'IsoSeqSim/bin/isoseqsim')
    util_dir = os.path.join(src_dir, 'IsoSeqSim/utilities/')
    res = subprocess.run([isoseqsim, '-g', str(args.genome),
                         '-a', str(args.annot), '--expr', str(args.expr),
                         '--c5', os.path.join(util_dir, '5_end_completeness.PacBio-Sequel.tab'),
                         '--c3', os.path.join(util_dir, '3_end_completeness.PacBio-Sequel.tab'),
                         '-o', os.path.join(args.output, 'PacBio_simulated'),
                         '-t', os.path.join(args.output, 'PacBio_simulated.tsv'),
                         '--es 0.01731', '--ed 0.01090', '--ei 0.02204',
                         '--read_numer', str(args.read_count / 1000000.0),
                         '-m normal', '--cpu', str(args.cpus)                    
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
        res = subprocess.run(['tar -xzf', model_name + '.tar.gz'])
        os.chdir(cwd)
        if res.returncode != 0:
            logger.error('Unpacking NanoSim pre-trained model failed')

    logging.info('***Simulating ONT reads with NanoSim')
    cmd = [nanosim, 'transcriptome', '-rt', str(args.transcriptome),
           '-rg', str(args.genome), '-e', str(args.expr),
           '-c', str(model_dir + 'training'),
           '-o', os.path.join(args.output, 'ONT_simulated'),
           '-n', str(args.read_count), '-r', r_type,
           '-b guppy', '-t', str(args.cpus), '--fastq'
    ]

    if uracil:
        cmd.append('--uracil')

    res = subprocess.run(cmd)

    if res.returncode != 0:
        logging.error('***ERROR running NanoSim, contact developers for support')
        return

    logger.info('***Renaming and counting ONT reads')
    ref_trans = set()
    with open(args.annot, 'r') as f_in:
        for line in f_in:
            if not line.startswith('#'):
                line_split = line.split()
                feature = line_split[2]
                if feature == 'exon':
                    trans_id = line_split[line_split.index('transcript_id') + 1]
                    trans_id = trans_id.replace(';', '').replace('"', '')
                    ref_trans.add(trans_id)
    f_in.close()

    fastqs = [os.path.join(args.output, "ONT_simulated_aligned_reads.fastq"),
              os.path.join(args.output, "ONT_simulated_unaligned_reads.fastq")]

    n_reads = 0
    f_name = os.path.join(args.output, 'ONT_simulated.fastq')
    f_out = open(f_name, 'w')

    for f in fastqs:
        f_in = open(f, 'r')
        for line in f_in:
            if line.startswith('@'):
                line = line.lstrip('@')
                trans_id = line.split('_')[0]

    # TODO: end changing names

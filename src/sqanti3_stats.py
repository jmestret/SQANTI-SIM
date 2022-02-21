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
from collections import defaultdict


class myQueryIsoforms:
    '''Features of the query isoform and its associated reference'''
    def __init__(self, id=None, gene_id=None, str_class=None, genes=None, transcripts=None, junctions=set(), start = None, end = None, names=[], counts=0):
        self.id = id
        self.gene_id = gene_id
        self.str_class = str_class
        self.genes = genes
        self.transcripts = transcripts
        self.junctions = junctions
        self.start = start
        self.end = end
        self.names = names
        self.counts = counts
        

def run_sqanti3(args):
    logging.info('***Running SQANTI3')
    src_dir = os.path.dirname(os.path.realpath(__file__))
    sqanti3 = os.path.join(src_dir, 'SQANTI3/sqanti3_lrgasp.challenge1.py')

    cmd =[sqanti3, args.isoforms, args.gtf, args.genome,
                          '-o', args.output, '-d', args.dir, '--cpus', args.cores,
                          '--gtf', '--force_id_ignore']
    if args.min_ref_len != 0:
        cmd.append('--min_ref_len')
        cmd.append(args.min_ref_len)
    
    res = subprocess.run(cmd)

    if res.returncode != 0:
        logging.error('***ERROR running SQANTI3, contact developers for support')
        return
    
    logging.info('***SQANTI3 quality control done')

def stats(args, classification_file, junctions_file):
    # Get junctions from deleted reads
    ref_by_SC = defaultdict(lambda: [])
    
    with open(args.deleted, 'r') as f_del:
        skip = f_del.readline()
        for line in f_del:
            juncs = set()
            line = line.split()
            SC = line[2]
            TSS = line[5]
            TTS = line[6]

            donors = line[7].split(',')
            acceptors = line[8].split(',')
            if isinstance(donors[0], str):
                juncs = set()
            else:
                for d, a in zip(donors, acceptors):
                    d = int(d) + 1
                    a = int(a) -1
                    juncs.add(str(d))
                    juncs.add(str(a))
            
            ref_by_SC[SC].append(myQueryIsoforms(id=line[0], gene_id=line[1],
                                       str_class=line[2],
                                       genes=line[3].split('_'),
                                       transcripts=line[4].split('_'),
                                       junctions=juncs,
                                       start=TSS, end=TTS))
    f_del.close()

    # Count ocurrencies, no need if NanoSim simulate the number of reads that you damm ask
    '''
    if args.nanosim:
        with open(args.fsim, 'r') as f_sim:
            for line in f_sim:
                if line.startswith('@'):
                    line = line.lstrip('@')
                    id = line.split('_')[0]

                    for SC in ref_by_SC:
                        for i in range(len(ref_by_SC[SC])):
                            if id == ref_by_SC[SC][i].id.split('.')[0]:
                                ref_by_SC[SC][i].names.append(line)
                                ref_by_SC[SC][i].counts += 1
        f_sim.close()



    elif args.isoseqsim:
        with open(args.fsim, 'r') as f_sim:
            for line in f_sim:
                    line = line.split()
                    simid = line[0]
                    refid = line[1]

                    for SC in ref_by_SC:
                        for i in range(len(ref_by_SC[SC])):
                            if refid == ref_by_SC[SC][i].id:
                                ref_by_SC[SC][i].names.append(simid)
                                ref_by_SC[SC][i].counts += 1
        f_sim.close()  
    

    # Delete those simulated with low coverage (smaller than the threshhold used in the pipeline)
    threshold = 3
    for SC in ref_by_SC:
        print(SC, len(ref_by_SC[SC]))
        for rec in ref_by_SC[SC]:
            if rec.counts < threshold:
                ref_by_SC[SC].remove(rec)
        print(SC, len(ref_by_SC[SC]))
    '''


    # READ read to isoform id by the pipeline file!
    # for talon is read_annot.tsv       

    # Get SC and ref from query isoforms
    isos = defaultdict(lambda: myQueryIsoforms())
    with open(classification_file, 'r') as f_class:
        skip = f_class.readline()
        for line in f_class:
            line = line.split()
            iso = line[0]
            SC = line[5]
            ref_g = line[6].split('_')
            ref_t = line[7].split('_')
            TSS = line[47]
            TTS = line[48]
            isos[iso] = myQueryIsoforms(id=iso, gene_id=None,
                                                 str_class=SC,
                                                 genes=ref_g,
                                                 transcripts=ref_t,
                                                 start=TSS, end=TTS)

    f_class.close()
    
    # Get junctions from query isoforms
    juncs = defaultdict(lambda: set())
    with open(junctions_file, 'r') as f_junc:
        skip = f_junc.readline()
        for line in f_junc:
            line = line.split()
            iso = line[0]
            d = line[4]
            a = line[5]
            juncs[iso].add(d)
            juncs[iso].add(a)
    f_junc.close()

    iso_by_SC = defaultdict(lambda: [])
    for id, iso in isos.items():
        iso.junctions = juncs[id]
        iso_by_SC[iso.str_class].append(iso)

    # Get Stats
    metrics = defaultdict(lambda: [])
    for SC in ref_by_SC:
        TP = 0
        FP = 0
        for iso in iso_by_SC[SC]:            
            if len(iso.junctions) == 0:
                for rec in ref_by_SC[SC]:
                    if rec.junctions == iso.junctions and rec.start <= iso.end and iso.start <= rec.end:
                        TP += 1
                        break
                else:
                    FP += 1
            else:
                for rec in ref_by_SC[SC]:
                    if rec.junctions == iso.junctions:
                        TP += 1
                        break
                else:
                    FP += 1
        FN = len(ref_by_SC[SC]) - TP
        metrics[SC] = [str(TP), str(FN), str(FP), str(round(TP/(TP+FP),2)),
                       str(round(FP/(TP+FP),2)), str(round(TP/(TP+FN),2))]
        
    print('\033[94m_' * 79 + '\033[0m')
    print('\033[92mS Q A N T I - S I M\033[0m \U0001F4CA')
    print()
    print('Performance metrics \U0001F50E')
    print('\033[94m_' * 79 + '\033[0m')
    
    d_keys = list(metrics.keys())
    stats = ['TP', 'FN', 'FP', 'Precision', 'FDR', 'Sensitivity']
    print('\033[92m|\033[0m SC:', end='')
    print(*d_keys, sep='\t')
    for i in range(len(metrics[SC])):
        l = []
        for k in d_keys:
            l.append(metrics[k][i])
        
        print('\033[92m|\033[0m ' + stats[i] + ': ', end='')
        print(*l, sep='\t')
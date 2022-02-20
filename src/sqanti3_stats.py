#!/usr/bin/env python3
'''
Stats from sqantisim
Author: Jorge Mestre Tomas
Date: 15/02/2021
Last update: 15/02/2022
'''

import argparse
from collections import defaultdict
from operator import ge
from termios import TCSETS
from tracemalloc import start


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
        


def main():
    parser = argparse.ArgumentParser(prog='sqanti3_sim.py', description="SQANTI-SIM: a simulator of controlled novelty and degradation of transcripts sequence by long-reads")
    parser.add_argument('--classi', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('--junc', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('--deleted', default = False,  help = '\t\tFile with deleted trans', required=True)
    parser.add_argument('--fsim', default = False,  help = '\t\tFile with deleted trans', required=True)
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPath for output files')
    parser.add_argument('--nanosim', action='store_true', help = '\t\tIf used the program will only categorize the GTF file but skipping writing a new modified GTF')
    parser.add_argument('--isoseqsim', action='store_true', help = '\t\tIf used the program will only categorize the GTF file but skipping writing a new modified GTF')
    
    args = parser.parse_args()

    # Get junctions from deleted reads
    #ref_by_chr = defaultdict(lambda: myQueryIsoforms())
    #ref = []
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

    # Count ocurrencies
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


    # READ read to isoform id by the pipeline file!
    # for talon is read_annot.tsv       

    # Get SC and ref from query isoforms
    isos = defaultdict(lambda: myQueryIsoforms())
    with open(args.classi, 'r') as f_class:
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
    with open(args.junc, 'r') as f_junc:
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
    for SC in ref_by_SC:
        print(SC, 'stats')
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
        print('TP', TP)
        print('FP', FP)
        print('FN', FN)
        print('Precision', TP/(TP+FP))
        print('FDR', FP/(TP+FP))
        print('Sensitivity', TP/(TP+FN))



if __name__ == '__main__':
    main()

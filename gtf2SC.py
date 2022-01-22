#!/usr/bin/env python3
'''
gtf2SC.py
Classify transcripts in SQANTI3 SC if potentially deleted from GTF.
Given a GTF file as input, determine its potential SQANTI3 structural
category not taking into account himself in the reference.

Author: Jorge Mestre Tomas
Date: 19/01/2020
Last update: 19/01/2021 by Jorge Mestre
'''

__author__ = 'jormart2@alumni.uv.es'
__version__ = '0.0'

import os
import argparse
from time import time

#####################################
#                                   #
#          DEFINE CLASSES           #
#                                   #
#####################################

class summary_table:
    '''This objects aims to output a summary table of the characterization'''
    counts = {
            'FSM':0,
            'ISM':0,
            'NIC':0,
            'NNC':0,
            'Fusion':0,
            'Antisense':0,
            'Genic-genomic':0,
            'Genic-intron':0,
            'Intergenic':0,
            'Unclassified':0
        }

    def addCounts(self, d_gene):
        '''Add counts of each SC to the table'''
        for chr in d_gene.values():
            for gen in chr:
                for trans in gen.transcripts:
                    if trans.SC in list(self.counts.keys()):
                        self.counts[trans.SC] += 1
                    else:
                        self.counts[trans.SC] = 1
    
    def __repr__(self):
        # TODO: print the object properly
        return 'summary_table()'

    def __str__(self):
        print('\033[94m_' * 79 + '\033[0m')
        print('\033[92mS Q A N T I - S I M\033[0m \U0001F4CA')
        print()
        print('Summary Table \U0001F50E')
        print('\033[94m_' * 79 + '\033[0m')
        for k, v in self.counts.items():
            print('\033[92m|\033[0m ' + k + ': ' + str(v))

class sequence:
    def __init__(self, id, start, end):
        self.id = id
        self.start = start
        self.end = end
        self.genes = list() 

class gene:
    def __init__(self, id, chr, strand, start, end):
        self.id = id
        self.chr = chr
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = list()

class transcript:
    def __init__(self, id, gene_id, strand, exon_coords):
        self.id = id
        self.gene_id = gene_id 
        self.strand = strand
        self.TSS = exon_coords[0]
        self.TTS = exon_coords[len(exon_coords)-1]
        #self.SJ = self.getSJ()
        #self.exons = self.getExon(exon_coords)
        self.SC = 'Unclassified'
        self.ref = str()
        self.diff_TSS = 0
        self.diff_TTS = 0


#####################################
#                                   #
#         DEFINE FUNCTIONS          #
#                                   #
#####################################

def readgtf(gtf: str)-> list:
    '''
    Given a GTF file name it save all the transcripts classified by gene
    and region in the genome.

    Args:
        gtf (str) GTF file name

    Returns:
        res (list) list of "seqname" objetcs containing a list of "gene" objects
                   with all "transcripts"
    '''
    
    l_coords = list()
    res = list()
    prev_seqname = None
    prev_gene = None
    prev_trans = None
    prev_strand = None

    # Read GTF file line by line
    with open(gtf, 'r') as f_in:
        for line in f_in:
            if not line.startswith('#'):
                line_split = line.split()
                feature = line_split[2]

                # Get only features that are 'exon'
                if feature == 'exon':
                    seqname_id = line_split[0]
                    gene_id = line_split[line_split.index('gene_id') + 1]
                    gene_id = gene_id.replace(';', '').replace('"', '')
                    trans_id = line_split[line_split.index('transcript_id') + 1]
                    trans_id = trans_id.replace(';', '').replace('"', '')
                    
                    # Swap coordinates if negative strand
                    strand = line_split[6]
                    if strand == '+':
                        start = int(line_split[3])
                        end = int(line_split[4])
                    else:
                        start = int(line_split[4])
                        end = int(line_split[3])

                    # Reading first exon
                    if not prev_trans:
                        prev_trans = trans_id
                        prev_gene = gene_id
                        prev_seqname = seqname_id
                        prev_strand = strand
                        l_coords = [start, end]
                        g = gene(gene_id, seqname_id, strand, 0, 0) # TODO: start and end
                        region = sequence(seqname_id, 0, 0) # TODO: start and end
                        
                    # If reading same transcript add exon start and end
                    elif prev_trans == trans_id:
                        l_coords.append(start)
                        l_coords.append(end)
                    
                    else:
                        t = transcript(prev_trans, prev_gene, prev_strand, l_coords)
                        g.transcripts.append(t)
                        prev_trans = trans_id
                        l_coords = [start, end]
    
                        if prev_gene != gene_id:
                            region.genes.append(g)
                            g = gene(gene_id, seqname_id, strand, 0, 0) # TODO: start and end
                            prev_gene = gene_id

                            if prev_seqname != seqname_id:
                                res.append(region)
                                region = sequence(seqname_id, 0, 0) # TODO: start and end
                                prev_seqname = seqname_id

    # Save last transcript    
    t = transcript(prev_trans, prev_gene, prev_strand, l_coords)
    g.transcripts.append(t)
    region.genes.append(g)
    res.append(region)

    return res




def main():
    os.chdir('/home/jorge/Desktop')

    f_name = '/home/jorge/Desktop/simulacion/getSC/chr3.gencode.v38.annotation.gtf'

    # Read GTF file
    data = readgtf(f_name)



    




if __name__ == '__main__':
    t_ini = time()
    main()
    t_fin = time()
    print('[Execution time %s seconds]' %(t_fin-t_ini))
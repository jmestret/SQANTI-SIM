#!/usr/bin/env python3
'''
gtf2SC.py
Classify transcripts in SQANTI3 SC if potentially deleted from GTF.
Given a GTF file as input, determine its potential SQANTI3 structural
category not taking into account himself in the reference.

Author: Jorge Mestre Tomas
Date: 19/01/2020
Last update: 23/01/2021 by Jorge Mestre
'''

__author__ = 'jormart2@alumni.uv.es'
__version__ = '0.0'

from operator import le
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
            'GeneOverlap':0,
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
        self.SJ = self.getSJ(exon_coords)
        #self.exons = self.getExon(exon_coords)
        self.SC = 'Unclassified'
        self.ref = str()
        self.diff_TSS = 0
        self.diff_TTS = 0
        self.gene_hits = []
        self.trans_hits = []

    def getSJ(self, exon_coords):
        '''
        Define the splice junctions
        '''
        splice_sites = exon_coords[1:len(exon_coords)-1]
        if len(self.splice_sites) > 0:
            SJs = []
            for i in range(0, len(self.splice_sites), 2):
                SJs.append((self.splice_sites[i], self.splice_sites[i+1]))
            return(SJs)
        else:
            return([])

    def classify_trans(self, ref):
        '''
        Given a group of reference trasncripts elucidate the SC of the
        target transcript
        '''
        # TODO: implement for sequence and gene level (Faster but no that accuarate)

        #----------------------------------#
        #                                  #
        #      UNSPLICED TRANSCRIPT        #
        #                                  #
        #----------------------------------#
        if self.is_monoexon():
            self.get_trans_hits()

            # hits any monoexon?
            hits_monoexon = False
            for hit_index in self.trans_hits:
                hit = ref[hit_index[0]][hit_index[1]]
                if len(hit.SJ) == 0:
                    hits_monoexon = True
                    #break
            
            if hits_monoexon:
                # 1) If completly within the ref TSS and TTS: FSM
                # 2) If not hitting at all: Intergenic (TODO: propably this is useless)
                # 3) If partially hitting: Genic-genomic
                for hit_index in self.trans_hits:
                    hit = ref[hit_index[0]][hit_index[1]]
                    if len(hit.SJ) == 0:
                        if self.strand != hit.strand:
                            self.eval_new_SC('Antisense')

                        elif self.TSS >= hit.TSS and self.TTS <= hit.TTS:
                            self.SC = 'FSM'

                        elif self.TTS <= hit.TSS or self.TSS >= hit.TTS:
                            self.eval_new_SC('Intergenic')
                        
                        elif (self.TSS <= hit.TSS and self.TTS > hit.TSS) or \
                             (self.TSS < hit.TTS and self.TTS >= hit.TTS):
                             self.eval_new_SC('Genic-genomic')
                return

            else:
                for hit_index in self.trans_hits:
                    hit = ref[hit_index[0]][hit_index[1]]
                    if self.is_within_intron(hit):
                        self.eval_new_SC('Genic-intron')
                    pass
        
        #----------------------------------#
        #                                  #
        #       SPLICED TRANSCRIPT         #
        #                                  #
        #----------------------------------#
        else:
            pass # TODO


    def eval_new_SC(self, new_SC):
        SCrank ={
            'FSM':1, 'ISM':2, 'Fusion':3,
            'NIC': 4, 'NNC':5, 'Antisense': 6,
            'Genic-genomic':7, 'Genic-intron':8, 'Intergenic':9,
            'GeneOverlap':10, 'Unclassified':11
        }
        
        if SCrank[self.SC] > SCrank[new_SC]:
            self.SC = new_SC

    def is_monoexon(self):
        if len(self.SJ) == 0:
            return True
        return False
    
    def get_trans_hits(self, ref):
        for g_index, g in enumerate(ref.genes):
            for t_index, t in enumerate(g.transcripts):
                # Overlaping transcripts?
                if (self.TSS <= t.TSS and self.TTS > t.TSS) or \
                   (self.TSS < t.TTS and self.TTS >= t.TTS) or \
                   (self.TSS >= t.TSS and self.TTS <= t.TTS):

                   self.trans_hits.append((g_index, t_index))

    def is_within_intron(self, trans):
        for i in trans.SJ:
            if self.TSS > i[0] and self.TTS < i[1]:
                return True
        return False



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
    gene_start = None
    gene_end = None

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
                        gene_start = start
                        gene_end = end
                        l_coords = [start, end]
                        g = gene(gene_id, seqname_id, strand, 0, 0) # TODO: start and end
                        region = sequence(seqname_id, 0, 0) # TODO: start and end
                        
                    # If reading same transcript add exon start and end
                    elif prev_trans == trans_id:
                        l_coords.append(start)
                        l_coords.append(end)
                        gene_end = end
                    
                    else:
                        t = transcript(prev_trans, prev_gene, prev_strand, l_coords)
                        g.transcripts.append(t)
                        prev_trans = trans_id
                        l_coords = [start, end]
    
                        if prev_gene != gene_id:
                            region.genes.append(g)
                            g = gene(gene_id, seqname_id, strand, gene_start, gene_end) # TODO: start and end
                            prev_gene = gene_id
                            gene_start = start

                            if prev_seqname != seqname_id:
                                res.append(region)
                                region = sequence(seqname_id, 0, 0) # TODO: start and end
                                prev_seqname = seqname_id
                        gene_end = end

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

    # Classify transcripts in each different sequence
    for seq in range(len(data)):

        pass



    




if __name__ == '__main__':
    t_ini = time()
    main()
    t_fin = time()
    print('[Execution time %s seconds]' %(t_fin-t_ini))
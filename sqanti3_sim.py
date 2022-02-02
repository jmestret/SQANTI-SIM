#!/usr/bin/env python3
'''
sqanti3_sim.py
Classify transcripts in SQANTI3 SC if potentially deleted from GTF.
Given a GTF file as input, determine its potential SQANTI3 structural
category not taking into account himself in the reference.
Modify original GTF deleting transcripts to simulate reads.

Author: Jorge Mestre Tomas
Date: 19/01/2020
Last update: 02/02/2021 by Jorge Mestre
'''

__author__ = 'jormart2@alumni.uv.es'
__version__ = '0.0'

import os
import copy
import argparse
import random
import subprocess
from time import time
from tqdm import tqdm


#####################################
#                                   #
#          DEFINE CLASSES           #
#                                   #
#####################################
class summary_table:
    '''
    This objects aims to output a summary table of the characterization
    '''

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


    def addCounts(self, data: list):
        '''
        Add counts of each SC to the table
        '''

        for seq in data:
            for g in seq.genes:
                for t in g.transcripts:
                    if t.SC in list(self.counts.keys()):
                        self.counts[t.SC] += 1
                    else:
                        self.counts[t.SC] = 1


    def __str__(self):
        print('\033[94m_' * 79 + '\033[0m')
        print('\033[92mS Q A N T I - S I M\033[0m \U0001F4CA')
        print()
        print('Summary Table \U0001F50E')
        print('\033[94m_' * 79 + '\033[0m')
        for k, v in self.counts.items():
            print('\033[92m|\033[0m ' + k + ': ' + str(v))
        return ''


class sequence:
    '''
    Represents a region of a seqname_id of a GTF file where genes may overlap
    '''
    def __init__(self, id, start: int, end: int):
        self.seqname = id
        self.start = start
        self.end = end
        self.genes = list() 

    def classify_trans(self):
        '''
        Categorize all transcripts from the genes in this sequence according
        to its SQANTI3 structural category
        '''
        for g_index in range(len(self.genes)):
            for t_index in range(len(self.genes[g_index].transcripts)):
                ref = copy.deepcopy(self)
                del ref.genes[g_index].transcripts[t_index]
                if len(ref.genes[g_index].transcripts) == 0:
                    del ref.genes[g_index]
                self.genes[g_index].transcripts[t_index].get_SC(ref)


class gene:
    '''
    Saves the info of a given gene (unique id)
    '''
    def __init__(self, id: str, chr: str, strand: str, start: int, end: int):
        self.id = id
        self.chr = chr
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = list()


class transcript:
    '''
    Collect all the relevant information of a transcript to classify it properly
    '''
    def __init__(self, id: str, gene_id: str, strand: str, exon_coords: list):
        self.id = id
        self.gene_id = gene_id 
        self.strand = strand
        self.TSS = exon_coords[0]
        self.TTS = exon_coords[len(exon_coords)-1]
        self.SJ = self.getSJ(exon_coords)
        self.SC = 'Unclassified'
        self.subtype = ''
        self.match_type = ''
        self.ref_trans = str()
        self.ref_gene = str()
        self.diff_TSS = None
        self.diff_TTS = None
        self.gene_hits = set()
        self.trans_hits = []


    def getSJ(self, exon_coords: list)-> list:
        '''
        Define the splice junctions of the transcripts

        Args:
            exon_coords (list) A list with the TSS, the splice sites and the TTS
        
        Returns:
            SJs (list) a list of tupples with the splice junctions
        '''

        splice_sites = exon_coords[1:len(exon_coords)-1]
        if len(splice_sites) > 0:
            SJs = []
            for i in range(0, len(splice_sites), 2):
                SJs.append((splice_sites[i], splice_sites[i+1]))
            return(SJs)
        else:
            return([])
            

    def get_SC(self, ref: sequence):
        '''
        Given a group of reference trasncripts (sequence class) elucidate the SC
        of the target transcript

        Args:
            self (transcript) the target transcript to find its SC
            ref (sequence) a region of a sequence with all overlapping genes
                and its transcripts
        '''

        # If no overlaping genes don't continue, it just can be Intergenic
        if len(ref.genes) == 0:
            self.SC = 'Intergenic'
            return

        #----------------------------------#
        #                                  #
        #      UNSPLICED TRANSCRIPT        #
        #                                  #
        #----------------------------------#
        if self.is_monoexon():
            self.get_trans_hits(ref)

            # hits any monoexon?
            hits_monoexon = False
            for hit_index in self.trans_hits:
                hit = ref.genes[hit_index[0]].transcripts[hit_index[1]]
                if len(hit.SJ) == 0:
                    hits_monoexon = True
                    #break
            
            if hits_monoexon:
                # 1) If completly within the ref TSS and TTS: FSM
                # 2) If not hitting at all: Intergenic (TODO: propably this is useless)
                # 3) If partially hitting: Genic-genomic
                for hit_index in self.trans_hits:
                    hit = ref.genes[hit_index[0]].transcripts[hit_index[1]]
                    if len(hit.SJ) == 0:
                        if self.strand != hit.strand and self.hit_exon(hit):
                            self.eval_new_SC('Antisense', hit)

                        elif self.TSS >= hit.TSS and self.TTS <= hit.TTS:
                            self.eval_new_SC('FSM', hit)

                        elif self.TTS <= hit.TSS or self.TSS >= hit.TTS:
                            self.eval_new_SC('Intergenic')
                        
                        elif (self.TSS <= hit.TSS and self.TTS > hit.TSS) or \
                             (self.TSS < hit.TTS and self.TTS >= hit.TTS):
                             self.eval_new_SC('Genic-genomic', hit)

            else:
                for hit_index in self.trans_hits:
                    hit = ref.genes[hit_index[0]].transcripts[hit_index[1]]

                    if self.strand != hit.strand and self.hit_exon(hit):
                        self.eval_new_SC('Antisense', hit)
                    
                    if self.strand == hit.strand:
                        if self.is_within_intron(hit):
                            self.eval_new_SC('Genic-intron', hit)

                        elif self.is_within_exon(hit):
                            self.eval_new_SC('ISM', hit)

                        elif hit.TSS <= self.TSS < self.TTS <= hit.TTS:
                            if hit.in_intron(self.TSS) or hit.in_intron(self.TTS):
                                self.eval_new_SC('Genic-genomic', hit)
                            else:
                                self.eval_new_SC('NIC', hit)
                        else:
                            self.eval_new_SC('Genic-genomic', hit)
                
                if self.SC == 'Unclassified':
                    self.eval_new_SC('Intergenic')
        
        #----------------------------------#
        #                                  #
        #       SPLICED TRANSCRIPT         #
        #                                  #
        #----------------------------------#
        else:
            for g_index, g in enumerate(ref.genes):
                for t_index, trans in enumerate(g.transcripts):
                    if len(trans.SJ) == 0:
                        self.eval_new_SC('GeneOverlap', trans)
                        self.gene_hits.add(g_index) # TODO: this generates a lot of fusion genes

                    else:
                        match_type = self.compare_junctions(trans)

                        if match_type == 'exact':
                            self.eval_new_SC('FSM', trans)
                        
                        elif match_type == 'subset':
                            self.eval_new_SC('ISM', trans)
                        
                        elif match_type == 'partially':
                            self.match_type = 'partially'
                            self.gene_hits.add(g_index)
                            self.trans_hits.append((g_index, t_index))
                        
                        elif match_type == 'no_match' and not self.match_type:
                            self.match_type = match_type

            if self.SC not in ['FSM', 'ISM']:
                if self.match_type == 'partially':
                    gene_hits = list(self.gene_hits)
                    l_genes = []
                    for i in gene_hits:
                        l_genes.append(ref.genes[i])

                    #if len(gene_hits) > 1 and dont_overlap(l_genes): # TODO overlaping genes?
                    if len(gene_hits) > 1:
                        gene_names = []
                        for i in gene_hits:
                            gene_names.append(ref.genes[i].id)
                        self.eval_new_SC('Fusion', gene_names) # TODO: which reference I include?
                    
                    #elif len(gene_hits) > 1:
                    #    for hit_index in self.gene_hits:
                    #        if self.acceptor_subset(ref.genes[hit_index]) and self.donor_subset(ref.genes[hit_index]):
                    #            self.eval_new_SC('NIC', ref.genes[hit_index].transcripts[0]) # TODO: which reference I include?
                    #        else:
                    #            self.eval_new_SC('NNC', ref.genes[hit_index].transcripts[0]) # TODO: which reference I include?
                    elif len(gene_hits) == 1:
                        if self.acceptor_subset(ref.genes[gene_hits[0]]) and self.donor_subset(ref.genes[gene_hits[0]]):
                            self.eval_new_SC('NIC', ref.genes[gene_hits[0]].transcripts[0]) # TODO: which reference I include?
                        elif not self.hit_exon(trans):
                            self.eval_new_SC('Intergenic')
                        else:
                            self.eval_new_SC('NNC', ref.genes[gene_hits[0]].transcripts[0]) # TODO: which reference I include?
                    else:
                        print('Algo raro pasa aqui -.-')

                elif self.match_type == 'no_match':
                    gene_hits = self.hits_a_gene(ref)
                    if gene_hits:
                        for index in gene_hits:
                            if self.strand != ref.genes[index].strand:
                                for t in ref.genes[index].transcripts:
                                    #if self.hit_exon(t):
                                    self.eval_new_SC('Antisense', t)
                                        #break
                                    #else:
                                    #    self.eval_new_SC('Genic-genomic', t) 
                            else:
                                for t in ref.genes[index].transcripts:
                                    if self.is_within_intron(t) and self.SC != 'Genic-genomic':
                                        self.eval_new_SC('Genic-intron', t)
                                    else:
                                        self.eval_new_SC('Genic-genomic', t)
                    else:
                        self.eval_new_SC('Intergenic') # TODO: ref not correct it shouldnt have       
            
            if self.SC == 'GeneOverlap':
                for g_index, g in enumerate(ref.genes):
                    for t_index, trans in enumerate(g.transcripts):
                        if len(trans.SJ) == 0: # TODO: improve this classification
                            if self.strand == trans.strand:
                                if self.hit_exon(trans):
                                    self.eval_new_SC('FSM', trans)
                                else:
                                    self.eval_new_SC('Intergenic')
                            elif self.strand != trans.strand:
                                self.eval_new_SC('Antisense', trans)


    def eval_new_SC(self, new_SC: str, ref_trans = None):
        '''
        Following the SQANTi3 structural category ranking/hierarchy, evaluates
        if the SC has to be changed for a higher one in the ranking

        Args:
            new_SC (str) new structural category called for the target transcript
            ref_trans (transcript or gene deppending on the SC found) reference
        '''

        SCrank ={
            'FSM':1, 'ISM':2, 'Fusion':3,
            'NIC': 4, 'NNC':5, 'Antisense': 6,
            'Genic-genomic':7, 'Genic-intron':8, 'Intergenic':9,
            'GeneOverlap':10, 'Unclassified':11
        }
        
        if SCrank[self.SC] > SCrank[new_SC]:
            self.SC = new_SC
            if ref_trans:
                if self.SC in ['FSM', 'ISM']: # diff TSS and TTS not calculated for non-FSM/ISM
                    # TODO: see exactly how to get diff TSS and TTS
                    self.ref_trans = ref_trans.id
                    self.ref_gene = ref_trans.gene_id
                    self.diff_TSS, self.diff_TTS = self.get_diff_TSS_TTS(ref_trans)
                elif self.SC in ['NIC', 'NNC']:
                    self.ref_trans = 'novel'
                    self.ref_gene = ref_trans.gene_id
                elif self.SC == 'Intergenic':
                    self.ref_trans = 'NA'
                    self.ref_gene = 'NA'
                elif self.SC == 'Fusion':
                    self.ref_trans = 'NA'
                    self.ref_gene = '_'.join(ref_trans)
                else:
                    self.ref_trans = ref_trans.id
                    self.ref_gene = ref_trans.gene_id
                    
        
        elif SCrank[self.SC] == SCrank[new_SC]:
            if self.SC in ['FSM', 'ISM']:
                diff_TSS, diff_TTS = self.get_diff_TSS_TTS(ref_trans)
                if (self.diff_TSS + self.diff_TTS) > (diff_TSS + diff_TTS):
                    self.SC = new_SC
                    self.ref_trans = ref_trans.id
                    self.ref_gene = ref_trans.gene_id
                    self.diff_TSS = diff_TSS
                    self.diff_TTS = diff_TTS
            else:
                pass # TODO: the rest of SC


    def is_monoexon(self)-> bool:
        '''
        Check if the transcript in single-exon or multi-exon

        Returns:
            (bool) True if the transcript is mono-exon
                   False if it is multi-exon
        '''

        if len(self.SJ) == 0:
            return True
        return False

    
    def get_trans_hits(self, ref: sequence):
        '''
        Given a reference region (sequence object) find which transcripts overlap
        with the target one (comparing only TSS and TTS)

        Args:
            ref (sequence) region with its genes and transcripts associated
        '''

        for g_index, g in enumerate(ref.genes):
            for t_index, t in enumerate(g.transcripts):
                # Overlaping transcripts?
                if (self.TSS <= t.TSS and self.TTS > t.TSS) or \
                   (self.TSS < t.TTS and self.TTS >= t.TTS) or \
                   (self.TSS >= t.TSS and self.TTS <= t.TTS):

                   self.trans_hits.append((g_index, t_index))


    def is_within_intron(self, trans)-> bool:
        '''
        Check if the target transcript is completely within an intron of a
        reference transcript

        Args:
            trans (transcript) reference transcript
        
        Returns:
            (bool) True if completely within intron
                   False if not within intron
        '''

        # GTF files are 1-based and with clossed intervals, meaning start and end
        # of exons are included [start, end]
        for i in trans.SJ:
            if i[0] < self.TSS < self.TTS < i[1]:
                return True
        return False

    
    def is_within_exon(self, trans)-> bool:
        '''
        Check if the target transcript is completely within an exon of a
        reference transcript

        Args:
            trans (transcript) reference transcript
        
        Returns:
            (bool) True if completely within exon
                   False if not within exon
        '''

        coords = [st for SJ in trans.SJ for st in SJ]
        coords.insert(0, trans.TSS)
        coords.append(trans.TTS)

        exons=[]
        for i in range(0, len(coords), 2):
            exons.append((coords[i], coords[i+1]))
        
        for i in exons:
            if i[0] <= self.TSS < self.TTS <= i[1]:
                return True
        return False


    def in_intron(self, coord: int)-> bool:
        '''
        Check if a coordinate hits an intron in the target transcript

        Args:
            coord (int) genome coordinate
        
        Returns:
            (bool) True if hits intron
                   False if doesn't hit intron
        '''

        for i in self.SJ:
            if i[0] < coord < i[1]:
                return True
        return False

    
    def compare_junctions(self, trans)-> str:
        '''
        Compare the match between the junctions of the target transcript and a 
        reference one

        Args:
            trans (transcript) reference transcript

        Returns:
            (str) the match type
        '''

        ref_first_exon = (trans.TSS, trans.SJ[0][0])
        ref_last_exon = (trans.SJ[-1][-1], trans.TTS)
        self_first_exon = (self.TSS, self.SJ[0][0])
        self_last_exon = (self.SJ[-1][-1], self.TTS)
        
        if self.strand != trans.strand:
            return 'no_match'

        # Exactly the same splice junctions
        elif self.SJ == trans.SJ:
            return 'exact'

        # See if its a perfect subset
        # 1) Must be a subset of the SJ
        # 2) TSS and TTS cannot be hitting an exon
        # 3) No exon-skipping, with the requirement 1 you acomplish this one too

        #elif set(self.SJ).issubset(set(trans.SJ)) and \
        #   not trans.in_intron(self.TSS) and not trans.in_intron(self.TTS):
        #elif set(self.SJ).issubset(set(trans.SJ)) and (overlap(self_first_exon, ref_first_exon) == False or overlap(self_last_exon, ref_last_exon) == False):
        elif set(self.SJ).issubset(set(trans.SJ)):
            if self.intron_retention(trans):
                return 'partially' # TODO this is obviusly wrong
                '''
                donor = {trans.TTS}
                acceptor = {trans.TTS}

                for i in trans.SJ:
                    donor.add(i[0])
                    acceptor.add(i[1])
                                
                if self.TSS in acceptor and self.TTS in donor:
                    return 'subset'
                else:
                     return 'partially'
                '''

            elif (overlap(self_first_exon, ref_first_exon) == False or overlap(self_last_exon, ref_last_exon) == False):
                self.eval_new_SC('ISM', trans)
                            
            else:
                return 'partially'

            #return 'subset'
        
        elif self.TTS <= trans.TSS or self.TSS >= trans.TTS or self.is_within_intron(trans):
            return 'no_match'
        
        else:
            return 'partially'

    
    def genes_overlap(self, ref: sequence)-> bool:
        '''
        See if two genes overlap according to its starting and end coordinates

        Args:
            ref (sequence) region with all its associated genes
        
        Return:
            (bool) True if at least one gene overlaps
                   False no overlap
        '''

        # TODO: not using it right now
        for i in self.gene_hits:
            for j in self.gene_hits:
                if i != j:
                    if not (ref.genes[i].start >= ref.genes[j].end or ref.genes[i].end <= ref.genes[j].start):
                        return False
        return True

    
    def hits_a_gene(self, ref: sequence)-> list:
        '''
        Check if the target transcript overlaps a gene comparing start and end

        Args:
            ref (sequence) region

        Returns:
            genes (list) list of indexes of the genes that overlap
        '''

        genes = []
        for g_index, g in enumerate(ref.genes):
            if (self.TSS <= g.start <= self.TTS) or (self.TSS <= g.end <= self.TTS) or (g.start <= self.TSS < self.TTS <= g.end):
                genes.append(g_index)
        return genes

    
    def acceptor_subset(self, gene: gene)-> bool:
        '''
        Check if the acceptor splice sites of the target transcript are a subset
        of the rest of the transcripts of the gene

        Args:
            gene (gene) gene object
        
        Returns:
            (bool) True if its a subset
                   False if not
        '''

        acceptor_ref = set()
        for trans in gene.transcripts:
            for i in trans.SJ:
                acceptor_ref.add(i[1])
        
        acceptor_trans = set()
        for i in self.SJ:
            acceptor_trans.add(i[1])
        
        if acceptor_trans.issubset(acceptor_ref):
            return True
        else:
            return False

    
    def donor_subset(self, gene):
        '''
        Check if the donor splice sites of the target transcript are a subset
        of the rest of the transcripts of the gene

        Args:
            gene (gene) reference gene
        
        Returns:
            (bool) True if its a subset
                   False if not
        '''

        donor_ref = set()
        for trans in gene.transcripts:
            for i in trans.SJ:
                donor_ref.add(i[0])
        
        donor_trans = set()
        for i in self.SJ:
            donor_trans.add(i[0])
        
        if donor_trans.issubset(donor_ref):
            return True
        else:
            return False

    
    def hit_exon(self, trans)-> bool:
        '''
        Check if exons from the target transcript hits exons from the reference one

        Args:
            trans (transcript) reference transcript
        
        Returns:
            (bool) True if it hits exon
                   False if not
        '''

        coords = [st for SJ in trans.SJ for st in SJ]
        coords.insert(0, trans.TSS)
        coords.append(trans.TTS)

        exons_ref=[]
        for i in range(0, len(coords), 2):
            exons_ref.append((coords[i], coords[i+1]))
        
        coords = [st for SJ in self.SJ for st in SJ]
        coords.insert(0, self.TSS)
        coords.append(self.TTS)

        exons_self=[]
        for i in range(0, len(coords), 2):
            exons_self.append((coords[i], coords[i+1]))

        for i in exons_ref:
            for j in exons_self:
                if j[0] <= i[0] <= j[1] or j[0] <= i[1] <= j[1] or i[0] <= j[0] < j[1] <= i[1]:
                    return True
        return False


    def get_diff_TSS_TTS(self, trans)-> tuple:
        '''
        Calculates the distance from the TSS and the TTS of the target transcript
        to the closer splice site of the reference transcript. The smaller the
        distance the better match for ISM and FSM structural categories

        Args:
            trans (transcript) reference transcript

        Returns:
            (tuple) two integers with the absolute distance to the TSS and TTS
        '''

        diff_TSS = 999999999999
        diff_TTS = 999999999999 
        for i in trans.SJ:
            donor = i[0]
            acceptor = i[1]

            diff_TSS = min(diff_TSS, abs(self.TSS - acceptor))
            diff_TTS = min(diff_TTS, abs(self.TTS - donor))
        diff_TSS = min(diff_TSS, abs(self.TSS - trans.TSS))
        diff_TTS = min(diff_TTS, abs(self.TTS - trans.TTS))

        return diff_TSS, diff_TTS

    
    def intron_retention(self, trans)-> bool:
        '''
        Check if the target transcript has intron retention according to the
        reference transcript

        Args:
            trans (transcript) reference transcript
        
        Returns:
            (bool) True if there is intron retention
                   False if not
        '''

        coords = [st for SJ in self.SJ for st in SJ]
        coords.insert(0, self.TSS)
        coords.append(self.TTS)

        exons=[]
        for i in range(0, len(coords), 2):
            exons.append((coords[i], coords[i+1]))
        
        for (s,e) in exons:
            for i in trans.SJ:
                if s <= i[0] < i[1] <= e:
                    return True
        
        return False


#####################################
#                                   #
#         DEFINE FUNCTIONS          #
#                                   #
#####################################

#------------------------------------
# CLASSIFY TRANSCRIPTS IN SQANTI3 SC

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
    trans_id = None

    # Progress bar
    #num_lines = sum(1 for line in open(gtf,'r'))
    num_lines = int(subprocess.check_output('wc -l ' + gtf, shell=True).split()[0])
    # Read GTF file line by line
    with open(gtf, 'r') as f_in:
        for line in tqdm(f_in, total=num_lines):
        #for line in f_in:
            if not line.startswith('#'):
                line_split = line.split()
                feature = line_split[2]

                # Get only features that are 'exon'
                if feature == 'exon':
                    new_trans = line_split[line_split.index('transcript_id') + 1]
                    new_trans = new_trans.replace(';', '').replace('"', '')

                    if not trans_id:
                        trans_id = new_trans
                        gene_id = line_split[line_split.index('gene_id') + 1]
                        gene_id = gene_id.replace(';', '').replace('"', '')
                        seqname_id = line_split[0]
                        strand = line_split[6]
                        start = int(line_split[3])
                        end = int(line_split[4])
                        l_coords = [start, end]
                        g = gene(gene_id, seqname_id, strand, start, end) # TODO: start and end
                        res.append(sequence(seqname_id, start, end)) # TODO: start and end

                    elif new_trans == trans_id:
                        if strand == '+':
                            start = int(line_split[3])
                            end = int(line_split[4])
                        else:
                            start = int(line_split[4])
                            end = int(line_split[3])
                        l_coords.append(start)
                        l_coords.append(end)
                    
                    else:
                        if strand == '-':
                            t = transcript(trans_id, gene_id, strand, l_coords[::-1])
                        else:
                            t = transcript(trans_id, gene_id, strand, l_coords)
                        g.transcripts.append(t)
                        g.start = min(g.start, t.TSS)
                        g.end = max(g.end, t.TTS)
                        
                        new_gene = line_split[line_split.index('gene_id') + 1]
                        new_gene = new_gene.replace(';', '').replace('"', '')
                        trans_id = new_trans
                        strand = line_split[6]

                        if strand == '+':
                            start = int(line_split[3])
                            end = int(line_split[4])
                        else:
                            start = int(line_split[4])
                            end = int(line_split[3])
                        l_coords = [start, end]
                        
                        
                        if gene_id != new_gene:
                            for r_index, r in enumerate(res):
                                if r.seqname == seqname_id:
                                    if r.start <= g.start <= r.end or r.start <= g.end <= r.end:
                                        res[r_index].genes.append(g)
                                        res[r_index].start = min(r.start, g.start)
                                        res[r_index].end = max(r.end, g.end)
                                        break
                            else:
                                res.append(sequence(seqname_id, g.start, g.end))
                                res[-1].genes.append(g)
                            


                            seqname_id = line_split[0]
                            gene_id = new_gene
                            g = gene(gene_id, seqname_id, strand, min(start, end), max(start, end))
                    
    f_in.close()

    # Save last transcript 
    if strand == '-':
        t = transcript(trans_id, gene_id, strand, l_coords[::-1])
    else:
        t = transcript(trans_id, gene_id, strand, l_coords)
    g.transcripts.append(t)
    g.start = min(g.start, t.TSS)
    g.end = max(g.end, t.TTS)

    for r_index, r in enumerate(res):
        if r.seqname == seqname_id:
            if r.start <= g.start <= r.end or r.start <= g.end <= r.end:
                res[r_index].genes.append(g)
                res[r_index].start = min(r.start, g.start)
                res[r_index].end = max(r.end, g.end)
                break
    else:
        res.append(sequence(seqname_id, g.start, g.end))
        res[-1].genes.append(g)

    return res


def dont_overlap(l_genes: list)-> bool:
    '''
    Check if at least one gene of a list of genes overlap comparing start and end

    Args:
        l_genes (list) list of gene objects
    
    Returns:
        (bool) True if at least 2 genes overlap
               False if none
    '''

    for A in range(len(l_genes)):
        for B in range(len(l_genes)):
            if A != B:
                if l_genes[A].start >= l_genes[B].end or l_genes[A].end <= l_genes[B].start:
                    return True
    return False


def overlap(A: tuple, B:tuple)-> bool:
    '''
    Check if 2 intervals overlap

    Args:
        A (tuple) with start and end of the interval (int)
        B (tuple) with start and end of the interval (int)
    
    Returns:
        (bool) True if overlap
               False if not
    '''

    if A[1] >= B[0] and B[1] >= A[0]:
        return True
    return False
    

def write_SC_file(data: list, out_name: str):
    '''
    Writes the file with the structural category of each transcript and its reference

    Args:
        data (list) list of sequence objects with tanscripts classified
        out_name (str) out file name
    '''

    f_out = open(out_name, 'w')
    f_out.write('TransID\tGeneID\tSC\tRefGene\tRefTrans\n')

    for seq in data:
        for g in seq.genes:
            for t in g.transcripts:
                f_out.write(str(t.id) + '\t' + str(t.gene_id) + '\t' + str(t.SC) + '\t' + str(t.ref_gene) +'\t' + str(t.ref_trans) + '\n')
    
    f_out.close()


#------------------------------------
# MODIFY GTF FUNCTIONS

def target_trans(f_name: str, counts: dict)-> tuple:
    '''
    Choose those transcripts that will be deleted from the original GTF
    to generate the modified file to use as the reference annotation

    Args:
        f_name (str) name of the file with the GTF classification
        counts (dict) dictinary with the number of transcripts of each SC to be
                      deleted
    '''

    # TODO: use sets intead of lists and add the rest of SC
    trans4SC = {
    }

    target_trans = []
    target_genes = []
    ref_trans = []
    ref_genes = []

    # Build a list for each SC with all transcripts that were classified there
    with open(f_name, 'r') as gtf:
        header = gtf.readline()
        for line in gtf:
            line_split = line.split()
            SC = line_split[2]

            if SC in list(trans4SC.keys()):
                trans4SC[SC].append(tuple(line_split))
            else:
                trans4SC[SC] = [tuple(line_split)]
    
    gtf.close()

    # Select randomly the transcripts of each SC that are going to be deleted
    # It's important to make sure you don't delete its reference trans or gene
    for SC in list(counts.keys()):
        if counts[SC] > 0:
            SCtrans = trans4SC[SC]
            random.Random(1234).shuffle(SCtrans)
            for trans in SCtrans:
                trans_id = trans[0]
                gene_id = trans[1]
                SC = trans[2]

                if SC in ['FSM', 'ISM'] and counts[SC] > 0 and trans_id not in ref_trans and gene_id not in ref_genes:
                    ref_t = trans[4]
                    if ref_t not in target_trans:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_trans.append(ref_t)
                        counts[SC] -= 1

                elif SC in ['NIC', 'NNC'] and counts[SC] > 0 and gene_id not in ref_genes:
                    ref_g = trans[3]
                    if ref_g not in target_genes:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_genes.append(ref_g)
                        counts[SC] -= 1
                
                elif SC == 'Fusion' and counts[SC] > 0 and gene_id not in ref_genes:
                    ref_g = trans[3].split('_')
                    for i in ref_g:
                        if i in target_genes:
                            break
                    else:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_genes.extend(ref_g)
                        counts[SC] -= 1
                
                if counts[SC] <= 0:
                    break

    return target_trans, ref_genes, ref_trans


def getGeneID(line: str)-> str:
    '''
    Returns the gene_id of a GTF line

    Args:
        line (str) line readed from GTF file

    Returns:
        gene_id (str) gene_id from that feature
    '''

    line_split = line.split()
    gene_id = line_split[line_split.index('gene_id') + 1]
    gene_id = gene_id.replace(';', '').replace('"', '')
    
    return gene_id


def getTransID(line: str)-> str:
    '''
    Returns the transcript_id of a GTF line

    Args:
        line (str) line readed from GTF file

    Returns:
        trans_id (str) transcript_id from that feature
    '''

    try:
        line_split = line.split()
        trans_id = line_split[line_split.index('transcript_id') + 1]
        trans_id = trans_id.replace(';', "").replace('"', "")
    except:
        trans_id = None

    return trans_id
    

def at_least_one_trans(gene_info: list)-> bool:
    '''
    Given all the lines of a GTF associated with one gene,
    determine if after deleting the transcripts asked for is there any other
    feature associated to that gene or we need to delete it

    Args:
        gene_info (list) list of the lines of the gtf associated to that gene

    Returns:
        (bool) True there are more features associated to that gene
               False the gene is empty after the deletion of the transcripts
    '''

    for i in gene_info:
        line_split = i.split()
        if line_split[2] == 'exon':
            return True
    return False


def modifyGTF(f_name_in: str, f_name_out: str, target_trans: list, ref_genes: list):
    '''
    Modify the original GTF deleting target transcripts to simulate specific
    SQANTI3 structural categorires

    Args:
        f_name_in (str) file name of the reference annotation GTF
        f_name_out (str) file name of the modified GTF generated
        target_trans (list) list of transcripts that will be deleted
        ref_genes (list) list of genes that can't be deleted
    '''
    
    f_out = open(f_name_out, 'w')

    gene_info = []
    prev_gene = str()

    with open(f_name_in, 'r') as gtf_in:
        for line in gtf_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                gene_id = getGeneID(line)
                if not prev_gene:
                    prev_gene = gene_id
                
                if gene_id == prev_gene:
                    gene_info.append(line)

                elif prev_gene in ref_genes:
                    for i in gene_info:
                        f_out.write(i)
                    prev_gene = gene_id
                    gene_info =  [line]

                else:
                    tmp = []
                    for i in gene_info:
                        trans_id = getTransID(i)
                        if trans_id not in target_trans:
                            tmp.append(i)
                    
                    if at_least_one_trans(tmp):
                        for i in tmp:
                            f_out.write(i)
                            
                    prev_gene = gene_id
                    gene_info =  [line]
    gtf_in.close()

    if prev_gene in ref_genes:
        for i in gene_info:
            f_out.write(i)
            prev_gene = gene_id
            gene_info =  [line]

    else:
        tmp = []
        for i in gene_info:
            trans_id = getTransID(i)
            if trans_id not in target_trans:
                tmp.append(i)
                    
        if at_least_one_trans(tmp):
            for i in tmp:
                f_out.write(i)
    
    f_out.close()

    return


#####################################
#                                   #
#               MAIN                #
#                                   #
#####################################
def main():
    # Welcome
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
    
    # Arguments
    parser = argparse.ArgumentParser(prog='sqanti3_sim.py', description="SQANTI-SIM: a simulator of controlled novelty and degradation of transcripts sequence by long-reads")
    #group = parser.add_mutually_exclusive_group(required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--gtf', default = False,  help = '\t\tReference annotation in GTF format')
    group.add_argument('--cat', default = False,  help = '\t\tFile with transcripts structural categories generated with SQANTI-SIM')
    parser.add_argument('-o', '--output', default='sqanti_sim', help = '\t\tPrefix for output files')
    parser.add_argument('-d', '--dir', default='.', help = '\t\tDirectory for output files. Default: Directory where the script was run')
    parser.add_argument('--ISM', default='0', help = '\t\tNumber of incomplete-splice-matches to delete')
    parser.add_argument('--NIC', default='0', help = '\t\tNumber of novel-in-catalog to delete')
    parser.add_argument('--NNC', default='0', help = '\t\tNumber of novel-not-in-catalog to delete')
    parser.add_argument('-k', '--cores', default='1', help = '\t\tNumber of cores to run in parallel')
    parser.add_argument('-v', '--version', help='Display program version number.', action='version', version='SQANTI-SIM '+str(__version__))
    
    args = parser.parse_args()

    ref_gtf = args.gtf
    cat_in = args.cat
    out_name = args.output
    dir = args.dir

    cat_out = os.path.join(dir, (out_name + '_categories.txt'))
    gtf_modif = os.path.join(dir, (out_name + '_modified.gtf'))

    dir = '/home/jorge/Desktop/TFM'
    ref_gtf = '/home/jorge/Desktop/TFM/getSC/chr3.gencode.v38.annotation.gtf'
    #ref_gtf = '/home/jorge/Desktop/prueba.gtf'
    #ref_gtf = '/home/jorge/Desktop/simulation/ref/chr3.gencode.v38.annotation.gtf'
    out_name = 'prueba'
    cat_out = os.path.join(dir, (out_name + '_categories.txt'))
    gtf_modif = os.path.join(dir, (out_name + '_modified.gtf'))

    counts = {
        'FSM' : 0, 'ISM' : 100, 'NIC' : 0, 'NNC' : 0, 'Fusion' : 0
    }

    if ref_gtf:
        # Read GTF file
        print('Reading the GTF reference annotation file')
        data = readgtf(ref_gtf)
        print('COMPLETED\n')

        # Classify transcripts in each different sequence
        print('Classifying transcripts according to its SQANTI3 structural category')
        for seq in tqdm(range(len(data))):
            #print('\t-Classifying transcripts from sequence "%s" (%s/%s)' %(data[seq].id, seq+1, len(data)))
            data[seq].classify_trans()
            pass # TODO: run the comaprisson
        print('COMPLETED\n')

        # Write output file
        print("Writting structural category file")
        write_SC_file(data, cat_out)
        print('COMPLETED\n')

        # Show output
        print('Building summary table')
        terminal_output = summary_table()
        terminal_output.addCounts(data)
        print(terminal_output)

        cat_in = cat_out
    
    print('Writting modified GTF')
    target, ref_genes, ref_trans = target_trans(cat_in, counts)
    modifyGTF(ref_gtf, gtf_modif, target, ref_genes)
    print('COMPLETED\n')


if __name__ == '__main__':
    t_ini = time()
    main()
    t_fin = time()
    print('[Execution time %s seconds]' %(t_fin-t_ini))
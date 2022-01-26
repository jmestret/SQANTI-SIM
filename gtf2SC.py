#!/usr/bin/env python3
'''
gtf2SC.py
Classify transcripts in SQANTI3 SC if potentially deleted from GTF.
Given a GTF file as input, determine its potential SQANTI3 structural
category not taking into account himself in the reference.

Author: Jorge Mestre Tomas
Date: 19/01/2020
Last update: 25/01/2021 by Jorge Mestre
'''

__author__ = 'jormart2@alumni.uv.es'
__version__ = '0.0'

import os
import copy
import argparse
from time import time
from tqdm import tqdm


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

    def addCounts(self, data):
        '''Add counts of each SC to the table'''
        for seq in data:
            for g in seq.genes:
                for t in g.transcripts:
                    if t.SC in list(self.counts.keys()):
                        self.counts[t.SC] += 1
                    else:
                        self.counts[t.SC] = 1
    
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
        return ''


class sequence:
    def __init__(self, id, start, end):
        self.seqname = id
        self.start = start
        self.end = end
        self.genes = list() 

    def classify_trans(self):
        for g_index in range(len(self.genes)):
            for t_index in range(len(self.genes[g_index].transcripts)):
                ref = copy.deepcopy(self)
                del ref.genes[g_index].transcripts[t_index]
                if len(ref.genes[g_index].transcripts) == 0:
                    del ref.genes[g_index]
                self.genes[g_index].transcripts[t_index].get_SC(ref)


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
        self.match_type = ''
        self.ref_trans = str()
        self.ref_gene = str()
        self.diff_TSS = None
        self.diff_TTS = None
        self.gene_hits = set()
        self.trans_hits = []

    def getSJ(self, exon_coords):
        '''
        Define the splice junctions
        '''
        splice_sites = exon_coords[1:len(exon_coords)-1]
        if len(splice_sites) > 0:
            SJs = []
            for i in range(0, len(splice_sites), 2):
                SJs.append((splice_sites[i], splice_sites[i+1]))
            return(SJs)
        else:
            return([])

    def get_SC(self, ref):
        '''
        Given a group of reference trasncripts (sequence class) elucidate the SC of the
        target transcript
        '''

        # TODO: implement for sequence and gene level (Faster but no that accuarate)
        #if len(ref.genes) == 1 and len(ref.genes[0].transcripts) == 0:
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
                            self.SC = 'FSM'

                        elif self.TTS <= hit.TSS or self.TSS >= hit.TTS:
                            self.eval_new_SC('Intergenic', hit)
                        
                        elif (self.TSS <= hit.TSS and self.TTS > hit.TSS) or \
                             (self.TSS < hit.TTS and self.TTS >= hit.TTS):
                             self.eval_new_SC('Genic-genomic', hit)

            else:
                for hit_index in self.trans_hits:
                    hit = ref.genes[hit_index[0]].transcripts[hit_index[1]]
                    if self.is_within_intron(hit):
                        self.eval_new_SC('Genic-intron', hit)

                    elif self.strand != hit.strand and self.hit_exon(hit):
                        self.eval_new_SC('Antisense', hit)

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

                    if len(gene_hits) > 1 and dont_overlap(l_genes):
                        self.eval_new_SC('Fusion', ref.genes[gene_hits[0]].transcripts[0]) # TODO: which reference I include?

                    elif len(gene_hits) > 1:
                        for hit_index in self.gene_hits:
                            if self.acceptor_subset(ref.genes[hit_index]) and self.donor_subset(ref.genes[hit_index]):
                                self.eval_new_SC('NIC', ref.genes[hit_index].transcripts[0]) # TODO: which reference I include?
                            else:
                                self.eval_new_SC('NNC', ref.genes[hit_index].transcripts[0]) # TODO: which reference I include?
                    elif len(gene_hits) == 1:
                        if self.acceptor_subset(ref.genes[gene_hits[0]]) and self.donor_subset(ref.genes[gene_hits[0]]):
                            self.eval_new_SC('NIC', ref.genes[gene_hits[0]].transcripts[0]) # TODO: which reference I include?
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
                                    if self.hit_exon(t):
                                        self.eval_new_SC('Antisense', t)
                                        break
                                    else:
                                        self.eval_new_SC('Genic-genomic', t) 
                            else:
                                for t in ref.genes[index].transcripts:
                                    if self.is_within_intron(t) and self.SC != 'Genic-genomic':
                                        self.eval_new_SC('Genic-intron', t)
                                    else:
                                        self.eval_new_SC('Genic-genomic', t)
                    else:
                        self.eval_new_SC('Intergenic', ref.genes[0].transcripts[0]) # TODO: ref not correct it shouldnt have       
            #END

    def eval_new_SC(self, new_SC, ref_trans = None):
        SCrank ={
            'FSM':1, 'ISM':2, 'Fusion':3,
            'NIC': 4, 'NNC':5, 'Antisense': 6,
            'Genic-genomic':7, 'Genic-intron':8, 'Intergenic':9,
            'GeneOverlap':10, 'Unclassified':11
        }
        
        if SCrank[self.SC] > SCrank[new_SC]:
            self.SC = new_SC
            if ref_trans:
                self.ref_trans = ref_trans.id
                self.ref_gene = ref_trans.gene_id
            if self.SC in ['FSM', 'ISM']: # diff TSS and TTS not calculated for non-FSM/ISM
                self.diff_TSS = self.TSS - ref_trans.TSS
                self.diff_TTS = self.TTS -ref_trans.TTS
        
        if SCrank[self.SC] == SCrank[new_SC]:
            if self.SC in ['FSM', 'ISM'] and \
               (abs(self.diff_TSS) + abs(self.diff_TTS)) > (abs(self.TSS - ref_trans.TSS) + abs(self.TTS -ref_trans.TTS)):
                self.SC = new_SC
                if ref_trans:
                    self.ref_trans = ref_trans.id
                    self.ref_gene = ref_trans.gene_id
                self.diff_TSS = self.TSS - ref_trans.TSS
                self.diff_TTS = self.TTS -ref_trans.TTS
            else:
                pass # TODO: the rest of SC


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
        # GTF files are 1-based and with clossed intervals, meaning start and end
        # of exons are included [start, end]
        for i in trans.SJ:
            if i[0] < self.TSS < self.TTS < i[1]:
                return True
        return False
    
    def is_within_exon(self, trans):
        coords = [st for SJ in trans.SJ for st in SJ]
        coords.insert(0, trans.TSS)
        coords.append(trans.TTS)

        exons=[]
        for i in range(0, len(coords), 2):
            exons.append((coords[i], coords[i+1]))
        
        for i in exons:
            if i[0] <= trans.TSS < trans.TTS <= i[1]:
                return True
        return False

    def in_intron(self, coord):
        for i in self.SJ:
            if i[0] < coord < i[1]:
                return True
        return False
    
    def compare_junctions(self, trans):
        
        if self.strand != trans.strand:
            return 'no_match'

        # Exactly the same splice junctions
        elif self.SJ == trans.SJ:
            return 'exact'

        # See if its a perfect subset
        # 1) Must be a subset of the SJ
        # 2) TSS and TTS cannot be hitting an exon
        # 3) No exon-skipping, with the requirement 1 you acomplish this one too

        elif set(self.SJ).issubset(set(trans.SJ)) and \
           not trans.in_intron(self.TSS) and not trans.in_intron(self.TSS):
           return 'subset'
        
        elif self.TTS <= trans.TSS or self.TSS >= trans.TTS or self.is_within_intron(trans):
            return 'no_match'
        
        else:
            return 'partially'
    
    def genes_overlap(self, ref):
        # TODO: not using it right now
        for i in self.gene_hits:
            for j in self.gene_hits:
                if i != j:
                    if not (ref.genes[i].start >= ref.genes[j].end or ref.genes[i].end <= ref.genes[j].start):
                        return False
        return True
    
    def hits_a_gene(self, ref):
        genes = []
        for g_index, g in enumerate(ref.genes):
            if (self.TSS <= g.start <= self.TTS) or (self.TSS <= g.end <= self.TTS) or (g.start <= self.TSS < self.TTS <= g.end):
                genes.append(g_index)
        return genes
    
    def acceptor_subset(self, gene):
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
    
    def hit_exon(self, trans):
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
    trans_id = None

    # Progress bar
    num_lines = sum(1 for line in open(gtf,'r'))
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


def dont_overlap(l_genes):
    for A in range(len(l_genes)):
        for B in range(len(l_genes)):
            if A != B:
                if l_genes[A].start >= l_genes[B].end or l_genes[A].end <= l_genes[B].start:
                    return True
    return False

def write_output(data, out_name):
    f_out = open(out_name, 'w')
    f_out.write('trans_id\tSC\tref_trans\tref_gene')

    for seq in data:
        for g in seq.genes:
            for t in g.transcripts:
                f_out.write(str(t.id) + '\t' + str(t.SC) + '\t' + str(t.ref_trans) +'\t' + str(t.ref_gene))
                f_out.write('\n')
    
    f_out.close()


#####################################
#                                   #
#               MAIN                #
#                                   #
#####################################

def main():
    print(
        '''
        ######   #######     ###    ##    ## ######## ####          ######  #### ##     ## 
       ##    ## ##     ##   ## ##   ###   ##    ##     ##          ##    ##  ##  ###   ### 
       ##       ##     ##  ##   ##  ####  ##    ##     ##          ##        ##  #### #### 
        ######  ##     ## ##     ## ## ## ##    ##     ##  #######  ######   ##  ## ### ## 
             ## ##  ## ## ######### ##  ####    ##     ##                ##  ##  ##     ## 
       ##    ## ##    ##  ##     ## ##   ###    ##     ##          ##    ##  ##  ##     ## 
        ######   ##### ## ##     ## ##    ##    ##    ####          ######  #### ##     ## 
        '''
    )
    os.chdir('/home/jorge/Desktop')

    f_name = '/home/jorge/Desktop/simulacion/getSC/chr3.gencode.v38.annotation.gtf'
    #f_name = '/home/jorge/Desktop/prueba.gtf'
    #f_name = '/home/jorge/Desktop/simulation/ref/chr3.gencode.v38.annotation.gtf'
    out_name = 'prueba_gtf2SC.txt'

    # Read GTF file
    print('Reading the GTF reference annotation file')
    data = readgtf(f_name)
    print('COMPLETED\n')

    # Classify transcripts in each different sequence
    print('Classifying transcripts according to its SQANTI3 structural category')
    for seq in tqdm(range(len(data))):
        #print('\t-Classifying transcripts from sequence "%s" (%s/%s)' %(data[seq].id, seq+1, len(data)))
        data[seq].classify_trans()
        pass # TODO: run the comaprisson
    print('COMPLETED\n')

    # Write output file
    print("\nWritting output file")
    write_output(data, out_name)
    print('COMPLETED\n')

    # Show output
    print('Building summary table...')
    terminal_output = summary_table()
    terminal_output.addCounts(data)
    print(terminal_output)


if __name__ == '__main__':
    t_ini = time()
    main()
    t_fin = time()
    print('[Execution time %s seconds]' %(t_fin-t_ini))
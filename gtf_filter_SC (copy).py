#!/usr/bin/env python3
# Classify transcripts in possible SQANTI3 SC if deleted from GTF
# Author: Jorge Martínez Tomás


import os
import argparse
from time import time

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class summary_table:
    counts = {
            'FSM':0,
            'ISM':0,
            'NIC':0,
            'NNC':0,
            'Fusion':0,
            'Intergenic':0,
            'Genic-intron':0,
            'Antisense':0,
            'Genic-genomic':0
        }

    def addCounts(self, d_gene):
        for chr in d_gene.values():
            for gen in chr:
                for trans in gen.transcripts:
                    if trans.SC in list(self.counts.keys()):
                        self.counts[trans.SC] += 1
                    else:
                        self.counts[trans.SC] = 1
    
    def output(self):
        print(bcolors.OKBLUE + '_'*79 + bcolors.ENDC)
        print(bcolors.OKGREEN + 'S Q A N T I - S I M ' + bcolors.ENDC + '\U0001F4CA')
        print()
        print('Summary Table ' + '\U0001F50E')
        print(print(bcolors.OKBLUE + '_'*79 + bcolors.ENDC))
        for k, v in self.counts.items():
            print(bcolors.OKGREEN + '|' + bcolors.ENDC + ' ' + k + ': ' + str(v))
        

# Classes
class gene:
    def __init__(self, id, chr,strand, start, end):
        self.id = id
        self.chr = chr
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = list()  

    def trans_class(self):
        'Classify transcripts from GTF'
        # 1 trans
        if len(self.transcripts) == 1:
            self.transcripts[0].SC = 'Intergenic'
        
        # 2 trans
        elif len(self.transcripts) == 2:
            self.transcripts[0].SC = self.transcripts[0].pairwise_SC(self.transcripts[1])
            self.transcripts[0].ref = self.transcripts[1].id
            self.transcripts[1].SC = self.transcripts[1].pairwise_SC(self.transcripts[0])
            self.transcripts[1].ref = self.transcripts[0].id

            if not self.transcripts[0].SC:
                self.transcripts[0].SC = 'NNC'
            if not self.transcripts[1].SC:
                self.transcripts[1].SC = 'NNC'
    
        # 3 or more trans
        else:

            self.find_NNC()

            for i in range(len(self.transcripts)):
                for j in range(len(self.transcripts)):
                    if i != j:
                        prev_SC = self.transcripts[i].SC
                        self.transcripts[i].SC = self.transcripts[i].pairwise_SC(self.transcripts[j])
                        if self.transcripts[i].SC != prev_SC:
                            self.transcripts[i].ref = self.transcripts[j].id

                        
                            
            for i in range(len(self.transcripts)):
                if not self.transcripts[i].SC:
                    self.transcripts[i].SC = 'NIC'


    def find_NNC(self):
        d_vals = self.count_splice_sites()
        totTSS = set([trans.TSS for trans in self.transcripts])
        totTTS = set([trans.TTS for trans in self.transcripts])
        a = []
        res = {}
        for k, v in d_vals.items():
             if v == 1 and k not in totTSS and k not in totTTS:
                a.append(k)
        if a:
            for pos in a:
                for trans in self.transcripts:
                    if pos in trans.splice_sites:
                        trans.SC = 'NNC'

    
    def count_splice_sites(self):
        l_vals = [coord for trans in self.transcripts for coord in trans.splice_sites]
        d_vals = {coord : l_vals.count(coord) for coord in set(l_vals)}
        return(d_vals)

                

class transcript:
    def __init__(self, id, gene_id, strand, exon_coords):
        self.id = id
        self.gene_id = gene_id 
        self.strand = strand
        self.splice_sites = exon_coords[1:len(exon_coords)-1]
        self.TSS = exon_coords[0]
        self.TTS = exon_coords[len(exon_coords)-1]
        self.SJ = self.getSJ()
        self.donor, self.acceptor = self.donor_acceptor()
        self.exons = self.getExon(exon_coords)
        self.SC = str()
        self.ref = str()
        self.hierarchy = ["FSM", "ISM", "NIC", "NNC", "Fusion", "Intergenic", "Genic-intron", "Antisense", "Genic-genomic"]
    
    def getSJ(self):
        '''
        Define the splice junctions
        '''
        if len(self.splice_sites) > 0:
            SJs = []
            for i in range(0, len(self.splice_sites), 2):
                SJs.append((self.splice_sites[i], self.splice_sites[i+1]))
            return(SJs)
        else:
            return([])
    
    def getExon(self, coords):
        exons=[]
        for i in range(0, len(coords), 2):
            exons.append((coords[i], coords[i+1]))
        return(exons)

    def donor_acceptor(self):
        donor = []
        acceptor = []
        for i in self.SJ:
            donor.append(i[0])
            acceptor.append(i[1])
        return((donor, acceptor))

    def pairwise_SC(self, trans):
        # SC
        FSM = 'FSM'
        ISM = 'ISM'
        NIC = 'NIC'
        NNC = 'NNC'

        

        # Follow hierarchy
        #if self.SC == FSM:
        #    return(FSM)
        if 'FSM' in self.SC:
            return(self.SC)
        
        #if self.SC == ISM:
        #    if self.SJ == trans.SJ:
        #        return(FSM)
        #    else:
        #        return(ISM)
        if 'ISM' in self.SC:
            if self.SJ == trans.SJ:
                return('FSM')
            else:
                return(self.SC)

        # Check for monoexons
        if len(self.SJ) == 0:
            if len(trans.SJ) == 0:
                return(self.mono_mono(trans))
            else:
                return(self.monoexon_SC(trans))
        
        elif len(trans.SJ) == 0:
            if self.donor[0] >= trans.TSS and self.acceptor[len(self.acceptor)-1] <= trans.TTS:
                return('NNC')
            else:
                return('Genic-genomic')
        
        # Check for FSM
        if self.SJ == trans.SJ:
            return('FSM')
        
        elif set(self.SJ).issubset(set(trans.SJ)):
            checkTSS = 0
            for ex in trans.exons:
                if ex[1] == self.splice_sites[0]:
                    break
                checkTSS +=1

            checkTTS = 0
            for ex in trans.exons:
                if ex[0] == self.splice_sites[len(self.splice_sites)-1]:
                    break
                checkTTS += 1

            if checkTSS == 0 and checkTTS == len(trans.exons)-1:
                return('ISM')
            elif checkTSS == 0:
                if self.TTS <= trans.exons[checkTTS][1]:
                    return('ISM')
                else:
                    return('NIC')
            elif checkTTS == len(trans.exons)-1:
                if self.TSS >= trans.exons[checkTSS][0]:
                    return('ISM')
                else:
                    return('NIC')
            elif self.TSS >= trans.exons[checkTSS][0] and self.TTS <= trans.exons[checkTTS][1]:
                return('ISM')

            checkTSS = False
            for ex in trans.exons:
                if ex[0] <= self.TSS < ex[1]:
                    checkTSS = True
                    break

            checkTTS = False
            for ex in trans.exons:
                if ex[0] < self.TTS <= ex[1]:
                    checkTTS = True
                    break
            if checkTSS == True and checkTTS == True:
                return('ISM')
            else:
                return('NIC')
        
        elif set(self.donor).issubset(set(trans.donor)) and set(self.acceptor).issubset(set(trans.acceptor)):
            return('NIC')
        
        else:
            if not self.SC:
                return(str())
            else:
                return(self.SC)

        
    def monoexon_SC(self, trans):
        # SC
        genic_genomic = 'Genic-genomic'
        genic_intron = 'Genic-intron'
        ISM_mono = 'ISM'
        ISM_mono_intron = 'ISM'

        #if self.SC == 'FSM':
        #    return(self.SC)

        if 'FSM' in self.SC:
            return(self.SC)

        pos = trans.splice_sites.copy()
        pos.extend([trans.TSS, trans.TTS])

        # If TSS AND TTS hit
        if self.TSS in pos and self.TTS in pos:
            donorTTS = trans.donor
            donorTTS.append(trans.TTS)
            acceptorTSS = trans.acceptor
            acceptorTSS.append(trans.TSS)
            if self.TSS in trans.donor and self.TTS in donorTTS:
                return('Genic-genomic')
            elif self.TSS in trans.donor and self.TTS in trans.acceptor:
                if (self.TSS, self.TTS) in trans.SJ:
                    return('Genic-intron')
                else:
                    return('Genic-genomic')
            if self.TSS in acceptorTSS and self.TTS in donorTTS:
                if (self.TSS, self.TTS) in trans.exons:
                    return('ISM')
                else:
                    return(ISM_mono_intron)
            elif self.TSS in acceptorTSS and self.TTS in trans.acceptor:
                return(genic_genomic)
            else:
                return('Issue-with-mono-exon')

        # If ONLY TSS hit
        elif self.TSS in pos and self.TTS not in pos:
            if self.TSS == trans.TSS:
                if self.TTS < trans.donor[0]:
                    return('ISM')
                else:
                    return('Genic-genomic')
            elif self.TSS in trans.donor:
                x = trans.donor.index(self.TSS)
                if self.TTS < trans.acceptor[x]:
                    return('Genic-intron')
                else:
                    return('Genic-genomic')
            elif self.TSS in trans.acceptor:
                x = trans.acceptor.index(self.TSS)
                if x < len(trans.acceptor)-1:
                    if self.TTS < trans.donor[x+1]:
                        return('ISM')
                    else:
                        return('Genic-genomic')
                elif self.TTS < trans.TTS:
                    return('ISM')
                else:
                    return('Genic-genomic')
        
        # If ONLY TTS hit
        elif self.TSS not in pos and self.TTS in pos:
            if self.TTS == trans.TSS:
                return('Intergenic')
            elif self.TTS in trans.donor:
                x = trans.donor.index(self.TTS)
                if x > 0:
                    if self.TSS > trans.acceptor[x-1]:
                        return('ISM')
                    else:
                        return('Genic-genomic')
                elif self.TSS > trans.TSS:
                    return('ISM')
                else:
                    return('Genic-genomic')
            elif self.TTS in trans.acceptor:
                x = trans.acceptor.index(self.TTS)
                if self.TSS < trans.donor[x]:
                    return('Genic-genomic')
                else:
                    return('Genic-genomic')
            elif self.TTS == trans.TTS:
                if self.TSS >= trans.splice_sites[len(trans.splice_sites)-1]:
                    return('ISM')
                else:
                    return('Genic-genomic')
        
        # If 0 hits
        else:
            for i in trans.SJ:
                if self.TSS >= i[0] and self.TTS <= i[1]:
                    return('Genic-intron')
            
            for i in trans.exons:
                if self.TSS >= i[0] and self.TTS <= i[1]:
                    return('ISM')
            
            if self.TTS < trans.TSS or self.TSS > trans.TTS:
                return('Intergenic')
            
            else:
                return('Genic-genomic')

    def mono_mono(self, trans):
        if self.TSS >= trans.TSS and self.TTS <= trans.TTS:
            return('FSM')
        else:
            return('Intergenic') 


def single_trans_SC(d_gen):
    possible_fusion = {}
    for chr in d_gen.keys():
        possible_fusion = {}
        for i in range(len(d_gen[chr])):
            for j in range(len(d_gen[chr])):
                if i != j:
                    if len(d_gen[chr][i].transcripts) == 1:
                        #if d_gen[chr][i].transcripts[0].SC != 'Fusion':
                        if 'Fusion' not in d_gen[chr][i].transcripts[0].SC:
                            if d_gen[chr][i].strand != d_gen[chr][j].strand:
                                #if not d_gen[chr][i].transcripts[0].SC or d_gen[chr][i].transcripts[0].SC == 'Intergenic':
                                if not d_gen[chr][i].transcripts[0].SC or 'Intergenic' in d_gen[chr][i].transcripts[0].SC:
                                    if d_gen[chr][i].end > d_gen[chr][j].start:
                                        if d_gen[chr][i].start < d_gen[chr][j].start or d_gen[chr][i].start < d_gen[chr][j].end:
                                            d_gen[chr][i].transcripts[0].SC = 'Antisense'
                            else:
                                if d_gen[chr][i].start > d_gen[chr][j].start and d_gen[chr][i].start < d_gen[chr][j].end and d_gen[chr][i].end > d_gen[chr][j].end:
                                    if i in list(possible_fusion.keys()):
                                        possible_fusion[i][0] = True
                                        d_gen[chr][i].transcripts[0].SC = 'Genic-genomic'
                                    else:
                                        possible_fusion[i] = [True, False]
                                    
                                    if possible_fusion[i] == [True, True]:
                                        d_gen[chr][i].transcripts[0].SC = 'Fusion'
                                elif d_gen[chr][i].end > d_gen[chr][j].start and d_gen[chr][i].end < d_gen[chr][j].end and d_gen[chr][i].start < d_gen[chr][j].start:
                                    if i in list(possible_fusion.keys()):
                                        possible_fusion[i][1] = True
                                        d_gen[chr][i].transcripts[0].SC = 'Genic-genomic'
                                    else:
                                        possible_fusion[i] = [False, True]
                                    
                                    if possible_fusion[i] == [True, True]:
                                        d_gen[chr][i].transcripts[0].SC = 'Fusion'
                                
                                elif d_gen[chr][i].start >= d_gen[chr][j].start and d_gen[chr][i].end <= d_gen[chr][j].end:
                                    for trans in d_gen[chr][j].transcripts:
                                        d_gen[chr][i].transcripts[0].SC = d_gen[chr][i].transcripts[0].pairwise_SC(trans)
                                        if not d_gen[chr][i].transcripts[0].SC:
                                            d_gen[chr][i].transcripts[0].SC = 'NNC'


    return(d_gen)

    
def main():
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-gtf", "--gtf_file", help="Absolute path to GTF file")
    parser.add_argument("-o", "--out_dir", help="Output directory")
    parser.add_argument("-n", "--out_name", help="Output name of summary file")
    '''

    # Variables
    os.chdir('/home/jorge/Desktop')

    
    f_name = '/home/jorge/Desktop/simulacion/getSC/chr3.gencode.v38.annotation.gtf'
    
    terminal_output = summary_table()

    l_coords = list()
    res = dict()
    prev = None
    prev_trans = None

    # Antisense
    prev_chr = None
    to_check = []
    antisense = []


    with open(f_name, 'r') as f_in:
        for line in f_in:
            if line.startswith('#'):
                pass
            else:
                line_split = line.split()
                feature = line_split[2]


                if feature == 'exon':
                    trans_id = line_split[line_split.index('transcript_id') + 1]
                    trans_id = trans_id.replace(';', "")
                    trans_id = trans_id.replace('"', "")
                    if prev_trans == trans_id:
                        l_coords.append(int(line_split[3]))
                        l_coords.append(int(line_split[4]))
                    else:
                        if prev_trans:
                            #if trans_strand == '-':
                                #l_coords = l_coords[::-1]
                            t = transcript(prev_trans, gene_id, trans_strand, l_coords)
                            gen.transcripts.append(t)
                        trans_strand = line_split[6]
                        l_coords = []
                        l_coords.append(int(line_split[3]))
                        l_coords.append(int(line_split[4]))
                        prev_trans = trans_id

                    prev = feature
                
                if feature == 'gene':
                    if prev:
                        #if trans_strand == '-':
                            #l_coords = l_coords[::-1]
                        t = transcript(trans_id, gene_id, trans_strand, l_coords)
                        prev_trans = None
                        gen.transcripts.append(t)
                        gen.trans_class()
                        if gen.chr in list(res.keys()):
                            res[gen.chr].append(gen)
                        else:
                            res[gen.chr] = [gen]
                    gene_id = line_split[line_split.index('gene_id') + 1]
                    gene_id = gene_id.replace(';', '')
                    gene_id = gene_id.replace('"', '')
                    gen = gene(gene_id, line_split[0], line_split[6], int(line_split[3]), int(line_split[4]))

                    prev = feature

    
    #if trans_strand == '-':
    #    l_coords = l_coords[::-1]
    t = transcript(trans_id, gene_id, trans_strand, l_coords)
    gen.transcripts.append(t)
    gen.trans_class()
    if gen.chr in list(res.keys()):
        res[gen.chr].append(gen)
    else:
        res[gen.chr] = [gen]

    # Check for Antisense
    #print('Polishing the classification')
    res = single_trans_SC(res)
    terminal_output.addCounts(res)
            
    f_out = open('gencode.v.38.out.txt', 'w')
    f_out.write('trans_id\tSC\tref_id\tgene_id')
    for chr in res.values():
        for i in chr:
            for j in i.transcripts:
                if j.ref:
                    f_out.write(str(j.id) + '\t' + str(j.SC) + '\t' + str(j.ref) +'\t' + str(i.id))
                else:
                    f_out.write(str(j.id) + '\t' + str(j.SC) + '\t' + 'NA' +'\t' + str(i.id))
                f_out.write('\n')
    
    terminal_output.output()

if __name__ == "__main__":
    t_ini = time()
    main()
    t_fin = time()
    print(t_fin - t_ini)
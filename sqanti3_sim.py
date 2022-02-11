#!/usr/bin/env python3
'''
sqanti3_sim_sqantibased.py
Classify transcripts in SQANTI3 SC if potentially deleted from GTF.
Given a GTF file as input, determine its potential SQANTI3 structural
category not taking into account himself in the reference.
Modify original GTF deleting transcripts to simulate reads.

Author: Jorge Mestre Tomas
Modified from original SQANTI3 Quality Control Script
(https://github.com/ConesaLab/SQANTI3/blob/master/sqanti3_qc.py)
Date: 10/02/2022
'''

import os
import sys
import distutils.spawn
from collections import defaultdict, Counter, namedtuple
import bisect
import copy
import argparse
import random
import re
import subprocess
import itertools
from time import time
from tqdm import tqdm

try:
    from bx.intervals import Interval, IntervalTree
except ImportError:
    print("Unable to import bx-python! Please make sure bx-python is installed.", file=sys.stderr)
    sys.exit(-1)

try:
    from cupcake.tofu.compare_junctions import compare_junctions
    from cupcake.tofu.filter_away_subset import read_count_file
    from cupcake.io.BioReaders import GMAPSAMReader
    from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
except ImportError:
    print("Unable to import cupcake.tofu! Please make sure you install cupcake.", file=sys.stderr)
    sys.exit(-1)

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
#GTF2GENEPRED_PROG = os.path.join(utilitiesPath,"gtfToGenePred")
GTF2GENEPRED_PROG = '/home/jorge/Desktop/SQANTI3/utilities/gtfToGenePred'

if distutils.spawn.find_executable(GTF2GENEPRED_PROG) is None:
    print("Cannot find executable {0}. Abort!".format(GTF2GENEPRED_PROG), file=sys.stderr)
    sys.exit(-1)

class genePredReader(object):
    def __init__(self, filename):
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return genePredRecord.from_line(line)

class genePredRecord(object):
    def __init__(self, id, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene=None):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.txStart = txStart         # 1-based start
        self.txEnd = txEnd             # 1-based end
        self.cdsStart = cdsStart       # 1-based start
        self.cdsEnd = cdsEnd           # 1-based end
        self.exonCount = exonCount
        self.exonStarts = exonStarts   # 0-based starts # TODO Jorge: why 0-based?
        self.exonEnds = exonEnds       # 1-based ends
        self.gene = gene

        self.length = 0
        self.exons = []

        for s,e in zip(exonStarts, exonEnds):
            self.length += e-s
            self.exons.append(Interval(s, e))

        # junctions are stored (1-based last base of prev exon, 1-based first base of next exon)
        self.junctions = [(self.exonEnds[i],self.exonStarts[i+1]) for i in range(self.exonCount-1)]

    @property
    def segments(self):
        return self.exons


    @classmethod
    def from_line(cls, line):
        raw = line.strip().split('\t')
        return cls(id=raw[0],
                  chrom=raw[1],
                  strand=raw[2],
                  txStart=int(raw[3]),
                  txEnd=int(raw[4]),
                  cdsStart=int(raw[5]),
                  cdsEnd=int(raw[6]),
                  exonCount=int(raw[7]),
                  exonStarts=[int(x) for x in raw[8][:-1].split(',')],  #exonStarts string has extra , at end
                  exonEnds=[int(x) for x in raw[9][:-1].split(',')],     #exonEnds string has extra , at end
                  gene=raw[11] if len(raw)>=12 else None,
                  )

    def get_splice_site(self, genome_dict, i):
        """
        Return the donor-acceptor site (ex: GTAG) for the i-th junction
        :param i: 0-based junction index
        :param genome_dict: dict of chrom --> SeqRecord
        :return: splice site pattern, ex: "GTAG", "GCAG" etc
        """
        assert 0 <= i < self.exonCount-1

        d = self.exonEnds[i]
        a = self.exonStarts[i+1]

        seq_d = genome_dict[self.chrom].seq[d:d+2]
        seq_a = genome_dict[self.chrom].seq[a-2:a]

        if self.strand == '+':
            return (str(seq_d)+str(seq_a)).upper()
        else:
            return (str(seq_a.reverse_complement())+str(seq_d.reverse_complement())).upper()

class myQueryTranscripts:
    def __init__(self, id, gene_id, tss_diff, tts_diff, num_exons, length, str_class, subtype=None,
                 genes=None, transcripts=None, chrom=None, strand=None, bite ="NA",
                 RT_switching ="????", canonical="NA", min_cov ="NA",
                 min_cov_pos ="NA", min_samp_cov="NA", sd ="NA", FL ="NA", FL_dict={},
                 nIndels ="NA", nIndelsJunc ="NA", proteinID=None,
                 ORFlen="NA", CDS_start="NA", CDS_end="NA",
                 CDS_genomic_start="NA", CDS_genomic_end="NA", 
                 ORFseq="NA",
                 is_NMD="NA",
                 isoExp ="NA", geneExp ="NA", coding ="non_coding",
                 refLen ="NA", refExons ="NA",
                 refStart = "NA", refEnd = "NA",
                 q_splicesite_hit = 0,
                 q_exon_overlap = 0,
                 FSM_class = None, percAdownTTS = None, seqAdownTTS=None,
                 dist_cage='NA', within_cage='NA',
                 dist_polya_site='NA', within_polya_site='NA',
                 polyA_motif='NA', polyA_dist='NA', ratio_TSS='NA'):

        self.id  = id
        self.gene_id = gene_id # By Jorge
        self.tss_diff    = tss_diff   # distance to TSS of best matching ref
        self.tts_diff    = tts_diff   # distance to TTS of best matching ref
        self.tss_gene_diff = 'NA'     # min distance to TSS of all genes matching the ref
        self.tts_gene_diff = 'NA'     # min distance to TTS of all genes matching the ref
        self.genes 		 = genes if genes is not None else []
        self.AS_genes    = set()   # ref genes that are hit on the opposite strand
        self.transcripts = transcripts if transcripts is not None else []
        self.num_exons = num_exons
        self.length      = length
        self.str_class   = str_class  	# structural classification of the isoform
        self.chrom       = chrom
        self.strand 	 = strand
        self.subtype 	 = subtype
        self.RT_switching= RT_switching
        self.canonical   = canonical
        self.min_samp_cov = min_samp_cov
        self.min_cov     = min_cov
        self.min_cov_pos = min_cov_pos
        self.sd 	     = sd
        self.proteinID   = proteinID
        self.ORFlen      = ORFlen
        self.ORFseq      = ORFseq
        self.CDS_start   = CDS_start
        self.CDS_end     = CDS_end
        self.coding      = coding
        self.CDS_genomic_start = CDS_genomic_start  # 1-based genomic coordinate of CDS start - strand aware
        self.CDS_genomic_end = CDS_genomic_end      # 1-based genomic coordinate of CDS end - strand aware
        self.is_NMD      = is_NMD                   # (TRUE,FALSE) for NMD if is coding, otherwise "NA"
        self.FL          = FL                       # count for a single sample
        self.FL_dict     = FL_dict                  # dict of sample -> FL count
        self.nIndels     = nIndels
        self.nIndelsJunc = nIndelsJunc
        self.isoExp      = isoExp
        self.geneExp     = geneExp
        self.refLen      = refLen
        self.refExons    = refExons
        self.refStart    = refStart
        self.refEnd      = refEnd
        self.q_splicesite_hit = q_splicesite_hit
        self.q_exon_overlap = q_exon_overlap
        self.FSM_class   = FSM_class
        self.bite        = bite
        self.percAdownTTS = percAdownTTS
        self.seqAdownTTS  = seqAdownTTS
        self.dist_cage   = dist_cage
        self.within_cage = within_cage
        self.within_polya_site = within_polya_site
        self.dist_polya_site   = dist_polya_site    # distance to the closest polyA site (--polyA_peak, BEF file)
        self.polyA_motif = polyA_motif
        self.polyA_dist  = polyA_dist               # distance to the closest polyA motif (--polyA_motif_list, 6mer motif list)
        self.ratio_TSS = ratio_TSS

    def get_total_diff(self):
        return abs(self.tss_diff)+abs(self.tts_diff)

    def modify(self, ref_transcript, ref_gene, tss_diff, tts_diff, refLen, refExons):
        self.transcripts = [ref_transcript]
        self.genes = [ref_gene]
        self.tss_diff = tss_diff
        self.tts_diff = tts_diff
        self.refLen = refLen
        self.refExons = refExons

    def geneName(self):
        geneName = "_".join(set(self.genes))
        return geneName

    def ratioExp(self):
        if self.geneExp == 0 or self.geneExp == "NA":
            return "NA"
        else:
            ratio = float(self.isoExp)/float(self.geneExp)
        return(ratio)

    def CDSlen(self):
        if self.coding == "coding":
            return(str(int(self.CDS_end) - int(self.CDS_start) + 1))
        else:
            return("NA")

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.strand,
                                                                                                                                                           str(self.length), str(self.num_exons),
                                                                                                                                                           str(self.str_class), "_".join(set(self.genes)),
                                                                                                                                                           self.id, str(self.refLen), str(self.refExons),
                                                                                                                                                           str(self.tss_diff), str(self.tts_diff),
                                                                                                                                                           self.subtype, self.RT_switching,
                                                                                                                                                           self.canonical, str(self.min_samp_cov),
                                                                                                                                                           str(self.min_cov), str(self.min_cov_pos),
                                                                                                                                                           str(self.sd), str(self.FL), str(self.nIndels),
                                                                                                                                                           str(self.nIndelsJunc), self.bite, str(self.isoExp),
                                                                                                                                                           str(self.geneExp), str(self.ratioExp()),
                                                                                                                                                           self.FSM_class, self.coding, str(self.ORFlen),
                                                                                                                                                           str(self.CDSlen()), str(self.CDS_start), str(self.CDS_end),
                                                                                                                                                           str(self.CDS_genomic_start), str(self.CDS_genomic_end), str(self.is_NMD),
                                                                                                                                                           str(self.percAdownTTS),
                                                                                                                                                           str(self.seqAdownTTS),
                                                                                                                                                           str(self.dist_cage),
                                                                                                                                                           str(self.within_cage),
                                                                                                                                                           str(self.dist_polya_site),
                                                                                                                                                           str(self.within_polya_site),
                                                                                                                                                           str(self.polyA_motif),
                                                                                                                                                           str(self.polyA_dist), str(self.ratio_TSS))


    def as_dict(self):
        d = {'isoform': self.id,
         'chrom': self.chrom,
         'strand': self.strand,
         'length': self.length,
         'exons': self.num_exons,
         'structural_category': self.str_class,
         'associated_gene': "_".join(set(self.genes)),
         'associated_transcript': "_".join(set(self.transcripts)),
         'ref_length': self.refLen,
         'ref_exons': self.refExons,
         'diff_to_TSS': self.tss_diff,
         'diff_to_TTS': self.tts_diff,
         'diff_to_gene_TSS': self.tss_gene_diff,
         'diff_to_gene_TTS': self.tts_gene_diff,
         'subcategory': self.subtype,
         'RTS_stage': self.RT_switching,
         'all_canonical': self.canonical,
         'min_sample_cov': self.min_samp_cov,
         'min_cov': self.min_cov,
         'min_cov_pos': self.min_cov_pos,
         'sd_cov': self.sd,
         'FL': self.FL,
         'n_indels': self.nIndels,
         'n_indels_junc': self.nIndelsJunc,
         'bite': self.bite,
         'iso_exp': self.isoExp,
         'gene_exp': self.geneExp,
         'ratio_exp': self.ratioExp(),
         'FSM_class': self.FSM_class,
         'coding': self.coding,
         'ORF_length': self.ORFlen,
         'ORF_seq': self.ORFseq,
         'CDS_length': self.CDSlen(),
         'CDS_start': self.CDS_start,
         'CDS_end': self.CDS_end,
         'CDS_genomic_start': self.CDS_genomic_start,
         'CDS_genomic_end': self.CDS_genomic_end,
         'predicted_NMD': self.is_NMD,
         'perc_A_downstream_TTS': self.percAdownTTS,
         'seq_A_downstream_TTS': self.seqAdownTTS,
         'dist_to_cage_peak': self.dist_cage,
         'within_cage_peak': self.within_cage,
         'dist_to_polya_site': self.dist_polya_site,
         'within_polya_site': self.within_polya_site,
         'polyA_motif': self.polyA_motif,
         'polyA_dist': self.polyA_dist,
         'ratio_TSS' : self.ratio_TSS
         }
        for sample,count in self.FL_dict.items():
            d["FL."+sample] = count
        return d


#####################################
#                                   #
#         DEFINE FUNCTIONS          #
#                                   #
#####################################

def gtf_parser(gtf_name):
    """
    'isoform_parser()' from SQANTI3
    Parse input isoforms (GTF) to dict (chr --> sorted list)
    """
    global queryFile
    queryFile = os.path.splitext(gtf_name)[0] +".genePred"

    print("**** Parsing Isoforms....", file=sys.stderr)

    # gtf to genePred
    cmd = GTF2GENEPRED_PROG + " {0} {1} -genePredExt -allErrors -ignoreGroupsWithoutExons".format(\
        gtf_name, queryFile)
    if subprocess.check_call(cmd, shell=True)!=0:
        print("ERROR running cmd: {0}".format(cmd), file=sys.stderr)
        sys.exit(-1)


    isoforms_list = defaultdict(lambda: []) # chr --> list to be sorted later
    # will convert the sets to sorted list later
    junctions_by_chr = defaultdict(lambda: {'donors': set(), 'acceptors': set(), 'da_pairs': set()})
    # dict of gene name --> set of junctions (don't need to record chromosome)
    junctions_by_gene = defaultdict(lambda: set())
    # dict of gene name --> list of known begins and ends (begin always < end, regardless of strand)
    known_5_3_by_gene = defaultdict(lambda: {'begin':set(), 'end': set()})

    for r in genePredReader(queryFile):
        isoforms_list[r.chrom].append(r)
        known_5_3_by_gene[r.gene]['begin'].add(r.txStart)
        known_5_3_by_gene[r.gene]['end'].add(r.txEnd)

        if r.exonCount >= 2:
            for d, a in r.junctions:
                junctions_by_chr[r.chrom]['donors'].add(d)
                junctions_by_chr[r.chrom]['acceptors'].add(a)
                junctions_by_chr[r.chrom]['da_pairs'].add((d,a))
                junctions_by_gene[r.gene].add((d,a))

    for k in junctions_by_chr:
        junctions_by_chr[k]['donors'] = list(junctions_by_chr[k]['donors'])
        junctions_by_chr[k]['donors'].sort()
        junctions_by_chr[k]['acceptors'] = list(junctions_by_chr[k]['acceptors'])
        junctions_by_chr[k]['acceptors'].sort()
        junctions_by_chr[k]['da_pairs'] = list(junctions_by_chr[k]['da_pairs'])
        junctions_by_chr[k]['da_pairs'].sort()

    for k in isoforms_list:
        isoforms_list[k].sort(key=lambda r: r.txStart)

    return isoforms_list, dict(junctions_by_chr), dict(junctions_by_gene), dict(known_5_3_by_gene)


def transcript_classification(trans_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene):
    res = defaultdict(lambda: [])
    for chrom, records in trans_by_chr.items(): #TODO Do by region instead of chromosome to make easier searchs
        for trans in tqdm(range(len(records))):
            trans = records[trans]

            # Esto es extremadamente lento!
            # will convert the sets to sorted list later
            junctions_by_chr = defaultdict(lambda: {'donors': set(), 'acceptors': set(), 'da_pairs': set()})
            # dict of gene name --> set of junctions (don't need to record chromosome)
            junctions_by_gene = defaultdict(lambda: set())
            # dict of gene name --> list of known begins and ends (begin always < end, regardless of strand)
            known_5_3_by_gene = defaultdict(lambda: {'begin':set(), 'end': set()})

            for r in records:
                if trans.id == r.id:
                    continue
                known_5_3_by_gene[r.gene]['begin'].add(r.txStart)
                known_5_3_by_gene[r.gene]['end'].add(r.txEnd)
                if r.exonCount >= 2:
                    for d, a in r.junctions:
                        junctions_by_chr[r.chrom]['donors'].add(d)
                        junctions_by_chr[r.chrom]['acceptors'].add(a)
                        junctions_by_chr[r.chrom]['da_pairs'].add((d,a))
                        junctions_by_gene[r.gene].add((d,a))
            
            for k in junctions_by_chr:
                junctions_by_chr[k]['donors'] = list(junctions_by_chr[k]['donors'])
                junctions_by_chr[k]['donors'].sort()
                junctions_by_chr[k]['acceptors'] = list(junctions_by_chr[k]['acceptors'])
                junctions_by_chr[k]['acceptors'].sort()
                junctions_by_chr[k]['da_pairs'] = list(junctions_by_chr[k]['da_pairs'])
                junctions_by_chr[k]['da_pairs'].sort()

            isoform_hit = transcriptsKnownSpliceSites(trans, records, start_ends_by_gene)

            if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
                # not FSM or ISM --> see if it is NIC, NNC, or fusion
                isoform_hit = novelIsoformsKnownGenes(isoform_hit, trans, junctions_by_chr, junctions_by_gene, start_ends_by_gene)
            elif isoform_hit.str_class in ("", "geneOverlap"):
                # possibly NNC, genic, genic intron, anti-sense, or intergenic
                isoform_hit = associationOverlapping(isoform_hit, trans, junctions_by_chr)

            # Save trans classification
            res[isoform_hit.chrom].append(isoform_hit)
    
    return res

def transcriptsKnownSpliceSites(trec, ref_chr, start_ends_by_gene):
    def calc_overlap(s1, e1, s2, e2):
        if s1=='NA' or s2=='NA': return 0
        if s1 > s2:
            s1, e1, s2, e2 = s2, e2, s1, e1
        return max(0, min(e1,e2)-max(s1,s2))

    def gene_overlap(ref1, ref2):
        if ref1==ref2: return True  # same gene, diff isoforms
        # return True if the two reference genes overlap
        s1, e1 = min(start_ends_by_gene[ref1]['begin']), max(start_ends_by_gene[ref1]['end'])
        s2, e2 = min(start_ends_by_gene[ref2]['begin']), max(start_ends_by_gene[ref2]['end'])
        if s1 <= s2:
            return e1 <= s2
        else:
            return e2 <= s1

    def calc_splicesite_agreement(query_exons, ref_exons):
        q_sites = {}
        for e in query_exons:
            q_sites[e.start] = 0
            q_sites[e.end] = 0
        for e in ref_exons:
            if e.start in q_sites: q_sites[e.start] = 1
            if e.end in q_sites: q_sites[e.end] = 1
        return sum(q_sites.values())

    def calc_exon_overlap(query_exons, ref_exons):
        q_bases = {}
        for e in query_exons:
            for b in range(e.start, e.end): q_bases[b] = 0

        for e in ref_exons:
            for b in range(e.start, e.end):
                if b in q_bases: q_bases[b] = 1
        return sum(q_bases.values())

    def get_diff_tss_tts(trec, ref):
        if trec.strand == '+':
            diff_tss = trec.txStart - ref.txStart
            diff_tts = ref.txEnd - trec.txEnd
        else:
            diff_tts = trec.txStart - ref.txStart
            diff_tss = ref.txEnd - trec.txEnd
        return diff_tss, diff_tts


    def get_gene_diff_tss_tts(isoform_hit):
        # now that we know the reference (isoform) it hits
        # add the nearest start/end site for that gene (all isoforms of the gene)
        nearest_start_diff, nearest_end_diff = float('inf'), float('inf')
        for ref_gene in isoform_hit.genes:
            for x in start_ends_by_gene[ref_gene]['begin']:
                d = trec.txStart - x
                if abs(d) < abs(nearest_start_diff):
                    nearest_start_diff = d
            for x in start_ends_by_gene[ref_gene]['end']:
                d = trec.txEnd - x
                if abs(d) < abs(nearest_end_diff):
                    nearest_end_diff = d

        if trec.strand == '+':
            isoform_hit.tss_gene_diff = nearest_start_diff if nearest_start_diff!=float('inf') else 'NA'
            isoform_hit.tts_gene_diff = nearest_end_diff if nearest_end_diff!=float('inf') else 'NA'
        else:
            isoform_hit.tss_gene_diff = -nearest_end_diff if nearest_start_diff!=float('inf') else 'NA'
            isoform_hit.tts_gene_diff = -nearest_start_diff if nearest_end_diff!=float('inf') else 'NA'

    def categorize_incomplete_matches(trec, ref):
        """
        intron_retention --- at least one trec exon covers at least two adjacent ref exons
        complete --- all junctions agree and is not IR
        5prime_fragment --- all junctions agree but trec has less 5' exons. The isoform is a 5' fragment of the reference transcript
        3prime_fragment --- all junctions agree but trec has less 3' exons. The isoform is a 3' fragment of the reference transcript
        internal_fragment --- all junctions agree but trec has less 5' and 3' exons
        """
        # check intron retention
        ref_exon_tree = IntervalTree()
        for i,e in enumerate(ref.exons): ref_exon_tree.insert(e.start, e.end, i)
        for e in trec.exons:
            if len(ref_exon_tree.find(e.start, e.end)) > 1: # multiple ref exons covered
                return "intron_retention"

        agree_front = trec.junctions[0]==ref.junctions[0]
        agree_end   = trec.junctions[-1]==ref.junctions[-1]
        if agree_front:
            if agree_end:
                return "complete"
            else: # front agrees, end does not
                return ("5prime_fragment" if trec.strand=='+' else '3prime_fragment')
        else:
            if agree_end: # front does not agree, end agrees
                return ("3prime_fragment" if trec.strand=='+' else '5prime_fragment')
            else:
                return "internal_fragment"

    isoform_hit = myQueryTranscripts(id=trec.id, gene_id=trec.gene, tts_diff="NA", tss_diff="NA",\
                                    num_exons=trec.exonCount,
                                    length=trec.length,
                                    str_class="",
                                    chrom=trec.chrom,
                                    strand=trec.strand,
                                    subtype="no_subcategory")
    
    cat_ranking = {'full-splice_match': 5, 'incomplete-splice_match': 4, 'anyKnownJunction': 3, 'anyKnownSpliceSite': 2,
                   'geneOverlap': 1, '': 0}

    #----------------------------------#
    #       SPLICED TRANSCRIPT         #
    #----------------------------------#
    if trec.exonCount >= 2:
        hits_by_gene = defaultdict(lambda: [])  # gene --> list of hits
        best_by_gene = {}  # gene --> best isoform_hit

        for ref in ref_chr:
            if trec.id != ref.id: # to not match with itself
                if hits_exon(trec, ref):
                    hits_by_gene[ref.gene].append(ref)
        
        if len(hits_by_gene) == 0: return isoform_hit

        for ref_gene in hits_by_gene:
            isoform_hit = myQueryTranscripts(id=trec.id, gene_id=trec.gene, tts_diff="NA", tss_diff="NA",\
                                            num_exons=trec.exonCount,
                                            length=trec.length,
                                            str_class="",
                                            chrom=trec.chrom,
                                            strand=trec.strand,
                                            subtype="no_subcategory")

            for ref in hits_by_gene[ref_gene]:
                if trec.strand != ref.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
            
                #--MONO-EXONIC REFERENCE--#
                if ref.exonCount == 1:
                    if ref.exonCount == 1: # mono-exonic reference, handle specially here
                        if calc_exon_overlap(trec.exons, ref.exons) > 0 and cat_ranking[isoform_hit.str_class] < cat_ranking["geneOverlap"]:
                            isoform_hit = myQueryTranscripts(trec.id, trec.gene, "NA", "NA", trec.exonCount, trec.length,
                                                                "geneOverlap",
                                                                subtype="mono-exon",
                                                                chrom=trec.chrom,
                                                                strand=trec.strand,
                                                                genes=[ref.gene],
                                                                transcripts=[ref.id],
                                                                refLen=ref.length,
                                                                refExons=ref.exonCount,
                                                                refStart=ref.txStart,
                                                                refEnd=ref.txEnd,
                                                                q_splicesite_hit=0,
                                                                q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons))
                #--MULTI-EXONIC REFERENCE--#
                else:
                    match_type = compare_junctions(trec, ref, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999)
                    
                    if match_type not in ('exact', 'subset', 'partial', 'concordant', 'super', 'nomatch'):
                        raise Exception("Unknown match category {0}!".format(match_type))

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)


                    # #############################
                    # SQANTI's full-splice_match
                    # #############################
                    if match_type == "exact":
                        subtype = "multi-exon"
                        # assign as a new hit if
                        # (1) no prev hits yet
                        # (2) this one is better (prev not FSM or is FSM but worse tss/tts)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["full-splice_match"] or \
                                                    abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                            # subcategory for matching 5' and matching 3'
                            if abs(diff_tss) <= 50 and abs(diff_tts) <= 50:
                                    subtype = 'reference_match'
                            # subcategory for matching 5' and non-matching 3'
                            if abs(diff_tss) <= 50 and abs(diff_tts) > 50:
                                subtype = 'alternative_3end'
                            # subcategory for matching 3' and non-matching 5'
                            if abs(diff_tss) > 50 and abs(diff_tts) <= 50:
                                subtype = 'alternative_5end'
                            # subcategory for non-matching 3' and non-matching 5'
                            if abs(diff_tss) > 50 and abs(diff_tts) > 50:
                                subtype = 'alternative_3end5end'
                            isoform_hit = myQueryTranscripts(trec.id, trec.gene, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                              str_class="full-splice_match",
                                                              subtype=subtype,
                                                              chrom=trec.chrom,
                                                              strand=trec.strand,
                                                              genes=[ref.gene],
                                                              transcripts=[ref.id],
                                                              refLen = ref.length,
                                                              refExons= ref.exonCount,
                                                              refStart=ref.txStart,
                                                              refEnd=ref.txEnd,
                                                              q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                              q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons))
                    
                    # #######################################################
                    # SQANTI's incomplete-splice_match
                    # (only check if don't already have a FSM match)
                    # #######################################################
                    elif match_type == "subset":
                        subtype = categorize_incomplete_matches(trec, ref)
                        # assign as a new (ISM) hit if
                        # (1) no prev hit
                        # (2) prev hit not as good (is ISM with worse tss/tts or anyKnownSpliceSite)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["incomplete-splice_match"] or \
                            (isoform_hit.str_class=='incomplete-splice_match' and abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff()):
                            isoform_hit = myQueryTranscripts(trec.id, trec.gene, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                             str_class="incomplete-splice_match",
                                                             subtype=subtype,
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=[ref.id],
                                                             refLen = ref.length,
                                                             refExons= ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                             q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons))

                    
                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ('partial', 'concordant', 'super'):
                        q_sp_hit = calc_splicesite_agreement(trec.exons, ref.exons)
                        q_ex_overlap = calc_exon_overlap(trec.exons, ref.exons)
                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["anyKnownJunction"] or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit > isoform_hit.q_splicesite_hit) or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit==isoform_hit.q_splicesite_hit and q_ex_overlap > isoform_hit.q_exon_overlap) or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit==isoform_hit.q_splicesite_hit and q_exon_d < abs(trec.exonCount-isoform_hit.refExons)):
                            isoform_hit = myQueryTranscripts(trec.id, trec.gene, "NA", "NA", trec.exonCount, trec.length,
                                                             str_class="anyKnownJunction",
                                                             subtype="no_subcategory",
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=["novel"],
                                                             refLen=ref.length,
                                                             refExons=ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                             q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons))
                    
                    else: # must be nomatch
                        assert match_type == 'nomatch'
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["anyKnownSpliceSite"] and calc_splicesite_agreement(trec.exons, ref.exons) > 0:
                            isoform_hit = myQueryTranscripts(trec.id, trec.gene, "NA", "NA", trec.exonCount, trec.length,
                                                             str_class="anyKnownSpliceSite",
                                                             subtype="no_subcategory",
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=["novel"],
                                                             refLen=ref.length,
                                                             refExons=ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                             q_exon_overlap=calc_exon_overlap(trec.exons,
                                                                                              ref.exons))

                        if isoform_hit.str_class=="": # still not hit yet, check exonic overlap
                            if cat_ranking[isoform_hit.str_class] < cat_ranking["geneOverlap"] and calc_exon_overlap(trec.exons, ref.exons) > 0:
                                isoform_hit = myQueryTranscripts(trec.id, trec.gene, "NA", "NA", trec.exonCount, trec.length,
                                                                 str_class="geneOverlap",
                                                                 subtype="no_subcategory",
                                                                 chrom=trec.chrom,
                                                                 strand=trec.strand,
                                                                 genes=[ref.gene],
                                                                 transcripts=["novel"],
                                                                 refLen=ref.length,
                                                                 refExons=ref.exonCount,
                                                                 refStart=ref.txStart,
                                                                 refEnd=ref.txEnd,
                                                                 q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                                 q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons))
            best_by_gene[ref_gene] = isoform_hit

        # now we have best_by_gene:
        # start with the best scoring one (FSM is best) --> can add other genes if they don't overlap
        #if trec.id.startswith('PB.1252.'):
        #    pdb.set_trace()
        geneHitTuple = namedtuple('geneHitTuple', ['score', 'rStart', 'rEnd', 'rGene', 'iso_hit'])
        best_by_gene = [geneHitTuple(cat_ranking[iso_hit.str_class],iso_hit.refStart,iso_hit.refEnd,ref_gene,iso_hit) for ref_gene,iso_hit in best_by_gene.items()]
        best_by_gene = list(filter(lambda x: x.score > 0, best_by_gene))
        if len(best_by_gene) == 0: # no hit
            return isoform_hit

        
        best_by_gene.sort(key=lambda x: (x.score,x.iso_hit.q_splicesite_hit+(x.iso_hit.q_exon_overlap)*1./sum(e.end-e.start for e in trec.exons)+calc_overlap(x.rStart,x.rEnd,trec.txStart,trec.txEnd)*1./(x.rEnd-x.rStart)-abs(trec.exonCount-x.iso_hit.refExons)), reverse=True)  # sort by (ranking score, overlap)
        isoform_hit = best_by_gene[0].iso_hit
        cur_start, cur_end = best_by_gene[0].rStart, best_by_gene[0].rEnd
        for t in best_by_gene[1:]:
            if t.score==0: break
            if calc_overlap(cur_start, cur_end, t.rStart, t.rEnd) <= 0:
                isoform_hit.genes.append(t.rGene)
                cur_start, cur_end = min(cur_start, t.rStart), max(cur_end, t.rEnd)

    #----------------------------------#
    #       UNSPLICED TRANSCRIPT       #
    #----------------------------------#
    else:
        for ref in ref_chr:
            if trec.id != ref.id: # to not match with itself
                if hits_exon(trec, ref) and ref.exonCount == 1:
                    if ref.strand != trec.strand:
                        # opposite strand, just record it in AS_genes
                        isoform_hit.AS_genes.add(ref.gene)
                        continue
                
                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                    if isoform_hit.str_class == "": # no match so far
                        isoform_hit = myQueryTranscripts(trec.id, trec.gene, diff_tss, diff_tts, trec.exonCount, trec.length, "full-splice_match",
                                                                subtype="mono-exon",
                                                                chrom=trec.chrom,
                                                                strand=trec.strand,
                                                                genes=[ref.gene],
                                                                transcripts=[ref.id],
                                                                refLen=ref.length,
                                                                refExons = ref.exonCount)
                    elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                        isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)
        
        if isoform_hit.str_class == "":
            for ref in ref_chr:
                if trec.id != ref.id: # to not match with itself
                    if hits_exon(trec, ref) and ref.exonCount >= 2:
                        if calc_exon_overlap(trec.exons, ref.exons) == 0:   # no exonic overlap, skip!
                            continue

                        if ref.strand != trec.strand:
                            # opposite strand, just record it in AS_genes
                            isoform_hit.AS_genes.add(ref.gene)
                            continue

                        diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                        for e in ref.exons:
                            if e.start <= trec.txStart < trec.txEnd <= e.end:
                                isoform_hit.str_class = "incomplete-splice_match"
                                isoform_hit.subtype = "mono-exon"
                                isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)
                                # this is as good a match as it gets, we can stop the search here
                                get_gene_diff_tss_tts(isoform_hit)
                                return isoform_hit

                        # if we haven't exited here, then ISM hit is not found yet
                        # instead check if it's NIC by intron retention
                        # but we don't exit here since the next gene could be a ISM hit
                        if ref.txStart <= trec.txStart < trec.txEnd <= ref.txEnd:
                            isoform_hit.str_class = "novel_in_catalog"
                            isoform_hit.subtype = "mono-exon"
                            # check for intron retention
                            if len(ref.junctions) > 0:
                                for (d,a) in ref.junctions:
                                    if trec.txStart < d < a < trec.txEnd:
                                        isoform_hit.subtype = "mono-exon_by_intron_retention"
                                        break
                            isoform_hit.modify("novel", ref.gene, 'NA', 'NA', ref.length, ref.exonCount)
                            get_gene_diff_tss_tts(isoform_hit)
                            return isoform_hit

                        # if we get to here, means neither ISM nor NIC, so just add a ref gene and categorize further later
                        isoform_hit.genes.append(ref.gene)

    get_gene_diff_tss_tts(isoform_hit)
    isoform_hit.genes.sort(key=lambda x: start_ends_by_gene[x]['begin'])
    return isoform_hit



def novelIsoformsKnownGenes(isoforms_hit, trec, junctions_by_chr, junctions_by_gene, start_ends_by_gene):
    """
    At this point: definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    :return isoforms_hit: updated isoforms hit (myQueryTranscripts object)
    """
    def has_intron_retention():
        for e in trec.exons:
            m = bisect.bisect_left(junctions_by_chr[trec.chrom]['da_pairs'], (e.start, e.end))
            if m < len(junctions_by_chr[trec.chrom]['da_pairs']) and e.start <= junctions_by_chr[trec.chrom]['da_pairs'][m][0] < junctions_by_chr[trec.chrom]['da_pairs'][m][1] < e.end:
                return True
        return False

    ref_genes = list(set(isoforms_hit.genes))

    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update refLen or refExon
    # or tss_diff and tts_diff (always set to "NA" for non-FSM/ISM matches)
    #
    isoforms_hit.transcripts = ["novel"]
    if len(ref_genes) == 1:
        # hits exactly one gene, must be either NIC or NNC
        ref_gene_junctions = junctions_by_gene[ref_genes[0]]
        # 1. check if all donors/acceptor sites are known (regardless of which ref gene it came from)
        # 2. check if this query isoform uses a subset of the junctions from the single ref hit
        all_junctions_known = True
        all_junctions_in_hit_ref = True
        for d,a in trec.junctions:
            all_junctions_known = all_junctions_known and (d in junctions_by_chr[trec.chrom]['donors']) and (a in junctions_by_chr[trec.chrom]['acceptors'])
            all_junctions_in_hit_ref = all_junctions_in_hit_ref and ((d,a) in ref_gene_junctions)
        if all_junctions_known:
            isoforms_hit.str_class="novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoforms_hit.subtype = "combination_of_known_junctions"
            else:
                isoforms_hit.subtype = "combination_of_known_splicesites"
        else:
            isoforms_hit.str_class="novel_not_in_catalog"
            isoforms_hit.subtype = "at_least_one_novel_splicesite"
    else: # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        # NOTE: some ref genes could be mono-exonic so no junctions
        all_ref_junctions = list(itertools.chain(junctions_by_gene[ref_gene] for ref_gene in ref_genes if ref_gene in junctions_by_gene))

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = dict((i, all_ref_junctions.count(junc)) for i,junc in enumerate(trec.junctions))

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.values()) > 1:
            isoforms_hit.str_class = "moreJunctions"
        else:
            isoforms_hit.str_class = "fusion"
            isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    if has_intron_retention():
        isoforms_hit.subtype = "intron_retention"

    return isoforms_hit


def associationOverlapping(isoforms_hit, trec, junctions_by_chr):
    # at this point: definitely not FSM or ISM or NIC or NNC
    # possibly (in order of preference assignment):
    #  - antisense  (on opp strand of a known gene)
    #  - genic (overlaps a combination of exons and introns), ignore strand
    #  - genic_intron  (completely within an intron), ignore strand
    #  - intergenic (does not overlap a gene at all), ignore strand

    isoforms_hit.str_class = "intergenic"
    isoforms_hit.transcripts = ["novel"]
    isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    if len(isoforms_hit.genes) == 0:
        # completely no overlap with any genes on the same strand
        # check if it is anti-sense to a known gene, otherwise it's genic_intron or intergenic
        if len(isoforms_hit.AS_genes) == 0:
            if trec.chrom in junctions_by_chr:
                # no hit even on opp strand
                # see if it is completely contained within a junction
                da_pairs = junctions_by_chr[trec.chrom]['da_pairs']
                i = bisect.bisect_left(da_pairs, (trec.txStart, trec.txEnd))
                while i < len(da_pairs) and da_pairs[i][0] <= trec.txStart:
                    if da_pairs[i][0] <= trec.txStart <= trec.txStart <= da_pairs[i][1]:
                        isoforms_hit.str_class = "genic_intron"
                        break
                    i += 1
            else:
                pass # remain intergenic
        else:
            # hits one or more genes on the opposite strand
            isoforms_hit.str_class = "antisense"
            isoforms_hit.genes = ["novelGene_{g}_AS".format(g=g) for g in isoforms_hit.AS_genes]
    else:
        # (Liz) used to put NNC here - now just genic
        isoforms_hit.str_class = "genic"
        # overlaps with one or more genes on the same strand
        #if trec.exonCount >= 2:
        #    # multi-exon and has a same strand gene hit, must be NNC
        #    isoforms_hit.str_class = "novel_not_in_catalog"
        #    isoforms_hit.subtype = "at_least_one_novel_splicesite"
        #else:
        #    # single exon, must be genic
        #    isoforms_hit.str_class = "genic"

    return isoforms_hit

def hits_exon(r1, r2):
    '''
    Check if any exon of r2 is overlaped by r1 start and end
    IntervalTree.find(end, start): Return a sorted list of all intervals overlapping [start,end)
    Overlap, not within that's the point
    '''
        
    for e in r2.exons:
        if r1.txStart <= e.end and e.start <= r1.txEnd:
            return True
    return False


def summary_table(data: dict):
    counts = defaultdict(lambda: 0, {
        'full-splice_match': 0,
        'incomplete-splice_match':0,
        'novel_in_catalog':0,
        'novel_not_in_catalog':0,
        'fusion' : 0,
        'antisense': 0,
        'genic_intron': 0,
        'genic' :0,
        'intergenic':0
    })
    for chrom in data.values():
        for trans in chrom:
            counts[trans.str_class] += 1
    
    print('\033[94m_' * 79 + '\033[0m')
    print('\033[92mS Q A N T I - S I M\033[0m \U0001F4CA')
    print()
    print('Summary Table \U0001F50E')
    print('\033[94m_' * 79 + '\033[0m')
    for k, v in counts.items():
        print('\033[92m|\033[0m ' + k + ': ' + str(v))


def write_category_file(data: dict, out_name: str):
    '''
    Writes the file with the structural category of each transcript and its reference

    Args:
        data (dict) list of sequence objects with tanscripts classified
        out_name (str) out file name
    '''

    f_out = open(out_name, 'w')
    f_out.write('TransID\tGeneID\tSC\tRefGene\tRefTrans\n')

    for chrom in data.values():
        for trans in chrom:
            f_out.write('%s\t%s\t%s\t%s\t%s\n' %(trans.id, trans.gene_id, trans.str_class, '_'.join(trans.genes), '_'.join(trans.transcripts)))

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

                if SC in ['full-splice_match', 'incomplete-splice_match', 'antisense', 'genic', 'genic_intron'] and counts[SC] > 0 and trans_id not in ref_trans and gene_id not in ref_genes:
                    ref_t = trans[4]
                    if ref_t not in target_trans:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_trans.append(ref_t)
                        counts[SC] -= 1

                elif SC in ['novel_in_catalog', 'novel_not_in_catalog'] and counts[SC] > 0 and gene_id not in ref_genes:
                    ref_g = trans[3]
                    if ref_g not in target_genes:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_genes.append(ref_g)
                        counts[SC] -= 1
                
                elif SC == 'fusion' and counts[SC] > 0 and gene_id not in ref_genes:
                    ref_g = trans[3].split('_')
                    for i in ref_g:
                        if i in target_genes:
                            break
                    else:
                        target_trans.append(trans_id)
                        target_genes.append(gene_id)
                        ref_genes.extend(ref_g)
                        counts[SC] -= 1
                
                elif SC == 'intergenic' and counts[SC] > 0 and gene_id not in ref_genes:
                    target_trans.append(trans_id)
                    target_genes.append(gene_id)
                
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
                
                    '''
                    elif prev_gene in ref_genes:
                        for i in gene_info:
                            f_out.write(i)
                        prev_gene = gene_id
                        gene_info =  [line]
                    '''

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

    # Input data
    dir = '/home/jorge/Desktop/ConesaLab/SQANTI-SIM'
    ref_gtf = '/home/jorge/Desktop/simulation/ref/chr3.gencode.v38.annotation.gtf'
    out_name = 'prueba'
    cat_out = os.path.join(dir, (out_name + '_categories.txt'))
    gtf_modif = os.path.join(dir, (out_name + '_modified.gtf'))


    # Parsing transcripts from GTF
    trans_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene = gtf_parser(ref_gtf)

    

    # Classify transcripts
    trans_info = transcript_classification(trans_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene)

    # Print summary table
    summary_table(trans_info)

    # Write category file
    write_category_file(trans_info, cat_out)

    # Write modified GTF
    cat_in = cat_out
    counts = defaultdict(lambda: 0, {
        'full-splice_match': 0,
        'incomplete-splice_match':0,
        'novel_in_catalog':0,
        'novel_not_in_catalog':0,
        'fusion' : 100,
        'antisense': 0,
        'genic_intron': 0,
        'genic' :0,
        'intergenic':0
    })
    print('Writting modified GTF')
    target, ref_genes, ref_trans = target_trans(cat_in, counts)
    modifyGTF(ref_gtf, gtf_modif, target, ref_genes)
    print('COMPLETED\n')



if __name__ == '__main__':
    t_ini = time()
    main()
    t_fin = time()
    print('[Execution time %s seconds]' %(t_fin-t_ini))







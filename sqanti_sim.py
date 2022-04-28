#!/usr/bin/env python3
"""
sqanti3_sim.py

Wrapper for long-read RNA-seq simulators (NanoSim and IsoSeqSim) to simulate
controlled novelty and degradation of transcripts based on SQANTI3 structural
categories

@author Jorge Mestre Tomas (jormart2@alumni.uv.es)
@date 19/01/2022
"""

__version__ = "1.0b1"

import argparse
import numpy
import os
import random
import sys
from collections import defaultdict
from src import classif_gtf
from src import pb_ont_sim
from src import sim_preparatory
from src import sqanti3_stats
from time import strftime


def classif(input: list):
    """Classify transcripts in SQANTI3 structural categories

    Given a GTF annotation generates an index file with the potential SQANTI3
    structural category of each transcript

    Args:
        input (list): arguments to parse
    """

    parser = argparse.ArgumentParser( prog="sqanti_sim.py classif", description="sqanti_sim.py classif parse options", )
    parser.add_argument( "--gtf", type=str, help="\t\tReference annotation in GTF format", required=True, )
    parser.add_argument( "-o", "--output", default="sqanti_sim", help="\t\tPrefix for output file", )
    parser.add_argument( "-d", "--dir", default=".", help="\t\tDirectory for output files (default: .)", )
    parser.add_argument( "-k", "--cores", default=1, type=int, help="\t\tNumber of cores to run in parallel", )

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print(
            "[SQANTI-SIM] classif mode unrecognized arguments: {}\n".format(
                " ".join(unknown)
            ),
            file=sys.stderr,
        )

    print("\n[SQANTI-SIM] Running with the following parameters:")
    print("[SQANTI-SIM] - Ref GTF:", str(args.gtf))
    print("[SQANTI-SIM] - Out prefix:", str(args.output))
    print("[SQANTI-SIM] - Out dir:", str(args.dir))
    print("[SQANTI-SIM] - N threads:", str(args.cores))

    print("\n[SQANTI-SIM][%s] Classifying transcripts in structural categories" %(strftime("%d-%m-%Y %H:%M:%S")))
    trans_info = classif_gtf.classify_gtf(args)
    
    print("[SQANTI-SIM] Summary table from categorization")
    classif_gtf.summary_table_cat(trans_info)

    print("[SQANTI-SIM][%s] classif step finished" %(strftime("%d-%m-%Y %H:%M:%S")))


def preparatory(input: list):
    """Modifies reference annotation GTF and builds expression matrix

    Given the novel and known transcripts to simulate and its counts, generates
    the expression matrix to give as input to the long-read RNA-seq simulators
    and generetes the modified GTF to use as reference annotation in your tool

    Args:
        input (list): arguments to parse
    """
    parser = argparse.ArgumentParser(prog="sqanti_sim.py preparatory", description="sqanti_sim.py preparatory parse options", )
    subparsers = parser.add_subparsers(dest="mode", description="\t\tDifferent modes to generate the expression matrix: equal (simulate with equal coverage for all reads), custom (simulate with diferent negative binomial distributions for novel and known transcripts) or sample (simulate using a real sample)")

    parser_e = subparsers.add_parser("equal", help="\t\tRun in equal mode")
    parser_e.add_argument( "-i", "--trans_index", type=str, help="\t\tFile with transcript information generated with SQANTI-SIM", required=True, )
    parser_e.add_argument( "--gtf", type=str, help="\t\tReference annotation in GTF format", required=True, )
    parser_e.add_argument( "-o", "--output", default=str(), help="\t\tPrefix for output files" )
    parser_e.add_argument( "-d", "--dir", default=".", help="\t\tDirectory for output files (default: .)", )
    parser_e.add_argument( "-nt", "--trans_number", default=10000, type=int, help="\t\tNumber of different transcripts to simulate", )
    parser_e.add_argument( "--read_count", default=50000, type=int, help="\t\tNumber of reads to simulate", )
    parser_e.add_argument( "--ISM", default=0, type=int, help="\t\tNumber of incomplete-splice-matches to delete", )
    parser_e.add_argument( "--NIC", default=0, type=int, help="\t\tNumber of novel-in-catalog to delete", )
    parser_e.add_argument( "--NNC", default=0, type=int, help="\t\tNumber of novel-not-in-catalog to delete", )
    parser_e.add_argument( "--Fusion", default=0, type=int, help="\t\tNumber of Fusion to delete" )
    parser_e.add_argument( "--Antisense", default=0, type=int, help="\t\tNumber of Antisense to delete", )
    parser_e.add_argument( "--GG", default=0, type=int, help="\t\tNumber of Genic-genomic to delete", )
    parser_e.add_argument( "--GI", default=0, type=int, help="\t\tNumber of Genic-intron to delete", )
    parser_e.add_argument( "--Intergenic", default=0, type=int, help="\t\tNumber of Intergenic to delete", )
    parser_e.add_argument( "-k", "--cores", default=1, type=int, help="\t\tNumber of cores to run in parallel", )
    parser_e.add_argument( "-s", "--seed", help="\t\tRandomizer seed", default=None, type=int )

    parser_c = subparsers.add_parser("custom", help="\t\tRun in custom mode")
    parser_c.add_argument( "-i", "--trans_index", type=str, help="\t\tFile with transcript information generated with SQANTI-SIM", required=True, )
    parser_c.add_argument( "--gtf", type=str, help="\t\tReference annotation in GTF format", required=True, )
    parser_c.add_argument( "-o", "--output", default=str(), help="\t\tPrefix for output files" )
    parser_c.add_argument( "-d", "--dir", default=".", help="\t\tDirectory for output files (default: .)", )
    parser_c.add_argument( "-nt", "--trans_number", default=10000, type=int, help="\t\tNumber of different transcripts to simulate", )
    parser_c.add_argument( "--nbn_known", default=15, type=float, help="\t\tAverage read count per known transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)", )
    parser_c.add_argument( "--nbp_known", default=0.5, type=float, help="\t\tThe parameter 'p' of the Negative Binomial distribution for known transcripts", )
    parser_c.add_argument( "--nbn_novel", default=5, type=float, help="\t\tAverage read count per novel transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)", )
    parser_c.add_argument( "--nbp_novel", default=0.5, type=float, help="\t\tThe parameter 'p' of the Negative Binomial distribution for novel transcripts", )
    parser_c.add_argument( "--ISM", default=0, type=int, help="\t\tNumber of incomplete-splice-matches to delete", )
    parser_c.add_argument( "--NIC", default=0, type=int, help="\t\tNumber of novel-in-catalog to delete", )
    parser_c.add_argument( "--NNC", default=0, type=int, help="\t\tNumber of novel-not-in-catalog to delete", )
    parser_c.add_argument( "--Fusion", default=0, type=int, help="\t\tNumber of Fusion to delete" )
    parser_c.add_argument( "--Antisense", default=0, type=int, help="\t\tNumber of Antisense to delete", )
    parser_c.add_argument( "--GG", default=0, type=int, help="\t\tNumber of Genic-genomic to delete", )
    parser_c.add_argument( "--GI", default=0, type=int, help="\t\tNumber of Genic-intron to delete", )
    parser_c.add_argument( "--Intergenic", default=0, type=int, help="\t\tNumber of Intergenic to delete", )
    parser_c.add_argument( "-k", "--cores", default=1, type=int, help="\t\tNumber of cores to run in parallel", )
    parser_c.add_argument( "-s", "--seed", help="\t\tRandomizer seed", default=None, type=int )

    parser_s = subparsers.add_parser("sample", help="\t\tRun in sample mode")
    parser_s.add_argument( "-i", "--trans_index", type=str, help="\t\tFile with transcript information generated with SQANTI-SIM", required=True, )
    parser_s.add_argument( "--gtf", type=str, help="\t\tReference annotation in GTF format", required=True, )
    parser_s.add_argument( "-o", "--output", default=str(), help="\t\tPrefix for output files" )
    parser_s.add_argument( "-d", "--dir", default=".", help="\t\tDirectory for output files (default: .)", )
    parser_s.add_argument( "-nt", "--trans_number", default=None, type=int, help="\t\tNumber of different transcripts to simulate", )
    parser_s.add_argument( "--genome", default=str(), type=str, help="\t\tReference genome FASTA", required=True, )
    group = parser_s.add_mutually_exclusive_group()
    group.add_argument( "--pb_reads", default=str(), type=str, help="\t\tInput PacBio reads for quantification", )
    group.add_argument( "--ont_reads", default=str(), type=str, help="\t\tInput ONT reads for quantification", )
    parser_s.add_argument( "--diff_exp", action="store_true", help="\t\tIf used the program will simulate different expression values for novel and known transcripts", )
    parser_s.add_argument( "--low_prob", default=0.25, type=float, help="\t\tLow value of prob vector (if --diff_exp)", )
    parser_s.add_argument( "--high_prob", default=0.75, type=float, help="\t\tHigh value of prob vector (if --diff_exp)", )
    parser_s.add_argument( "--ISM", default=0, type=int, help="\t\tNumber of incomplete-splice-matches to delete", )
    parser_s.add_argument( "--NIC", default=0, type=int, help="\t\tNumber of novel-in-catalog to delete", )
    parser_s.add_argument( "--NNC", default=0, type=int, help="\t\tNumber of novel-not-in-catalog to delete", )
    parser_s.add_argument( "--Fusion", default=0, type=int, help="\t\tNumber of Fusion to delete" )
    parser_s.add_argument( "--Antisense", default=0, type=int, help="\t\tNumber of Antisense to delete", )
    parser_s.add_argument( "--GG", default=0, type=int, help="\t\tNumber of Genic-genomic to delete", )
    parser_s.add_argument( "--GI", default=0, type=int, help="\t\tNumber of Genic-intron to delete", )
    parser_s.add_argument( "--Intergenic", default=0, type=int, help="\t\tNumber of Intergenic to delete", )
    parser_s.add_argument( "-k", "--cores", default=3, type=int, help="\t\tNumber of cores to run in parallel", )
    parser_s.add_argument( "-s", "--seed", help="\t\tRandomizer seed", default=None, type=int )
    
    args, unknown = parser.parse_known_args(input)

    if unknown:
        print(
            "[SQANTI-SIM] preparatory mode unrecognized arguments: {}\n".format(
                " ".join(unknown)
            ),
            file=sys.stderr,
        )

    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    print("\n[SQANTI-SIM] Running with the following parameters:")
    if args.mode == "equal":
        print("[SQANTI-SIM] - Mode: equal")
        print("[SQANTI-SIM] - Ref GTF:", str(args.gtf))
        print("[SQANTI-SIM] - Out prefix:", str(args.output))
        print("[SQANTI-SIM] - Out dir:", str(args.dir))
        print("[SQANTI-SIM] - N transcripts:", str(args.trans_number))
        print("[SQANTI-SIM] - N reads:", str(args.read_count))

    elif args.mode == "custom":
        print("[SQANTI-SIM] - Mode: custom")
        print("[SQANTI-SIM] - Ref GTF:", str(args.gtf))
        print("[SQANTI-SIM] - Out prefix:", str(args.output))
        print("[SQANTI-SIM] - Out dir:", str(args.dir))
        print("[SQANTI-SIM] - N transcripts:", str(args.trans_number))
        print("[SQANTI-SIM] - Known NB mean count:", str(args.nbn_known))
        print("[SQANTI-SIM] - Known NB probability:", str(args.nbp_known))
        print("[SQANTI-SIM] - Novel NB mean count:", str(args.nbn_novel))
        print("[SQANTI-SIM] - Novel NB probability:", str(args.nbp_novel))

    elif args.mode == "sample":
        print("[SQANTI-SIM] - Mode: sample")
        print("[SQANTI-SIM] - Ref GTF:", str(args.gtf))
        print("[SQANTI-SIM] - Ref genome:", str(args.genome))
        print("[SQANTI-SIM] - Out prefix:", str(args.output))
        print("[SQANTI-SIM] - Out dir:", str(args.dir))
        print("[SQANTI-SIM] - N transcripts:", str(args.trans_number))
        if args.pb_reads:
            print("[SQANTI-SIM] - PacBio reads:", str(args.pb_reads))
        else:
            print("[SQANTI-SIM] - ONT reads:", str(args.ont_reads))
        print("[SQANTI-SIM] - N threads:", str(args.cores))

    if not args.seed():
        args.seed = int.from_bytes(os.urandom(4), 'big')
    random.seed(args.seed)
    numpy.random.seed(args.seed)
    print("[SQANTI-SIM] - Seed:", str(args.seed))

    print("[SQANTI-SIM]\tISM\tNIC\tNNC\tFusion\tAS\tGG\tGI\tIntergenic")
    print("[SQANTI-SIM]\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(
        str(args.ISM), str(args.NIC), str(args.NNC), str(args.Fusion),
        str(args.Antisense), str(args.GG), str(args.GI), str(args.Intergenic)
    ))

    if not args.output:
        output = os.path.basename(args.trans_index).split("_")
        args.output = "_".join(output[:-1])

    # Modify GTF
    print("\n[SQANTI-SIM][%s] Generating modified GTF" %(strftime("%d-%m-%Y %H:%M:%S")))
    counts_end = sim_preparatory.simulate_gtf(args)

    counts_ini = defaultdict(
        lambda: 0,
        {
            "full-splice_match": 0,
            "incomplete-splice_match": args.ISM,
            "novel_in_catalog": args.NIC,
            "novel_not_in_catalog": args.NNC,
            "fusion": args.Fusion,
            "antisense": args.Antisense,
            "genic_intron": args.GI,
            "genic": args.GG,
            "intergenic": args.Intergenic,
        },
    )
    sim_preparatory.summary_table_del(counts_ini, counts_end)

    # Generate expression matrix
    print("[SQANTI-SIM][%s] Generating expression matrix" %(strftime("%d-%m-%Y %H:%M:%S")))
    index_file = os.path.join(args.dir, (args.output + "_index.tsv"))

    if args.mode == "equal":
        sim_preparatory.create_expr_file_fixed_count(index_file, args)

    elif args.mode == "custom":
        sim_preparatory.create_expr_file_nbinom(index_file, args)

    elif args.mode == "sample":
        if args.pb_reads:
            sim_preparatory.create_expr_file_sample(index_file, args, "pb")
        else:
            sim_preparatory.create_expr_file_sample(index_file, args, "ont")

    print("[SQANTI-SIM][%s] preparatory step finished" %(strftime("%d-%m-%Y %H:%M:%S")))
    

def sim(input: list):
    """Simulate reads

    Simulate PacBio and/or ONT reads using the IsoSeqSim or NanoSim pipelines.
    It can also simulate Illumina reads using the polyester pipeline

    Args:
        input (list): arguments to parse
    """

    parser = argparse.ArgumentParser( prog="sqanti_sim.py sim", description="sqanti_sim.py sim parse options" )
    parser.add_argument( "--gtf", type=str, help="\t\tReference annotation in GTF format", required=True, )
    parser.add_argument( "--genome", default=False, help="\t\tReference genome FASTA", required=True, )
    parser.add_argument( "-i", "--trans_index", type=str, help="\t\tFile with transcript information generated with SQANTI-SIM", required=True, )
    parser.add_argument( "--read_type", default="dRNA", type=str, help="\t\tRead type for NanoSim simulation", )
    parser.add_argument( "-d", "--dir", default=".", help="\t\tDirectory for output files (default: .)", )
    parser.add_argument( "-k", "--cores", default=1, type=int, help="\t\tNumber of cores to run in parallel", )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument( "--pb", action="store_true", help="\t\tIf used the program will simulate PacBio reads with IsoSeqSim", )
    group.add_argument( "--ont", action="store_true", help="\t\tIf used the program will simulate ONT reads with NanoSim", )
    parser.add_argument( "--illumina", action="store_true", help="\t\tIf used the program will simulate Illumina reads with RSEM", )
    parser.add_argument( "--long_count", default=None, type=int, help="\t\tNumber of long reads to simulate (if not given it will use the counts of the given expression file)", )
    parser.add_argument( "--short_count", default=None, type=int, help="\t\tNumber of short reads to simulate (if not given it will use the counts of the given expression file)", )
    parser.add_argument( "-s", "--seed", help="\t\tRandomizer seed", default=None, type=int )

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print(
            "[SQANTI-SIM] sim mode unrecognized arguments: {}\n".format(" ".join(unknown)),
            file=sys.stderr,
        )
    
    print("\n[SQANTI-SIM] Running with the following parameters:")
    print("[SQANTI-SIM] - Ref GTF:", str(args.gtf))
    print("[SQANTI-SIM] - Ref genome:", str(args.genome))
    print("[SQANTI-SIM] - Index file:", str(args.trans_index))
    print("[SQANTI-SIM] - Out dir:", str(args.dir))

    if args.ont:
        print("[SQANTI-SIM] - Platform: ONT")
        print("[SQANTI-SIM] - Read type:", str(args.read_type))
    else:
        print("[SQANTI-SIM] - Platform: PacBio")
    
    if args.long_count:
        print("[SQANTI-SIM] - Long reads:", str(args.long_count))
    else:
        print("[SQANTI-SIM] - Long reads: requested_counts from index file")
    
    if args.illumina:
        print("[SQANTI-SIM] - Platform: Illumina")
        if args.short_count:
            print("[SQANTI-SIM] - Short reads:", str(args.short_count))
        else:
            print("[SQANTI-SIM] - Short reads: requested_counts from index file")

    print("[SQANTI-SIM] - N threads:", str(args.cores))

    if not args.seed:
        args.seed = int.from_bytes(os.urandom(4), 'big')
    random.seed(args.seed)
    numpy.random.seed(args.seed)
    print("[SQANTI-SIM] - Seed:", str(args.seed))

    # Simulation with IsoSeqSim, NanoSim and/or Polyester
    if args.pb:
        print("\n[SQANTI-SIM][%s] Simulating PacBio reads" %(strftime("%d-%m-%Y %H:%M:%S")))
        pb_ont_sim.pb_simulation(args)
    if args.ont:
        print("\n[SQANTI-SIM][%s] Simulating ONT reads" %(strftime("%d-%m-%Y %H:%M:%S")))
        pb_ont_sim.ont_simulation(args)
    if args.illumina:
        print("\n[SQANTI-SIM][%s] Simulating Illumina reads" %(strftime("%d-%m-%Y %H:%M:%S")))
        pb_ont_sim.illumina_simulation(args)

    print("[SQANTI-SIM][%s] sim step finished" %(strftime("%d-%m-%Y %H:%M:%S")))


def eval(input: list):
    """Generates SQANTI-SIM report

    Run SQANTI3 with the reconstructed transcripts retrieved by your pipeline
    and generates the SQANTI-SIM report with the evaluation metrics

    Args:
        input (list): arguments to parse
    """

    parser = argparse.ArgumentParser( prog="sqanti_sim.py eval", description="sqanti_sim.py eval parse options", )
    parser.add_argument( "--isoforms", default=str(), help="\t\tGTF with trancriptome reconstructed with your pipeline", required=True, )
    parser.add_argument( "--gtf", type=str, help="\t\tReference annotation in GTF format", required=True, )
    parser.add_argument( "--genome", default=False, help="\t\tReference genome FASTA", required=True, )
    parser.add_argument( "-i", "--trans_index", type=str, help="\t\tFile with transcript information generated with SQANTI-SIM", required=True, )
    parser.add_argument( "-o", "--output", default="sqanti_sim", help="\t\tPrefix for output files", )
    parser.add_argument( "-d", "--dir", default=".", help="\t\tDirectory for output files (default: .)", )
    parser.add_argument( "--short_reads", default=None, help="\t\tFile Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.", required=False, )
    parser.add_argument( "--cage_peak", default=None,help="\t\tFANTOM5 Cage Peak (BED format, optional)" )
    parser.add_argument( "--min_support", default=3, type=int, help="\t\tMinimum number of supporting reads for an isoform", )
    parser.add_argument( "-k", "--cores", default=1, type=int, help="\t\tNumber of cores to run in parallel", )

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print(
            "[SQANTI-SIM] sim mode unrecognized arguments: {}\n".format(" ".join(unknown)),
            file=sys.stderr,
        )

    print("\n[SQANTI-SIM] Running with the following parameters:")
    print("[SQANTI-SIM] - Reconstructed transcripts:", str(args.isoforms))
    print("[SQANTI-SIM] - Modified ref GTF:", str(args.gtf))
    print("[SQANTI-SIM] - Ref genome:", str(args.genome))
    print("[SQANTI-SIM] - Index file:", str(args.trans_index))
    print("[SQANTI-SIM] - Out prefix:", str(args.output))
    print("[SQANTI-SIM] - Out dir:", str(args.dir))

    if args.short_reads:
        print("[SQANTI-SIM] - Short reads:", str(args.short_reads))
    if args.cage_peak:
        print("[SQANTI-SIM] - CAGE Peak:", str(args.cage_peak))
    
    print("[SQANTI-SIM] - Min support:", str(args.min_support))
    print("[SQANTI-SIM] - N threads:", str(args.cores))
    print()

    sqanti3_stats.sqanti3_stats(args)

    print("[SQANTI-SIM][%s] eval step finished" %(strftime("%d-%m-%Y %H:%M:%S")))


#####################################
#                                   #
#               MAIN                #
#                                   #
#####################################

print(
    """                                                                      
      _____  ____            _   _ _______ _____      _____ _____ __  __  
     / ____|/ __ \     /\   | \ | |__   __|_   _|    / ____|_   _|  \/  | 
    | (___ | |  | |   /  \  |  \| |  | |    | |_____| (___   | | | \  / | 
     \___ \| |  | |  / /\ \ | . ` |  | |    | |______\___ \  | | | |\/| | 
     ____) | |__| | / ____ \| |\  |  | |   _| |_     ____) |_| |_| |  | | 
    |_____/ \___\_\/_/    \_\_| \_|  |_|  |_____|   |_____/|_____|_|  |_| 
                                                                          
              A SIMULATOR OF CONTROLLED NOVELTY AND DEGRADATION           
                    OF TRANSCRIPTS SEQUENCED BY LONG-READS                
    """
)

if len(sys.argv) < 2:
    print("[SQANTI-SIM] usage: python sqanti_sim.py <mode> --help\n", file=sys.stderr)
    print("[SQANTI-SIM] modes: classif, preparatory, sim, eval\n", file=sys.stderr)
    sys.exit(1)

else:
    mode = sys.argv[1].lower()
    input = sys.argv[2:]

if mode == "classif":
    print("[SQANTI-SIM] CLASSIF MODE")
    res = classif(input)

elif mode == "preparatory":
    print("[SQANTI-SIM] PREPARATORY MODE")
    res = preparatory(input)

elif mode == "sim":
    print("[SQANTI-SIM] SIM MODE")
    res = sim(input)

elif mode == "eval":
    print("[SQANTI-SIM] EVAL MODE")
    res = eval(input)

elif mode in ["--version", "-v"]:
    print("[SQANTI-SIM] SQANTI-SIM %s\n" %(__version__))

else:
    print("[SQANTI-SIM] usage: python sqanti_sim.py <mode> --help\n", file=sys.stderr)
    print("[SQANTI-SIM] modes: classif, preparatory, sim, eval\n", file=sys.stderr)
    sys.exit(1)

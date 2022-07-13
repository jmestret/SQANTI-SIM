#!/usr/bin/env python3
"""
sqantisim.py

Wrapper for long-read RNA-seq simulators (NanoSim and IsoSeqSim) to simulate
controlled novelty and degradation of transcripts based on SQANTI3 structural
categories

Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)
"""

__version__ = "0.1.0-beta"

import argparse
import numpy
import os
import random
import sys
from collections import defaultdict
from src import classify_gtf
from src import simulate_reads
from src import design_simulation
from src import evaluation_metrics
from time import strftime


def classif(input: list):
    """Classify transcripts in SQANTI3 structural categories

    Given a GTF annotation generates an index file with the potential SQANTI3
    structural category of each transcript

    Args:
        input (list): arguments to parse
    """

    parser = argparse.ArgumentParser( prog="sqantisim.py classif", description="sqantisim.py classif parse options", )
    parser.add_argument("--gtf", type=str, required=True, help="\t\tReference annotation in GTF format", )
    parser.add_argument("-o", "--output", type=str, default="sqantisim", help="\t\tPrefix for output file", )
    parser.add_argument("-d", "--dir", type=str, default=".", help="\t\tDirectory for output files (default: .)", )
    parser.add_argument("-k", "--cores", type=int, default=1, help="\t\tNumber of cores to run in parallel", )

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print("[SQANTISIM] classif mode unrecognized arguments: {}\n".format(" ".join(unknown)),file=sys.stderr)
    
    if not os.path.exists(args.gtf):
        print("[SQANTISIM] ERROR: --gtf file does not exist. Provide a valid path", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)

    print("\n[SQANTISIM] Running with the following parameters:")
    print("[SQANTISIM] - Ref GTF:", str(args.gtf))
    print("[SQANTISIM] - Out prefix:", str(args.output))
    print("[SQANTISIM] - Out dir:", str(args.dir))
    print("[SQANTISIM] - N threads:", str(args.cores))

    print("\n[SQANTISIM][%s] Classifying transcripts in structural categories" %(strftime("%d-%m-%Y %H:%M:%S")))
    trans_info = classify_gtf.classify_gtf(args)
    
    print("[SQANTISIM] Summary table from categorization")
    classify_gtf.summary_table_cat(trans_info)

    print("[SQANTISIM][%s] classif step finished" %(strftime("%d-%m-%Y %H:%M:%S")))


def design(input: list):
    """Modifies reference annotation GTF and builds expression matrix

    Given the novel and known transcripts to simulate and its counts, generates
    the expression matrix to give as input to the long-read RNA-seq simulators
    and generetes the modified GTF to use as reference annotation in your tool

    Args:
        input (list): arguments to parse
    """
    parser = argparse.ArgumentParser(prog="sqantisim.py design", description="sqantisim.py design parse options", )
    subparsers = parser.add_subparsers(dest="mode", description="\t\tDifferent modes to generate the expression matrix: equal (simulate with equal coverage for all reads), custom (simulate with diferent negative binomial distributions for novel and known transcripts) or sample (simulate using a real sample)")

    parser_e = subparsers.add_parser("equal", help="\t\tRun in equal mode")
    parser_e.add_argument("-i", "--trans_index", type=str, required=True, help="\t\tFile with transcript information generated with SQANTISIM", )
    parser_e.add_argument("--gtf", type=str, required=True, help="\t\ttComplete reference annotation in GTF format", )
    parser_e.add_argument("-o", "--output", type=str, default=str(), help="\t\tPrefix for output files" )
    parser_e.add_argument("-d", "--dir", type=str, default=".", help="\t\tDirectory for output files (default: .)", )
    parser_e.add_argument("-nt", "--trans_number", type=int, default=10000, help="\t\tNumber of different transcripts to simulate", )
    parser_e.add_argument("--read_count", default=50000, type=int, help="\t\tNumber of reads to simulate", )
    parser_e.add_argument("--ISM", type=int, default=0, help="\t\tNumber of incomplete-splice-matches to delete", )
    parser_e.add_argument("--NIC", type=int, default=0, help="\t\tNumber of novel-in-catalog to delete", )
    parser_e.add_argument("--NNC", type=int, default=0, help="\t\tNumber of novel-not-in-catalog to delete", )
    parser_e.add_argument("--Fusion", type=int, default=0, help="\t\tNumber of Fusion to delete" )
    parser_e.add_argument("--Antisense", type=int, default=0, help="\t\tNumber of Antisense to delete", )
    parser_e.add_argument("--GG", type=int, default=0, help="\t\tNumber of Genic-genomic to delete", )
    parser_e.add_argument("--GI", type=int, default=0, help="\t\tNumber of Genic-intron to delete", )
    parser_e.add_argument("--Intergenic", type=int, default=0, help="\t\tNumber of Intergenic to delete", )
    parser_e.add_argument("-k", "--cores", type=int, default=1, help="\t\tNumber of cores to run in parallel", )
    parser_e.add_argument("-s", "--seed", type=int, default=None, help="\t\tRandomizer seed", )

    parser_c = subparsers.add_parser("custom", help="\t\tRun in custom mode")
    parser_c.add_argument("-i", "--trans_index", type=str, required=True, help="\t\tFile with transcript information generated with SQANTISIM", )
    parser_c.add_argument("--gtf", type=str, required=True, help="\t\ttComplete reference annotation in GTF format", )
    parser_c.add_argument("-o", "--output", type=str, default=str(), help="\t\tPrefix for output files" )
    parser_c.add_argument("-d", "--dir", type=str, default=".", help="\t\tDirectory for output files (default: .)", )
    parser_c.add_argument("-nt", "--trans_number", type=int, default=10000, help="\t\tNumber of different transcripts to simulate", )
    parser_c.add_argument("--nbn_known", type=float, default=15, help="\t\tAverage read count per known transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)", )
    parser_c.add_argument("--nbp_known", type=float, default=0.5, help="\t\tThe parameter 'p' of the Negative Binomial distribution for known transcripts", )
    parser_c.add_argument("--nbn_novel", type=float, default=5, help="\t\tAverage read count per novel transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution)", )
    parser_c.add_argument("--nbp_novel", type=float, default=0.5, help="\t\tThe parameter 'p' of the Negative Binomial distribution for novel transcripts", )
    parser_c.add_argument("--ISM", type=int, default=0, help="\t\tNumber of incomplete-splice-matches to delete", )
    parser_c.add_argument("--NIC", type=int, default=0, help="\t\tNumber of novel-in-catalog to delete", )
    parser_c.add_argument("--NNC", type=int, default=0, help="\t\tNumber of novel-not-in-catalog to delete", )
    parser_c.add_argument("--Fusion", type=int, default=0, help="\t\tNumber of Fusion to delete" )
    parser_c.add_argument("--Antisense", type=int, default=0, help="\t\tNumber of Antisense to delete", )
    parser_c.add_argument("--GG", type=int, default=0, help="\t\tNumber of Genic-genomic to delete", )
    parser_c.add_argument("--GI", type=int, default=0, help="\t\tNumber of Genic-intron to delete", )
    parser_c.add_argument("--Intergenic", type=int, default=0, help="\t\tNumber of Intergenic to delete", )
    parser_c.add_argument("-k", "--cores", type=int, default=1, help="\t\tNumber of cores to run in parallel", )
    parser_c.add_argument("-s", "--seed", type=int, default=None, help="\t\tRandomizer seed", )

    parser_s = subparsers.add_parser("sample", help="\t\tRun in sample mode")
    parser_s.add_argument("-i", "--trans_index", type=str, required=True, help="\t\tFile with transcript information generated with SQANTISIM", )
    parser_s.add_argument("--gtf", type=str, required=True, help="\t\tComplete reference annotation in GTF format", )
    parser_s.add_argument("-o", "--output", type=str, default=str(), help="\t\tPrefix for output files" )
    parser_s.add_argument("-d", "--dir", type=str, default=".", help="\t\tDirectory for output files (default: .)", )
    parser_s.add_argument("-nt", "--trans_number", type=int, default=None, help="\t\tNumber of different transcripts to simulate", )
    parser_s.add_argument("--genome", type=str, required=True, help="\t\tReference genome FASTA", )
    group = parser_s.add_mutually_exclusive_group()
    group.add_argument("--pb_reads", type=str, default=str(), help="\t\tInput PacBio reads for quantification in FASTA or FASTQ format", )
    group.add_argument("--ont_reads", type=str, default=str(), help="\t\tInput ONT reads for quantification in FASTA or FASTQ format", )
    group.add_argument("--mapped_reads", type=str, default=str(), help="\t\tAligned reads in SAM format", )
    parser_s.add_argument("--iso_complex", action="store_true", help="\t\tIf used the program will simulate the expressed isoform complexity (number of isoforms per gene)", )
    parser_s.add_argument("--diff_exp", action="store_true", help="\t\tIf used the program will simulate different expression values for novel and known transcripts", )
    parser_s.add_argument("--low_prob", type=float, default=0.1, help="\t\tLow value of prob vector (if --diff_exp)", )
    parser_s.add_argument("--high_prob", type=float, default=0.9, help="\t\tHigh value of prob vector (if --diff_exp)", )
    parser_s.add_argument("--ISM", type=int, default=0, help="\t\tNumber of incomplete-splice-matches to delete", )
    parser_s.add_argument("--NIC", type=int, default=0, help="\t\tNumber of novel-in-catalog to delete", )
    parser_s.add_argument("--NNC", type=int, default=0, help="\t\tNumber of novel-not-in-catalog to delete", )
    parser_s.add_argument("--Fusion", type=int, default=0, help="\t\tNumber of Fusion to delete" )
    parser_s.add_argument("--Antisense", type=int, default=0, help="\t\tNumber of Antisense to delete", )
    parser_s.add_argument("--GG", type=int, default=0, help="\t\tNumber of Genic-genomic to delete", )
    parser_s.add_argument("--GI", type=int, default=0, help="\t\tNumber of Genic-intron to delete", )
    parser_s.add_argument("--Intergenic", type=int, default=0, help="\t\tNumber of Intergenic to delete", )
    parser_s.add_argument("-k", "--cores", type=int, default=1, help="\t\tNumber of cores to run in parallel", )
    parser_s.add_argument("-s", "--seed", type=int, default=None, help="\t\tRandomizer seed", )
    
    args, unknown = parser.parse_known_args(input)

    if unknown:
        print("[SQANTISIM] design mode unrecognized arguments: {}\n".format(" ".join(unknown)), file=sys.stderr)

    total_novel = sum([args.ISM, args.NIC, args.NNC, args.Fusion, args.Antisense, args.GG, args.GI, args.Intergenic])
    if args.trans_number is not None and total_novel > args.trans_number:
        print("[SQANTISIM] WARNING: -nt is lower than the novel transcripts to simulate, only novel transcripts will be simulated", file=sys.stderr)

    if not os.path.isdir(args.dir):
        os.makedirs(args.dir)
    
    if not args.output:
        output = os.path.basename(args.trans_index).split("_")
        args.output = "_".join(output[:-1])

    if not args.seed:
        args.seed = int.from_bytes(os.urandom(1), 'big')
    random.seed(args.seed)
    numpy.random.seed(args.seed)

    print("\n[SQANTISIM] Running with the following parameters:")
    if args.mode == "equal":
        print("[SQANTISIM] - Mode: equal")
        print("[SQANTISIM] - Ref GTF:", str(args.gtf))
        print("[SQANTISIM] - Out prefix:", str(args.output))
        print("[SQANTISIM] - Out dir:", str(args.dir))
        print("[SQANTISIM] - N transcripts:", str(args.trans_number))
        print("[SQANTISIM] - N reads:", str(args.read_count))

    elif args.mode == "custom":
        if args.nbn_known < 0 or args.nbn_novel < 0:
            print("[SQANTISIM] ERROR: --nbn_known and --nbn_novel must be greater than 0", file=sys.stderr)
            sys.exit(1)
        if args.nbp_known < 0 or args.nbp_known > 1 or args.nbp_novel < 0 or args.nbp_novel > 1:
            print("[SQANTISIM] ERROR: --nbp_known and --nbp_novel must be in the interval [0,1]", file=sys.stderr)
            sys.exit(1)

        print("[SQANTISIM] - Mode: custom")
        print("[SQANTISIM] - Ref GTF:", str(args.gtf))
        print("[SQANTISIM] - Out prefix:", str(args.output))
        print("[SQANTISIM] - Out dir:", str(args.dir))
        print("[SQANTISIM] - N transcripts:", str(args.trans_number))
        print("[SQANTISIM] - Known NB mean count:", str(args.nbn_known))
        print("[SQANTISIM] - Known NB probability:", str(args.nbp_known))
        print("[SQANTISIM] - Novel NB mean count:", str(args.nbn_novel))
        print("[SQANTISIM] - Novel NB probability:", str(args.nbp_novel))

    elif args.mode == "sample":
        if args.low_prob < 0 or args.low_prob > 1 or args.high_prob < 0 or args.high_prob > 1:
            print("[SQANTISIM] ERROR: --low_prob and --high_prob must be in the interval [0,1]", file=sys.stderr)
            sys.exit(1)

        print("[SQANTISIM] - Mode: sample")
        print("[SQANTISIM] - Ref GTF:", str(args.gtf))
        print("[SQANTISIM] - Ref genome:", str(args.genome))
        print("[SQANTISIM] - Out prefix:", str(args.output))
        print("[SQANTISIM] - Out dir:", str(args.dir))
        print("[SQANTISIM] - N transcripts:", str(args.trans_number))
        if args.pb_reads:
            print("[SQANTISIM] - PacBio reads:", str(args.pb_reads))
        else:
            print("[SQANTISIM] - ONT reads:", str(args.ont_reads))
        print("[SQANTISIM] - N threads:", str(args.cores))

    print("[SQANTISIM] - Seed:", str(args.seed))
    print("[SQANTISIM]\tISM\tNIC\tNNC\tFusion\tAS\tGG\tGI\tIntergenic")
    print("[SQANTISIM]\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(
        str(args.ISM), str(args.NIC), str(args.NNC), str(args.Fusion),
        str(args.Antisense), str(args.GG), str(args.GI), str(args.Intergenic)
    ))

    # Modify GTF
    print("\n[SQANTISIM][%s] Generating modified GTF" %(strftime("%d-%m-%Y %H:%M:%S")))
    counts_end = design_simulation.simulate_gtf(args)

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
    design_simulation.summary_table_del(counts_ini, counts_end)

    # Generate expression matrix
    print("[SQANTISIM][%s] Generating expression matrix" %(strftime("%d-%m-%Y %H:%M:%S")))
    index_file = os.path.join(args.dir, (args.output + "_index.tsv"))

    if args.mode == "equal":
        design_simulation.create_expr_file_fixed_count(index_file, args)

    elif args.mode == "custom":
        design_simulation.create_expr_file_nbinom(index_file, args)

    elif args.mode == "sample":
        if args.pb_reads:
            design_simulation.create_expr_file_sample(index_file, args, "pb")
        else:
            design_simulation.create_expr_file_sample(index_file, args, "ont")

    print("[SQANTISIM][%s] design step finished" %(strftime("%d-%m-%Y %H:%M:%S")))
    

def sim(input: list):
    """Simulate reads

    Simulate PacBio and/or ONT reads using the IsoSeqSim or NanoSim pipelines.
    It can also simulate Illumina reads using the polyester pipeline

    Args:
        input (list): arguments to parse
    """

    parser = argparse.ArgumentParser( prog="sqantisim.py sim", description="sqantisim.py sim parse options" )
    parser.add_argument("--gtf", type=str, required=True, help="\t\tReduced reference annotation in GTF format", )
    parser.add_argument("--genome", type=str, required=True, help="\t\tReference genome FASTA", )
    parser.add_argument("-i", "--trans_index", type=str, required=True, help="\t\tFile with transcript information generated with SQANTISIM", )
    parser.add_argument("--read_type", type=str, default="dRNA", help="\t\tRead type for NanoSim simulation", )
    parser.add_argument("-d", "--dir", type=str, default=".", help="\t\tDirectory for output files (default: .)", )
    parser.add_argument("-k", "--cores", type=int, default=1, help="\t\tNumber of cores to run in parallel", )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pb", action="store_true", help="\t\tIf used the program will simulate PacBio reads with IsoSeqSim", )
    group.add_argument("--ont", action="store_true", help="\t\tIf used the program will simulate ONT reads with NanoSim", )
    parser.add_argument("--illumina", action="store_true", help="\t\tIf used the program will simulate Illumina reads with RSEM", )
    parser.add_argument("--long_count", type=int, default=None, help="\t\tNumber of long reads to simulate (if not given it will use the counts of the given expression file)", )
    parser.add_argument("--short_count", type=int, default=None, help="\t\tNumber of short reads to simulate (if not given it will use the counts of the given expression file)", )
    parser.add_argument("-s", "--seed", type=int, default=None, help="\t\tRandomizer seed", )

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print("[SQANTISIM] sim mode unrecognized arguments: {}\n".format(" ".join(unknown)), file=sys.stderr)

    if not args.seed:
        args.seed = int.from_bytes(os.urandom(1), 'big')
    random.seed(args.seed)
    numpy.random.seed(args.seed)
    
    print("\n[SQANTISIM] Running with the following parameters:")
    print("[SQANTISIM] - Ref GTF:", str(args.gtf))
    print("[SQANTISIM] - Ref genome:", str(args.genome))
    print("[SQANTISIM] - Index file:", str(args.trans_index))
    print("[SQANTISIM] - Out dir:", str(args.dir))

    if args.ont:
        print("[SQANTISIM] - Platform: ONT")
        print("[SQANTISIM] - Read type:", str(args.read_type))
    else:
        print("[SQANTISIM] - Platform: PacBio")
    
    if args.long_count:
        print("[SQANTISIM] - Long reads:", str(args.long_count))
    else:
        print("[SQANTISIM] - Long reads: requested_counts from index file")
    
    if args.illumina:
        print("[SQANTISIM] - Platform: Illumina")
        if args.short_count:
            print("[SQANTISIM] - Short reads:", str(args.short_count))
        else:
            print("[SQANTISIM] - Short reads: requested_counts from index file")

    print("[SQANTISIM] - N threads:", str(args.cores))
    print("[SQANTISIM] - Seed:", str(args.seed))

    # Simulation with IsoSeqSim, NanoSim and/or Polyester
    if args.pb:
        print("\n[SQANTISIM][%s] Simulating PacBio reads" %(strftime("%d-%m-%Y %H:%M:%S")))
        simulate_reads.pb_simulation(args)
    if args.ont:
        print("\n[SQANTISIM][%s] Simulating ONT reads" %(strftime("%d-%m-%Y %H:%M:%S")))
        simulate_reads.ont_simulation(args)
    if args.illumina:
        print("\n[SQANTISIM][%s] Simulating Illumina reads" %(strftime("%d-%m-%Y %H:%M:%S")))
        simulate_reads.illumina_simulation(args)

    print("[SQANTISIM][%s] sim step finished" %(strftime("%d-%m-%Y %H:%M:%S")))


def eval(input: list):
    """Generates SQANTISIM report

    Run SQANTI3 with the reconstructed transcripts retrieved by your pipeline
    and generates the SQANTISIM report with the evaluation metrics

    Args:
        input (list): arguments to parse
    """

    parser = argparse.ArgumentParser( prog="sqantisim.py eval", description="sqantisim.py eval parse options", )
    parser.add_argument("--isoforms", type=str, required=True, help="\t\tLong-read trancriptome reconstructed with your pipeline in GTF, FASTA or FASTQ format", )
    parser.add_argument("--gtf", type=str, required=True, help="\t\tReference annotation in GTF format", )
    parser.add_argument("--genome", type=str, required=True, help="\t\tReference genome FASTA", )
    parser.add_argument("-i", "--trans_index", type=str, required=True, help="\t\tFile with transcript information generated with SQANTISIM", )
    parser.add_argument("-o", "--output", type=str, default="sqantisim", help="\t\tPrefix for output files", )
    parser.add_argument("-d", "--dir", type=str, default=".", help="\t\tDirectory for output files (default: .)", )
    parser.add_argument('-c','--coverage', help='\t\tJunction coverage files (provide a single file, comma-delmited filenames, or a file pattern, ex: "mydir/*.junctions").', required=False)
    parser.add_argument('--SR_bam' , help='\t\t Directory or fofn file with the sorted bam files of Short Reads RNA-Seq mapped against the genome', required=False)
    parser.add_argument("--short_reads", type=str, default=None, help="\t\tFile Of File Names (fofn, space separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq. If expression or coverage files are not provided, Kallisto (just for pair-end data) and STAR, respectively, will be run to calculate them.",)
    parser.add_argument("--CAGE_peak", type=str, default=None,help="\t\tFANTOM5 Cage Peak (BED format, optional)" )
    parser.add_argument("--fasta", action="store_true", help="\t\tUse when running SQANTI by using as input a FASTA/FASTQ with the sequences of isoforms", )
    parser.add_argument("--aligner_choice", type=str, default="minimap2",help="\t\tIf --fasta used, choose the aligner to map your isoforms", choices=["minimap2","deSALT","gmap","uLTRA"])
    parser.add_argument("--min_support", type=int, default=3, help="\t\tMinimum number of supporting reads for an isoform", )
    parser.add_argument("-k", "--cores", type=int, default=1, help="\t\tNumber of cores to run in parallel", )

    args, unknown = parser.parse_known_args(input)

    if unknown:
        print(
            "[SQANTISIM] sim mode unrecognized arguments: {}\n".format(" ".join(unknown)),
            file=sys.stderr,
        )

    print("\n[SQANTISIM] Running with the following parameters:")
    print("[SQANTISIM] - Reconstructed transcripts:", str(args.isoforms))
    print("[SQANTISIM] - Modified ref GTF:", str(args.gtf))
    print("[SQANTISIM] - Ref genome:", str(args.genome))
    print("[SQANTISIM] - Index file:", str(args.trans_index))
    print("[SQANTISIM] - Out prefix:", str(args.output))
    print("[SQANTISIM] - Out dir:", str(args.dir))

    if args.coverage:
        print("[SQANTISIM] - Coverage:", str(args.coverage))
    if args.SR_bam:
        print("[SQANTISIM] - Short-read BAM files:", str(args.SR_bam))
    if args.short_reads:
        print("[SQANTISIM] - Short reads:", str(args.short_reads))
    if args.CAGE_peak:
        print("[SQANTISIM] - CAGE Peak:", str(args.CAGE_peak))
    
    print("[SQANTISIM] - Min support:", str(args.min_support))
    print("[SQANTISIM] - N threads:", str(args.cores))
    print()

    evaluation_metrics.sqanti3_stats(args)

    print("[SQANTISIM][%s] eval step finished" %(strftime("%d-%m-%Y %H:%M:%S")))


#####################################
#                                   #
#               MAIN                #
#                                   #
#####################################

print(
    """                                                                      
      _____  ____            _   _ _______ _____  _____ _____ __  __  
     / ____|/ __ \     /\   | \ | |__   __|_   _|/ ____|_   _|  \/  | 
    | (___ | |  | |   /  \  |  \| |  | |    | | | (___   | | | \  / | 
     \___ \| |  | |  / /\ \ | . ` |  | |    | |  \___ \  | | | |\/| | 
     ____) | |__| | / ____ \| |\  |  | |   _| |_ ____) |_| |_| |  | | 
    |_____/ \___\_\/_/    \_\_| \_|  |_|  |_____|_____/|_____|_|  |_| 
                                                                          
            A SIMULATOR OF CONTROLLED NOVELTY AND DEGRADATION           
                 OF TRANSCRIPTS SEQUENCED BY LONG-READS                
    """
)

if len(sys.argv) < 2:
    print("[SQANTISIM] usage: python sqantisim.py <mode> --help\n", file=sys.stderr)
    print("[SQANTISIM] modes: classif, design, sim, eval\n", file=sys.stderr)
    sys.exit(1)

else:
    mode = sys.argv[1].lower()
    input = sys.argv[2:]

if mode == "classif":
    print("[SQANTISIM] CLASSIF MODE")
    res = classif(input)

elif mode == "design":
    print("[SQANTISIM] DESIGN MODE")
    res = design(input)

elif mode == "sim":
    print("[SQANTISIM] SIM MODE")
    res = sim(input)

elif mode == "eval":
    print("[SQANTISIM] EVAL MODE")
    res = eval(input)

elif mode in ["--version", "-v"]:
    print("[SQANTISIM] SQANTISIM %s\n" %(__version__))

else:
    print("[SQANTISIM] usage: python sqantisim.py <mode> --help\n", file=sys.stderr)
    print("[SQANTISIM] modes: classif, design, sim, eval\n", file=sys.stderr)
    sys.exit(1)

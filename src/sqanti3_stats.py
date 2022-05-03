#!/usr/bin/env python3
"""
sqanti3_stats.py

Generate SQANTI-SIM report and metrics

Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)
"""

import os
import pandas
import subprocess
import sys
from collections import defaultdict
from src.SQANTI3.utilities.short_reads import get_TSS_bed, get_ratio_TSS
from src.SQANTI3.sqanti3_qc import CAGEPeak, STARcov_parser
from time import strftime


def sqanti3_stats(args):
    """Runs SQANTI3 and generates SQANTI-SIM report

    Given the reconstructed transcripts in GTF format it runs the SQANTI3
    pipeline and computes the SQANTI-SIM metrics to evaluate how well it
    recovered the novel transcripts

    Args:
        args (list): arguments to parse
    """

    def write_whithin_cage(row):
        return within_cage_dict[row["transcript_id"]]

    def write_dist_cage(row):
        return dist_cage_dict[row["transcript_id"]]

    def write_ratio_TSS(row):
        if row["transcript_id"] in ratio_TSS_dict:
            return ratio_TSS_dict[row["transcript_id"]]["max_ratio_TSS"]
        else:
            return 1

    def write_SJ_cov(row):
        min_cov = "NA"
        if row["exons"] == 1:
            return min_cov
        d = row["donors"].split(",")
        a = row["acceptors"].split(",")
        for i in range(int(row["exons"]) - 1):
            sample_cov = SJcovInfo[row["chrom"], row["strand"]][(int(d[i]), (int(a[i])-1))] # make exon starts (SJ acceptors) 0 based
            total_coverage_unique = (
                sum(
                    [cov_uniq for (cov_uniq, cov_multi) in sample_cov.values()]
                )
                if SJcovInfo is not None
                else "NA"
            )
            if min_cov == "NA" or min_cov > total_coverage_unique:
                min_cov = total_coverage_unique
        return min_cov

    print("[SQANTI-SIM][%s] Running SQANTI3" %(strftime("%d-%m-%Y %H:%M:%S")))
    src_dir = os.path.dirname(os.path.realpath(__file__))
    sqanti3 = os.path.join(src_dir, "SQANTI3/sqanti3_qc.py")

    MIN_REF_LEN = 0
    cmd = [
        sqanti3,
        args.isoforms,
        args.gtf,
        args.genome,
        "-o",
        args.output,
        "-d",
        args.dir,
        "--cpus",
        str(args.cores),
        "--min_ref_len",
        str(MIN_REF_LEN),
        "--force_id_ignore",
        "--skipORF",
    ]

    if args.cage_peak:
        cmd.append("--cage_peak")
        cmd.append(args.cage_peak)

    if args.short_reads:
        cmd.append("--short_reads")
        cmd.append(args.short_reads)

    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("[SQANTI-SIM] ERROR running SQANTI3: {0}".format(cmd), file=sys.stderr)
        #sys.exit(1)

    trans_index = pandas.read_csv(args.trans_index, sep="\t", header=0, dtype={"chrom":str})
    if args.cage_peak:
        print("[SQANTI-SIM][%s] Parsing CAGE Peak data" %(strftime("%d-%m-%Y %H:%M:%S")))
        cage_peak_data = CAGEPeak(args.cage_peak)

        within_cage_dict = defaultdict(lambda: False)
        dist_cage_dict = defaultdict(lambda: False)
        with open(args.trans_index, "r") as index_file:
            header_names = index_file.readline()
            header_names = header_names.split()
            id_pos = header_names.index("transcript_id")
            chrom_pos = header_names.index("chrom")
            strand_pos = header_names.index("strand")
            start_pos = header_names.index("TSS_genomic_coord") # No need to swap start and end coordinates -> already swapped in the index file for negative strand
            for line in index_file:
                line = line.split()
                within_cage, dist_cage = cage_peak_data.find(
                    line[chrom_pos], line[strand_pos], (int(line[start_pos])-1)
                ) # 0 based TSS
                within_cage_dict[line[id_pos]] = within_cage
                dist_cage_dict[line[id_pos]] = dist_cage
        index_file.close()

        trans_index["dist_to_cage_peak"] = trans_index.apply(
            write_dist_cage, axis=1
        )
        trans_index["within_cage_peak"] = trans_index.apply(
            write_whithin_cage, axis=1
        )

    if args.short_reads:
        print("[SQANTI-SIM][%s] Parsing Short Read data" %(strftime("%d-%m-%Y %H:%M:%S")))
        star_out = os.path.join(args.dir, "STAR_mapping/")
        star_index = os.path.join(args.dir, "STAR_index/")

        # Short Read Coverage
        SJcovNames, SJcovInfo = STARcov_parser(star_out)
        trans_index["min_cov"] = trans_index.apply(write_SJ_cov, axis=1)

        # Short reads ratio TSS
        chr_order = os.path.join(star_index, "chrNameLength.txt")
        inside_bed, outside_bed = get_TSS_bed(args.gtf, chr_order)
        bams = []
        for filename in os.listdir(star_out):
            if filename.endswith(".bam"):
                bams.append(star_out + "/" + filename)
        ratio_TSS_dict = get_ratio_TSS(
            inside_bed, outside_bed, bams, chr_order
        )
        trans_index["ratio_TSS"] = trans_index.apply(write_ratio_TSS, axis=1)

    trans_index.to_csv(
        args.trans_index, sep="\t", na_rep="NA", header=True, index=False
    )

    print("[SQANTI-SIM][%s] Generating SQANTI-SIM report" %(strftime("%d-%m-%Y %H:%M:%S")))
    src_dir = os.path.dirname(os.path.realpath(__file__))
    classification_file = os.path.join(
        args.dir, (args.output + "_classification.txt")
    )
    junctions_file = os.path.join(args.dir, (args.output + "_junctions.txt"))

    cmd = [
        "Rscript",
        os.path.join(src_dir, "SQANTI_SIM_report.R"),
        classification_file,
        junctions_file,
        args.trans_index,
        args.min_support,
        src_dir,
    ]

    cmd = " ".join(cmd)
    if subprocess.check_call(cmd, shell=True) != 0:
        print(
            "[SQANTI-SIM] ERROR running SQANTI-SIM report generation: {0}".format(cmd),
            file=sys.stderr,
        )
        sys.exit(1)

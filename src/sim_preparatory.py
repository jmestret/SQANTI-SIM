#!/usr/bin/env python3
"""
sim_preparatory.py

Modifies GTF reference annotation and generate expression values for 
the simulation step

Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)
"""

from bisect import bisect_left
from calendar import c
from email.policy import default
import numpy
import os
import pandas
import pysam
import random
import subprocess
import sys
from bisect import bisect_left
from collections import defaultdict

MIN_SIM_LEN = 200 # Minimum length of transcripts to simulate

def target_trans(f_idx: str, f_idx_out: str, counts: dict) -> tuple:
    """
    Choose those transcripts that will be deleted from the original GTF
    to generate the modified file that will be used as the reference annotation

    Args:
        f_idx (str): name of the input transcript index file
        f_idx_out (str): name of the output transcript index file
        counts (dict): dictinary with the number of transcripts of each
                       structural category to be deleted
    Returns:
        final_target (set): all transcripts to be simulated (deleted from GTF)
    """

    def pick_sim_type(row):
        if row["transcript_id"] in target_trans:
            return "novel"
        else:
            return "known"

    trans_by_SC = defaultdict(lambda: [])
    trans_by_gene = defaultdict(lambda: [])

    target_trans = set()
    target_genes = set()
    ref_trans = set()
    ref_genes = set()

    # Build a dict with all transcripts classified in each structural category
    with open(f_idx, "r") as cat:
        col_names = cat.readline()
        col_names = col_names.split()
        for line in cat:
            line_split = line.split()
            gene = line_split[1]
            SC = line_split[2]

            trans_by_SC[SC].append(tuple(line_split))
            trans_by_gene[gene].append(tuple(line_split))

    cat.close()

    # Select randomly the transcripts of each SC that are going to be deleted
    # It's important to make sure you don't delete its reference trans or gene
    categories = list(counts.keys())
    weight_list = []
    for SC in categories:
        weight_list.append(len(trans_by_SC[SC]))
        random.shuffle(trans_by_SC[SC])

    while categories:
        SC = random.choices(categories, weights=weight_list, k=1)
        SC = SC[0]
        if counts[SC] <= 0 or len(trans_by_SC[SC]) == 0:
            i = categories.index(SC)
            del categories[i]
            del weight_list[i]
        else:
            trans = trans_by_SC[SC].pop()
            trans_id = trans[0]
            gene_id = trans[1]
            SC = trans[2]
            ref_g = trans[3]
            ref_t = trans[4]
            TSS = int(trans[col_names.index("TSS_genomic_coord")])
            TTS = int(trans[col_names.index("TTS_genomic_coord")])

            if abs(TSS - TTS) <= MIN_SIM_LEN: # Dont simulate small transcripts
                continue

            if SC in ["full-splice_match", "incomplete-splice_match"]:
                if (
                    trans_id not in ref_trans
                    and gene_id not in ref_genes
                    and ref_t not in target_trans
                ):
                    target_trans.add(trans_id)
                    target_genes.add(gene_id)
                    ref_trans.add(ref_t)
                    counts[SC] -= 1

            elif SC in [ "novel_not_in_catalog", "genic_intron"]:
                if (
                    trans_id not in ref_trans
                    and gene_id not in ref_genes
                    and gene_id not in target_genes
                    and ref_g not in target_genes
                ):
                    target_trans.add(trans_id)
                    target_genes.add(gene_id)
                    ref_genes.add(ref_g)
                    counts[SC] -= 1

            elif SC in ["novel_in_catalog", "fusion", "antisense", "genic"]:
                if (
                    trans_id not in ref_trans
                    and gene_id not in ref_genes
                    and gene_id not in target_genes
                ):
                    ref_g = trans[3].split("_")
                    for i in ref_g:
                        if i in target_genes:
                            break
                    else:
                        target_trans.add(trans_id)
                        target_genes.add(gene_id)
                        for i in ref_g:
                            ref_genes.add(i)
                        counts[SC] -= 1

            elif SC == "intergenic":
                if (
                    trans_id not in ref_trans
                    and gene_id not in ref_genes
                    and gene_id not in target_genes
                ):
                    target_trans.add(trans_id)
                    target_genes.add(gene_id)
                    counts[SC] -= 1

    # List of transcript and genes that will be deleted from reference
    # if all transcripts from a gene are being deleted the gene will be deleted too
    final_target = target_trans
    for gene in trans_by_gene:
        for trans in trans_by_gene[gene]:
            if trans[0] in target_trans:
                trans_by_gene[gene].remove(trans)
                if len(trans_by_gene[gene]) == 0:
                    final_target.add(gene)

    trans_index = pandas.read_csv(f_idx, sep="\t", header=0, dtype={"chrom":str})
    trans_index["sim_type"] = trans_index.apply(pick_sim_type, axis=1)
    trans_index["sim_type"] = trans_index["sim_type"].fillna("NA")
    trans_index.to_csv(
        f_idx_out, sep="\t", header=True, index=False, na_rep="NA"
    )

    return final_target


def getGeneID(line: str) -> str:
    """Returns the gene_id of a GTF line

    Args:
        line (str) line readed from GTF file

    Returns:
        gene_id (str) gene_id from that feature
    """

    line_split = line.split()
    gene_id = line_split[line_split.index("gene_id") + 1]
    gene_id = gene_id.replace(";", "").replace('"', "")

    return gene_id


def getTransID(line: str) -> str:
    """Returns the transcript_id of a GTF line

    Args:
        line (str) line readed from GTF file

    Returns:
        trans_id (str) transcript_id from that feature
    """

    try:
        line_split = line.split()
        trans_id = line_split[line_split.index("transcript_id") + 1]
        trans_id = trans_id.replace(";", "").replace('"', "")
    except:
        trans_id = None

    return trans_id


def modifyGTF(f_name_in: str, f_name_out: str, target: list):
    """
    Modify the original GTF deleting target transcripts to simulate specific
    SQANTI3 structural categorires

    Args:
        f_name_in (str) file name of the reference annotation GTF
        f_name_out (str) file name of the modified GTF generated
        target_trans (list) list of transcripts that will be deleted
    """

    f_out = open(f_name_out, "w")

    with open(f_name_in, "r") as gtf_in:
        for line in gtf_in:
            if line.startswith("#"):
                f_out.write(line)
            else:
                gene_id = getGeneID(line)
                trans_id = getTransID(line)
                if gene_id in target or trans_id in target:
                    pass
                else:
                    f_out.write(line)
    gtf_in.close()
    f_out.close()

    return


def simulate_gtf(args):
    """Generates the modified reference annotation"""

    print("[SQANTI-SIM] Writting modified GTF")
    counts = defaultdict(
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

    gtf_modif = os.path.join(args.dir, (args.output + "_modified.gtf"))
    f_idx_out = os.path.join(args.dir, (args.output + "_index.tsv"))

    target = target_trans(args.trans_index, f_idx_out, counts)
    modifyGTF(args.gtf, gtf_modif, target)

    return counts


def summary_table_del(counts_ini: dict, counts_end: dict):
    """Prints summary table of simulate_gtf

    Prints a summary table of the transcripts that were deleted from the
    original reference annotation file

    Args:
        counts_ini (dict) dictionary with all the transcripts associated to each
                          SQANTI3 structural category
        counts_end (dict) remaining classified transcripts after deletion
    """

    for sc in counts_end:
        counts_ini[sc] -= counts_end[sc]

    print("\033[94m_" * 79 + "\033[0m")
    print("\033[92mS Q A N T I - S I M\033[0m \U0001F4CA")
    print()
    print("Deleted transcripts from GTF \U0001F50E")
    print("\033[94m_" * 79 + "\033[0m")
    for k, v in counts_ini.items():
        print("\033[92m|\033[0m " + k + ": " + str(v))


def take_closest(my_list: list, number: int, bias: str)-> int:
    """Finds the clossest bottom value

    Given a list it finds the clossest value to the query interger always
    givin clossest LOWER or equal position, never a higher number (unless there
    is no lower value). Assumes list is sorted!

    Args:
        my_list (list) list with integers
        number (int) query integer
        bias (str) get higher "high" or lower "low" value
    
    Returns:
        pos (int) index of the clossest value
    """

    pos = bisect_left(my_list, number)
    if pos == 0:
        return pos
    if pos == len(my_list):
        return -1
    if bias == "high":
        return pos
    return pos-1 # return the value before

def create_expr_file_fixed_count(f_idx: str, args: list):
    """ Expression matrix - equal mode

    Modifies the index file adding the counts and TPM for the transcripts that
    will be simulated with a fixed count value

    Args:
        f_idx (str) index file name
        args (list) the number of transcripts and reads to be simulated
    """

    def fixed_coverage(row):
        if row["transcript_id"] in tot_trans:
            return coverage
        return 0

    novel_trans = []
    known_trans = []

    with open(f_idx, "r") as f_in:
        skip = f_in.readline()
        skip = skip.split()
        i = skip.index("sim_type")
        j = skip.index("transcript_id")
        for line in f_in:
            line = line.split()
            sim_type = line[i]
            if sim_type == "novel":
                novel_trans.append(line[j])
            else:
                known_trans.append(line[j])
    f_in.close()

    tot_trans = len(novel_trans) + len(known_trans)
    if args.trans_number > tot_trans:
        print(
            "[SQANTI-SIM] WARNING: A higher number than annotated transcripts was requested to simulate, only %s transcript will be simulated"
            % (tot_trans)
        )
        args.trans_number = tot_trans

    random.shuffle(known_trans)
    known_trans = known_trans[: (args.trans_number - len(novel_trans))]

    tot_trans = novel_trans + known_trans
    coverage = args.read_count // args.trans_number

    trans_index = pandas.read_csv(f_idx, sep="\t", header=0, dtype={"chrom":str})
    trans_index["requested_counts"] = trans_index.apply(fixed_coverage, axis=1)
    trans_index["requested_tpm"] = round(
        (
            (1000000.0 * trans_index["requested_counts"])
            / (trans_index["requested_counts"] * args.trans_number)
        ),
        2,
    ) # Not taking into account transcript length
    trans_index["requested_counts"] = trans_index["requested_counts"].fillna(0)
    trans_index["requested_tpm"] = trans_index["requested_tpm"].fillna(0)
    trans_index.to_csv(f_idx, sep="\t", header=True, index=False, na_rep="NA")


def create_expr_file_nbinom(f_idx: str, args: list):
    """ Expression matrix - custom mode

    Modifies the index file adding the counts and TPM for the transcripts that
    will be simulated from 2 different negative binomial distributions: one for
    known transcripts and the other for the novel ones

    Args:
        f_idx (str) index file name
        args (list) the number of transcripts to be simulated and parameters for
                    the negative binomial distributions
    """

    def nbinom_coverage(row):
        if row["transcript_id"] in novel_trans:
            coverage = nb_novel.pop()
        elif row["transcript_id"] in known_trans:
            coverage = nb_known.pop()
        else:
            coverage = 0
        return coverage

    novel_trans = []
    known_trans = []
    with open(f_idx, "r") as f_in:
        skip = f_in.readline()
        skip = skip.split()
        i = skip.index("sim_type")
        j = skip.index("transcript_id")
        for line in f_in:
            line = line.split()
            sim_type = line[i]
            if sim_type == "novel":
                novel_trans.append(line[j])
            else:
                known_trans.append(line[j])
    f_in.close()

    random.shuffle(known_trans)
    known_trans = known_trans[: (args.trans_number - len(novel_trans))]

    nb_known = numpy.random.negative_binomial(
        args.nbn_known, args.nbp_known, len(known_trans)
    ).tolist()
    nb_known = [
        1 if n == 0 else n for n in nb_known
    ]  # minimum one count per transcript
    nb_novel = numpy.random.negative_binomial(
        args.nbn_novel, args.nbp_novel, len(novel_trans)
    ).tolist()
    nb_novel = [
        1 if n == 0 else n for n in nb_novel
    ]  # minimum one count per transcript
    n_reads = sum(nb_known) + sum(nb_novel)

    trans_index = pandas.read_csv(f_idx, sep="\t", header=0, dtype={"chrom":str})
    trans_index["requested_counts"] = trans_index.apply(
        nbinom_coverage, axis=1
    )
    trans_index["requested_tpm"] = round(
        ((1000000.0 * trans_index["requested_counts"]) / n_reads), 2
    )
    trans_index["requested_counts"] = trans_index["requested_counts"].fillna(0)
    trans_index["requested_tpm"] = trans_index["requested_tpm"].fillna(0)
    trans_index.to_csv(f_idx, sep="\t", header=True, index=False, na_rep="NA")


def create_expr_file_sample(f_idx: str, args: list, tech: str):
    """ Expression matrix - sample mode

    Modifies the index file adding the counts and TPM for the transcripts that
    will be simulated using a real expression distribution

    Args:
        f_idx (str) index file name
        args (list) reference transcriptome and real reads
        tech (str) sequencing platform {pb, ont}
    """

    def sample_coverage(row):
        if row["transcript_id"] in novel_trans:
            coverage = novel_expr.pop()
        elif row["transcript_id"] in known_trans:
            coverage = known_expr.pop()
        else:
            coverage = 0
        return coverage

    # Extract fasta transcripts
    ref_t = os.path.join(os.path.dirname(args.genome), "sqanti_sim.transcripts.fa")
    if os.path.exists(ref_t):
        print("[SQANTI-SIM] WARNING: %s already exists. Overwritting!" %(ref_t))

    cmd = ["gffread", "-w", str(ref_t), "-g", str(args.genome), str(args.gtf)]
    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("[SQANTI-SIM] ERROR running gffread: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)

    # Read transcripts from index file
    novel_trans = []
    known_trans = []
    trans_to_gene = defaultdict(lambda: str())
    trans_by_gene = defaultdict(lambda: [])
    with open(f_idx, "r") as f_in:
        skip = f_in.readline()
        skip = skip.split()
        i = skip.index("sim_type")
        j = skip.index("transcript_id")
        k = skip.index("gene_id")
        for line in f_in:
            line = line.split()
            sim_type = line[i]
            if sim_type == "novel":
                novel_trans.append(line[j])
            else:
                known_trans.append(line[j])
            trans_to_gene[line[j]] = line[k]
            trans_by_gene[line[k]].append(line[j])
    f_in.close()

    # Align with minimap
    sam_file = "_".join(f_idx.split("_")[:-1]) + "_align_" + tech + ".sam"

    if tech == "pb":
        cmd = [
            "minimap2",
            ref_t,
            args.pb_reads,
            "-x",
            "map-pb",
            "-a",
            "--secondary=no",
            "-o",
            sam_file,
            "-t",
            str(args.cores),
        ]
    elif tech == "ont":
        cmd = [
            "minimap2",
            ref_t,
            args.ont_reads,
            "-x",
            "map-ont",
            "-a",
            "--secondary=no",
            "-o",
            sam_file,
            "-t",
            str(args.cores),
        ]

    cmd = " ".join(cmd)
    sys.stdout.flush()
    if subprocess.check_call(cmd, shell=True) != 0:
        print("[SQANTI-SIM] ERROR running minimap2: {0}".format(cmd), file=sys.stderr)
        sys.exit(1)

    # Raw counts -> Count only primary alignments
    trans_counts = defaultdict(lambda: 0)
    with pysam.AlignmentFile(sam_file, "r") as sam_file_in:
        for align in sam_file_in:
            trans_id = align.reference_name

            if (
                align.reference_id == -1
                or align.is_supplementary
                or align.is_secondary
            ):
                continue
            trans_counts[trans_id] += 1
    os.remove(sam_file)
    os.remove(ref_t)

    # Empirical expression values
    expr_distr = list(trans_counts.values())
    expr_distr.sort()

    if args.trans_number:
        n_trans = args.trans_number
    else:
        n_trans = len(expr_distr)

    if n_trans < len(novel_trans):
        print("[SQANTI-SIM] ERROR: -nt/--trans number must be higher than the novel transcripts to simulate")
        sys.exit(1)
    

    if n_trans > (len(novel_trans) + len(known_trans)):
        n_trans = (len(novel_trans) + len(known_trans))
        print(
            "[SQANTI-SIM] WARNING: A higher number than annotated transcripts was requested to simulate, only %s transcript will be simulated"
            % (n_trans)
        )
        expr_distr = expr_distr[-n_trans,]

    if args.iso_complex:
        gene_isoforms_counts = defaultdict(lambda: 0)
        for i in trans_counts:
            gene_id = trans_to_gene[i]
            if gene_id:
                gene_isoforms_counts[gene_id] += 1
        known_trans = []
        already_scaned = []
        complex_distr = list(gene_isoforms_counts.values())
        complex_distr.sort()
        for trans_id in novel_trans:
            gene_id = trans_to_gene[trans_id]
            if gene_id in already_scaned:
                continue
            n = len(trans_by_gene[gene_id])
            pos = take_closest(complex_distr, n, "low")
            diff_isos = complex_distr.pop(pos)
            if diff_isos > n:
                diff_isos = n
            for i in range(diff_isos):
                curr_trans = trans_by_gene[gene_id][i]
                if trans_id in novel_trans:
                    continue
                else:
                    known_trans.append(curr_trans)
            del trans_by_gene[gene_id]
            already_scaned.append(gene_id)
            if len(complex_distr) <= 0:
                complex_distr = list(gene_isoforms_counts.values())

        counts_to_gene = defaultdict(lambda: [])
        for i in trans_by_gene:
            counts_to_gene[len(trans_by_gene[i])].append(i)
    
        random.shuffle(complex_distr)
        while n_trans > (len(novel_trans) + len(known_trans)) and len(complex_distr) > 0:
            diff_isos = complex_distr.pop()
            pos = take_closest(list(counts_to_gene.keys()), diff_isos, "high")
            if counts_to_gene[pos]:
                gene_id = counts_to_gene[pos].pop()
                if pos > diff_isos:
                    for i in range(diff_isos):
                        curr_trans = trans_by_gene[gene_id][i]
                        known_trans.append(curr_trans)
                else:
                    for i in range(pos):
                        curr_trans = trans_by_gene[gene_id][i]
                        known_trans.append(curr_trans)
                del trans_by_gene[gene_id]
            
                if len(counts_to_gene[pos]) == 0:
                    del counts_to_gene[pos]
            if len(complex_distr) <= 0:
                complex_distr = list(gene_isoforms_counts.values())
                random.shuffle(complex_distr)

        if n_trans < (len(novel_trans) + len(known_trans)):
            known_trans[:n_trans]

    else:
        random.shuffle(known_trans)
        known_trans = known_trans[: (n_trans - len(novel_trans))]

    trans_index = pandas.read_csv(f_idx, sep="\t", header=0, dtype={"chrom":str})
    if args.diff_exp:
        # Generate a vector of inverse probabilities to assign lower TPM values to novel transcripts and higher to known transcripts
        prob = numpy.linspace(
            start=args.low_prob, stop=args.high_prob, num=int((max(expr_distr)+1) - min(expr_distr))
        )

        # Choose expression values for novel and known transcripts
        min_expr = min(expr_distr)
        novel_expr = []
        n_novel = 0
        while n_novel < len(novel_trans):
            s = random.choice(expr_distr)
            r = random.uniform(0, 1)
            if r > prob[(s - min_expr)]:
                novel_expr.append(s)
                n_novel += 1

        known_expr = []
        n_known = 0
        while n_known < len(known_trans):
            s = random.choice(expr_distr)
            r = random.uniform(0, 1)
            if r < prob[(s - min_expr)]:
                known_expr.append(s)
                n_known += 1

        trans_index["requested_counts"] = trans_index.apply(
            sample_coverage, axis=1
        )

    else: # No bias for novel/known expression distribution
        novel_expr = []
        n_novel = 0
        while n_novel < len(novel_trans):
            s = random.choice(expr_distr)
            novel_expr.append(s)
            n_novel += 1

        known_expr = []
        n_known = 0
        while n_known < len(known_trans):
            s = random.choice(expr_distr)
            known_expr.append(s)
            n_known += 1

        trans_index["requested_counts"] = trans_index.apply(
            sample_coverage, axis=1
        )

        
    n_reads = trans_index["requested_counts"].sum()
    trans_index["requested_tpm"] = round(
        ((1000000.0 * trans_index["requested_counts"]) / n_reads), 2
    )
    trans_index["requested_counts"] = trans_index["requested_counts"].fillna(0)
    trans_index["requested_tpm"] = trans_index["requested_tpm"].fillna(0)
    trans_index.to_csv(f_idx, sep="\t", header=True, index=False, na_rep="NA")

    print("[SQANTI-SIM] Mapped transcripts: %s" %(len(expr_distr)))
    print("[SQANTI-SIM] Mapped reads: %s" %(n_reads))

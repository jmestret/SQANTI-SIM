### Table of Contents

- [Introduction](#intro)
- [Usage](#use)
- [Arguments](#args)
- [Output explanation](#out)

## <a name="intro"></a>Introduction

`sqantisim.py sim` is a wrapper tool that simulates in a controlled way the reads for the experiment you designed with the SQANTISIM *design* mode. This mode takes as input the index file with the expression values, the original reference annotation in GTF format and the reference genome in FASTA format. It will simulate ONT reads ([NanoSim](https://github.com/bcgsc/NanoSim) pipeline), PacBio reads ([IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) pipeline) and/or Illumina reads ([Polyester](https://github.com/alyssafrazee/polyester) tool).

In the *sim* step, you can simulate in different sequencing depths changing the `--long_count` and `--short_count` parameters. If not used, it will use the values of the *requested_counts* from the index file assigned in the *design* step. NanoSim and IsoSeqSim pipelines have been slightly adapted so we can control exactly how many reads are simulated from each transcript.

To simulate ONT transcriptome reads (cDNA / dRNA) we use the NanoSim simulator. It is run using the default NanoSim options, with the basecaller guppy and not modelling intron retention. You can decide with the parameter `--read_type` whether to simulate cDNA or dRNA (default option). The error models used are *human_NA12878_dRNA_Bham1_guppy* and *human_NA12878_cDNA_Bham1_guppy*, taken from the pre-trained models by NanoSim The output reads are given in FASTQ format.

To simulate PacBio Full-length reads, we use the IsoSeqSim simulator. It is run in the normal mode with the 5' end and 3' end completeness of PacBio Sequel and the error rates of PacBio Sequel suggested in the IsoSeqSim GitHub (substitution 1.731%, deletion 1.090% and insertion 2.204%). The output reads are given in FASTA format.

Finally, to simulate Illumina RNA-seq reads, we use the Polyester simulator. It is run using the default options, generating paired-end read of 100 nucleotides length. The error model is uniformly distributed at an error rate of 0.5%.

The output of this tool is the index file with the number of reads simulated for each transcript and the FASTA or FASTQ files with the simulated reads. It also generates a file assigning the simulated read id to the reference transcript it comes from.

**IMPORTANT**. After this step and before running SQANTISIM *eval* you should use your transcript reconstruction pipeline to generate the transcriptome in GTF format. You must input your pipeline the simulated reads generated in this step and the modified reference annotation generated in the *design* step.

## <a name="use"></a>Usage

SQANTISIM *sim* mode usage:

```
sqantisim.py sim [-h] --gtf GTF --genome GENOME -i TRANS_INDEX
                         [--read_type READ_TYPE] [-d DIR] [-k CORES]
                         (--pb | --ont) [--illumina] [--long_count LONG_COUNT]
                         [--short_count SHORT_COUNT] [-s SEED]
```

With the `--help` option you can display a complete description of the arguments:

```
sqantisim.py sim parse options

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Reference annotation in GTF format
  --genome GENOME       Reference genome FASTA
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTISIM
  --read_type READ_TYPE
                        Read type for NanoSim simulation
  -d DIR, --dir DIR     Directory for output files (default: .)
  -k CORES, --cores CORES
                        Number of cores to run in parallel
  --pb                  If used the program will simulate PacBio reads with
                        IsoSeqSim
  --ont                 If used the program will simulate ONT reads with
                        NanoSim
  --illumina            If used the program will simulate Illumina reads with
                        RSEM
  --long_count LONG_COUNT
                        Number of long reads to simulate (if not given it will
                        use the counts of the given expression file)
  --short_count SHORT_COUNT
                        Number of short reads to simulate (if not given it
                        will use the counts of the given expression file)
  -s SEED, --seed SEED  Randomizer seed
```

Running the tool with the minimum input to simulate PacBio reads will look as follows:

```
(SQANTISIM.env)$ python sqantisim.py sim \
			--trans_index prefix_index.tsv \
			--gtf reference_annotation.gtf \
			--genome reference_genome.fasta \
			--pb
```

## <a name="args"></a>Arguments detailed explanation

### Required input

These are the minimum parameters you will need to run `sqantisim.py classif`:

- **Transcript index** file (`-i`): This file is the *prefix_index.tsv* file generated in the previous *design* step. New columns and information will be added to this file, so you should track it and use the most updated version of this file in the next steps of the SQANTISIM pipeline.
- **Reference annotation** in GTF format (`--gtf`): This file must be the same as the one used in the *classif* and *design* steps, so the transcript names and references match properly. An example of a reference transcriptome and its required format is [GENCODE](https://www.gencodegenes.org/).
- **Reference genome** in FASTA format (`--genome`): This is the reference genome in FASTA format.
- **Long-read sequencing platform** (`--pb/--ont`): Decide whether to simulate PacBio reads using the IsoSeqSim simulator or ONT reads with the NanoSim simulator. You can choose how many reads simulate with `--long_count`; if not used, it will use the values from the "requested_counts" column from the index file. If `--long_count` is used, the value of total simulate reads will not be exact; the number of total reads is used to multiply the requested TPM to decide how many reads simulate for each different transcript.

### Optional input

- **Illumina reads** simulation (`--illumina`): If used, you will also simulate Illumina RNA-seq reads with the Polyester simulator. You can choose how many reads simulate with `--short_count`; if not used, it will use the values from the "requested_counts" column from the index file.
- **Output directory**(`-d`): Output directory for output files. SQANTISIM will use the directory where the script was run as the default output directory.
- **Parallelization**(`-k`): Number of cores to run in parallel. Most SQANTISIM modes have code chunks that can be run in parallel. However the default option is to run SQANTISIM in one single thread.

## <a name="out"></a>Output explanation

#### ONT_simulated.fastq or PacBio_simulated.fasta

FASTQ or FASTA file of simulated long-reads. The header of the reads looks like the following:

```
# ONT reads
@ENST00000686945_ONT_simulated_read_1229

# PacBio reads
>ENST00000379236_PacBio_simulated_read_1229
```

The information before the first `_` in the reference transcript name from which the read comes. It is written the platform from which you are simulating (`ONT_simulated_read` or `PacBio_simulated_read`) and the total read count of simulated read, meaning this read was the simulated read number `1230` (the count starts in 0).

#### {ONT,PacBio}_simulated.read_to_isoform.tsv

This is a tab-separated file with two columns and no header. In the first column, we have all the read names, and in the second column, we have the reference transcript name from which that specific read was simulated. This information is redundant because we can find the reference transcript in the header of the read.

#### Illumina_simulated_{1,2}.fasta

Fasta file of simulated Illumina RNA-seq reads. The header format is the same as the one given by the Polyester simulator.

#### prefix_index.tsv

This is the main index file with the structural annotation of the reference transcripts. This file will be requested and modified in each other SQANTISIM mode. It must be parsed using `-i/--trans_index`. After running the SQANTISIM *classif* mode, this file will contain the potential structural category of each transcript. The columns of the *prefix_index.tsv* file are described below:

1. **transcript_id**: The transcript id of the query transcript taken from the reference annotation.
2. **gene_id**: The reference gene where this transcript comes from.
3. **structural_category**: The potential structural category in which this transcript could be classified if simulated as a novel.
4. **associated_gene**: The reference gene that confers the query transcript's structural category. If there is more than one associated gene, they will be written separated by "_", e.g., *Gene1_Gene2_Gene3*.
5. **associated_trans**: The reference transcript name that confers the structural category. If there is no specific associated transcript, "novel" will appear instead.
6. **chrom**: Chromosome.
7. **strand**: Strand.
8. **exons**: Number of exons.
9. **donors**: End of the junction (if - strand, it will be the start of the junction). Values are 1-based.
10. **acceptors**: Start of the junction (if - strand, it will be the end of the junction). Values are 1-based.
11. **TSS_genomic_coord**: Start site of the transcript. Values are 1-based.
12. **TTS_genomic_coord**: Termination site of the transcript. Values are 1-based.
13. **sim_type**: If the transcript is present in the modified reference annotation (*known*) or not (*novel*).
14. **requested_counts**: The requested reads to simulate in the *sim* step.
15. **requested_tpm**: The requested TPM value to simulate in the *sim* step.
16. **sim_counts**: The number of simulated long-reads.
17. **illumina_counts**: The number of simulated short-reads.

#### Other output generated by default by the simulators

- **ONT_simulated_aligned_error_profile**: Contains all the information of error introduced into each read.
- **ONT_simulated_aligned_reads.fastq**: Contains all aligned reads. In SQANTISIM, all simulated reads.
- **ONT_simulated_unaligned_reads.fastq**: Contains all unaligned reads. NanoSim was adapted in SQANTISIM, so all the simulated reads are aligned reads.
- **PacBio_simulated.tsv**: Structural annotation of each reference transcript and the number of reads simulated.
- **temp_isoseqsim**: Temporary directory for IsoSeqSim simulation.


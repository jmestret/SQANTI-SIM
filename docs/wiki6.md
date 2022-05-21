### Table of Contents

- [Introduction](#intro)
- [Usage](#use)
- [Arguments](#args)
- [Output explanation](#out)

## <a name="intro"></a>Introduction

`sqantisim.py design` creates the modified GTF to simulate the novel transcripts and generates the expression matrix for the transcripts simulated with the SQANTISIM *sim* mode. This tool inputs the reference annotation in GTF format that was used in the *classif* step and the transcript index file (*prefix_index.tsv*) generated. In this step, you will decide which SQANTI3 structural categories you want to simulate as novel and how many transcripts (it will simulate transcripts with a minimum length of 200bp).

Moreover, you will decide the expression values (counts and TPM) for those simulated transcripts. You can do this in 3 different modes:

- **equal**: This mode will assign the same TPM and count values for all the simulated transcripts, not considering whether they are novel or known.

- **custom**: With this mode, you can customise the differences in expression between the novel and known transcripts. For so, you can choose the parameters of 2 negative binomial distributions that will be assigned, one for the novel transcripts and the other for the known. Zero values samples from those distributions will be swapped by one, so all requested transcripts will always have at least the minimum expression value.

- **sample**: Finally, you can reproduce the expression values from a real sample. In this mode, a real sample is mapped using Minimap2 and the raw count distribution is used to assign the expression values to the requested transcripts to simulate. This mode, by default, gives the expression values to novel and known transcripts without bias by simply bootstrapping the empirical distribution. However, using the `--diff_exp` parameter, SQANTISIM will generate an inverted vector of probabilities to assign lower expression values to novel transcripts and higher to known. You can customise this probability vector by selecting the start (`--low_prob`) and end (`--high_prob`) probabilities of the range.

![design_modes](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/design_modes.png)

The output of this tool will be the modified reference annotation in GTF format and the index file with new columns referring to the expression values. **IMPORTANT**: The modified reference annotation must be used as reference when using your transcript identification pipeline.

## <a name="use"></a>Usage

SQANTISIM *design* mode usage:

```
python sqantisim.py design [-h] {equal,custom,sample} ...
```

With the `--help` option you can display a complete description of the arguments:

```
(SQANTISIM.env)$ python sqantisim.py design --help

sqantisim.py design parse options

optional arguments:
  -h, --help            show this help message and exit

subcommands:
  Different modes to generate the expression matrix: equal (simulate with
  equal coverage for all reads), custom (simulate with diferent negative
  binomial distributions for novel and known transcripts) or sample
  (simulate using a real sample)

  {equal,custom,sample}
    equal               Run in equal mode
    custom              Run in custom mode
    sample              Run in sample mode
```

#### equal mode

```
(SQANTISIM.env)$ python sqantisim.py design equal --help

usage: sqantisim.py design equal [-h] -i TRANS_INDEX --gtf GTF
                                       [-o OUTPUT] [-d DIR] [-nt TRANS_NUMBER]
                                       [--read_count READ_COUNT] [--ISM ISM]
                                       [--NIC NIC] [--NNC NNC]
                                       [--Fusion FUSION]
                                       [--Antisense ANTISENSE] [--GG GG]
                                       [--GI GI] [--Intergenic INTERGENIC]
                                       [-k CORES] [-s SEED]

optional arguments:
  -h, --help            show this help message and exit
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTISIM
  --gtf GTF             Reference annotation in GTF format
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  -nt TRANS_NUMBER, --trans_number TRANS_NUMBER
                        Number of different transcripts to simulate
  --read_count READ_COUNT
                        Number of reads to simulate
  --ISM ISM             Number of incomplete-splice-matches to delete
  --NIC NIC             Number of novel-in-catalog to delete
  --NNC NNC             Number of novel-not-in-catalog to delete
  --Fusion FUSION       Number of Fusion to delete
  --Antisense ANTISENSE
                        Number of Antisense to delete
  --GG GG               Number of Genic-genomic to delete
  --GI GI               Number of Genic-intron to delete
  --Intergenic INTERGENIC
                        Number of Intergenic to delete
  -k CORES, --cores CORES
                        Number of cores to run in parallel
  -s SEED, --seed SEED  Randomizer seed
```

Running the *equal* mode with the minimum input will look as follows:

```
(SQANTISIM.env)$ python sqantisim.py design equal \
			--gtf reference_annotation.gtf \
			--trans_index prefix_index.tsv \
			--read_count 50000
```

#### custom mode

```
(SQANTISIM.env)$ python sqantisim.py design custom --help

usage: sqantisim.py design custom [-h] -i TRANS_INDEX --gtf GTF
                                        [-o OUTPUT] [-d DIR]
                                        [-nt TRANS_NUMBER]
                                        [--nbn_known NBN_KNOWN]
                                        [--nbp_known NBP_KNOWN]
                                        [--nbn_novel NBN_NOVEL]
                                        [--nbp_novel NBP_NOVEL] [--ISM ISM]
                                        [--NIC NIC] [--NNC NNC]
                                        [--Fusion FUSION]
                                        [--Antisense ANTISENSE] [--GG GG]
                                        [--GI GI] [--Intergenic INTERGENIC]
                                        [-k CORES] [-s SEED]

optional arguments:
  -h, --help            show this help message and exit
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTISIM
  --gtf GTF             Reference annotation in GTF format
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  -nt TRANS_NUMBER, --trans_number TRANS_NUMBER
                        Number of different transcripts to simulate
  --nbn_known NBN_KNOWN
                        Average read count per known transcript to simulate
                        (i.e., the parameter 'n' of the Negative Binomial
                        distribution)
  --nbp_known NBP_KNOWN
                        The parameter 'p' of the Negative Binomial
                        distribution for known transcripts
  --nbn_novel NBN_NOVEL
                        Average read count per novel transcript to simulate
                        (i.e., the parameter 'n' of the Negative Binomial
                        distribution)
  --nbp_novel NBP_NOVEL
                        The parameter 'p' of the Negative Binomial
                        distribution for novel transcripts
  --ISM ISM             Number of incomplete-splice-matches to delete
  --NIC NIC             Number of novel-in-catalog to delete
  --NNC NNC             Number of novel-not-in-catalog to delete
  --Fusion FUSION       Number of Fusion to delete
  --Antisense ANTISENSE
                        Number of Antisense to delete
  --GG GG               Number of Genic-genomic to delete
  --GI GI               Number of Genic-intron to delete
  --Intergenic INTERGENIC
                        Number of Intergenic to delete
  -k CORES, --cores CORES
                        Number of cores to run in parallel
  -s SEED, --seed SEED  Randomizer seed
```

Running the *custom* mode with the minimum input will look as follows:

```
(SQANTISIM.env)$ python sqantisim.py design custom \
			--gtf reference_annotation.gtf \
			--trans_index prefix_index.tsv \
			--nbn_known 15 --nbp_known 0.5 \
			--nbn_novel 5 --nbp_novel 0.5
```

#### sample mode

```
(SQANTISIM.env)$ python sqantisim.py design sample --help

usage: sqantisim.py design sample [-h] -i TRANS_INDEX --gtf GTF
                                        [-o OUTPUT] [-d DIR]
                                        [-nt TRANS_NUMBER] --genome GENOME
                                        [--pb_reads PB_READS | --ont_reads ONT_READS]
                                        [--diff_exp] [--low_prob LOW_PROB]
                                        [--high_prob HIGH_PROB] [--ISM ISM]
                                        [--NIC NIC] [--NNC NNC]
                                        [--Fusion FUSION]
                                        [--Antisense ANTISENSE] [--GG GG]
                                        [--GI GI] [--Intergenic INTERGENIC]
                                        [-k CORES] [-s SEED]

optional arguments:
  -h, --help            show this help message and exit
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTISIM
  --gtf GTF             Reference annotation in GTF format
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  -nt TRANS_NUMBER, --trans_number TRANS_NUMBER
                        Number of different transcripts to simulate
  --genome GENOME       Reference genome FASTA
  --pb_reads PB_READS   Input PacBio reads for quantification
  --ont_reads ONT_READS
                        Input ONT reads for quantification
  --diff_exp            If used the program will simulate different expression
                        values for novel and known transcripts
  --low_prob LOW_PROB   Low value of prob vector (if --diff_exp)
  --high_prob HIGH_PROB
                        High value of prob vector (if --diff_exp)
  --ISM ISM             Number of incomplete-splice-matches to delete
  --NIC NIC             Number of novel-in-catalog to delete
  --NNC NNC             Number of novel-not-in-catalog to delete
  --Fusion FUSION       Number of Fusion to delete
  --Antisense ANTISENSE
                        Number of Antisense to delete
  --GG GG               Number of Genic-genomic to delete
  --GI GI               Number of Genic-intron to delete
  --Intergenic INTERGENIC
                        Number of Intergenic to delete
  -k CORES, --cores CORES
                        Number of cores to run in parallel
  -s SEED, --seed SEED  Randomizer seed
```

Running the *sample* mode with the minimum input and mode-specific parameters will look as follows:

```
(SQANTISIM.env)$ python sqantisim.py design sample \
			--gtf reference_annotation.gtf \
			--trans_index prefix_index.tsv \
			--genome reference_genome.fasta \
			--pb_reads reads.fastq
```

## <a name="args"></a>Arguments detailed explanation

### Required input

These are the minimum parameters you will need to run `sqantisim.py design`:

- **Design mode**: This is the first position argument right after "design" (`sqantisim.py design <mode>`). It can take three different values: *equal* (simulate with equal coverage for all reads), *custom* (simulate with different negative binomial distributions for novel and known transcripts) or *sample* (simulate using a real sample)
- **Transcript index** file (`-i`): This file is the *prefix_index.tsv* file generated in the previous *classif* step. New columns and information will be added to this file, so you should track it and use the most updated version of this file in the next steps of the SQANTISIM pipeline.
- **Reference annotation** in GTF format (`--gtf`): This file must be the same as the one used in the *classif* step, so the transcript names and references match correctly. From this file, SQANTISIM will generate a modified reference annotation to use as a reference in your transcript reconstruction pipeline. An example of a reference transcriptome and it required format is [GENCODE](https://www.gencodegenes.org/).
- **Reference genome** in FASTA format (`--genome`): Only required for *sample* mode, if not ignored. This is the reference genome in FASTA format. It will be used together with the reference annotation to extract the transcript sequences in FASTA format using GffRead.
- **Long-read mRNA sequences** (`--pb_reads/--ont_reads`): Only required for *sample* mode, if not ignore. A real sample of Long-read mRNA sequences in FASTA or FASTQ format to quantify and extract the expression values. The reads will be mapped with Minimap2, and the expression distribution will be generated from the raw counts of the alignment.

### Optional input

- **Novel transcripts to simulate**: The number of transcripts of each structural category to simulate as a novel. Use `--ISM` to simulate *Incomplete Splice Match* transcripts, `--NIC` for *Novel In Catalog*, `--NNC` for *Novel not In Catalog*, `--Fusion` for *Fusion*, `--Antisense` for *Antisense*, `--GG` for *Genig-genomic*, `--GI` for *Genic-intron* and `--Intergenic` for *Intergenic* transcripts. Use this parameter followed by the number of transcripts to simulate for that structural category.
- **Transcripts to simulate** (`-nt`): This is the total number to simulate; this means that this will be all the different transcripts that will be simulated. Reads will not come from a more significant number of transcripts than this one. The number of novel transcripts (those assigned to other structural categories) will be subtracted from this value. The rest will be simulated from known annotated transcripts kept in the modified reference annotation. Then, the number of different known transcripts to simulate will be the requested number of transcripts to simulate (`-nt`) minus the number of novel transcripts (sum of different structural categories)
- **Output prefix** (`-o`): The output prefix for the index file. SQANTISIM will use "sqantisim" as the default prefix.
- **Output directory** (`-d`): Output directory for output files. SQANTISIM will use the directory where the script was run as the default output directory.
- **Parallelization** (`-k`): Number of cores to run in parallel. Most SQANTISIM modes have code chunks that can be run in parallel. However the default option is to run SQANTISIM in one single thread.
- **Number of reads** (`--read_count`): Can only be used for *equal* mode, if not ignored. The number of total reads to simulate. It will be divided by the number of total transcripts to compute the requested counts for each transcript and the TPM values. All transcripts will have the same TPM and count values.
- **Isoform complexity** per gene (`--iso_complex`): Can only be used in *sample* mode. By default, SQANTISIM chooses the known transcripts randomly to simulate; however, if this argument is used, it will mimic the isoform expression complexity of the given sample, i.e., it will simulate transcripts so that the differently expressed isoforms per gene follow the same distribution as the real sample (`--pb_reads/--ont_reads`).
- **Negative binomial distribution** parameters: Can only be used in *custom* mode, if not ignored. When running in *custom* modes, expression values are taken from 2 different negative binomial distributions, one for novel transcripts and the other for known so that you can customise the differences in expression levels. To generate both negative binomial distributions you need to specify the parameter *n* successes (`--nbn_known` and `--nbn_novel`) and *p* probability of success (`--nbp_known` and `--nbp_novel`) where *n* is > 0 and *p* is in the interval [0,1].
- **Different expression** novel and known (`--diff_exp`): Can only be used in *sample* mode, if not ignored. *sample* mode by default simulates the same level of expression for novel and known transcripts. However, this parameter can be used to simulate different levels of expressions. It uses a probability vector in the interval [0,1] and its inverse to assign probability values for each count value of the expression distribution. You can decide the start and end of the range of the probability vector using `--low_prob` as the start probability and `--high_prob` as the end probability.

## <a name="out"></a>Output explanation

#### prefix_modified.gtf

This file contains the reference annotation in GTF format without those transcripts classified as novel in the index file. This file has the same content as the original reference annotation except for those novel transcripts and those genes left without a single annotated transcript when removing the novel transcripts. It is **mandatory** to use this file as the reference annotation in your transcript identification pipeline to ensure that the novel transcripts that will be simulated in the *sim* step are supposed to be novel for your pipeline.

#### prefix_index.tsv

This is the main index file with the structural annotation of the reference transcripts. This file will be requested and modified in each other SQANTISIM mode. It must be parsed using `-i/--trans_index`. After running the SQANTISIM *design* mode, this file will contain three extra columns referring to if it is novel or known and the transcript expression value. The columns of the *prefix_index.tsv* file are described below:

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


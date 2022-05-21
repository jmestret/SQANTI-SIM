### Table of Contents

- [Introduction](#intro)
- [Usage](#use)
- [Arguments](#args)
- [Output explanation](#out)

## <a name="intro"></a>Introduction

`sqantisim.py classif` allows you to classify all transcripts from a GTF reference annotation into its potential SQANTI3 structural category. Structural categories are assigned assuming that the query transcript cannot match itself as if it was a **novel** not annotated transcript.

![classif_example](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/classif_example.png)

Transcripts are classified using the same structural category decision tree as SQANTI3, so it has a **consistent transcript structural classification** with the different SQANTI3 tools. The output of this tool is a tab-separated index file with some structural annotations of all the transcripts from the reference annotation. It will also print a summary table with all the transcripts classified in each structural category.

![decission_tree](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/sqantisim_class_decision_tree.png)

The main purpose of this tool is to be used as part of the SQANTISIM pipeline. However, it can also be helpful for other interests, e.g., the *classif* mode can be used to identify all those redundant transcripts in the annotation that match as FSM with other annotated transcripts, meaning they have the same splice junctions.

## <a name="use"></a>Usage

SQANTISIM *classif* mode usage:

```
python sqantisim.py classif [-h] --gtf GTF [-o OUTPUT] [-d DIR] [-k CORES]
```

With the `--help` option you can display a complete description of the arguments:

```
(SQANTISIM.env)$ python sqantisim.py classif --help

sqantisim.py classif parse options

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Reference annotation in GTF format
  -o OUTPUT, --output OUTPUT
                        Prefix for output file
  -d DIR, --dir DIR     Directory for output files (default: .)
  -k CORES, --cores CORES
                        Number of cores to run in parallel
```

Running the tool with the minimum input will look as follows:

```
(SQANTISIM.env)$ python sqantisim.py classif --gtf reference_annotation.gtf
```

## <a name="args"></a>Arguments detailed explanation

### Required input

These are the minimum parameters you will need to run `sqantisim.py classif`:

- **Reference annotation** in GTF format (`--gtf`): Annotated transcripts from this file will be taken as input as will be classified in SQANTI3 structural categories. An example of reference transcriptome and it required format is [GENCODE](https://www.gencodegenes.org/).

### Optional input

- **Output prefix**(`-o`): The output prefix for the index file. SQANTISIM will use "sqantisim" as the default prefix.
- **Output directory**(`-d`): Output directory for output files. SQANTISIM will use the directory where the script was run as the default output directory.
- **Parallelization**(`-k`): Number of cores to run in parallel. Most SQANTISIM modes have code chunks that can be run in parallel. However the default option is to run SQANTISIM in one single thread.

## <a name="out"></a>Output explanation

#### prefix_index.tsv

This is the main index file with the structural annotation of the reference transcripts. This file will be requested and modified in each other SQANTISIM mode. It must be parsed using `-i/--trans_index`. After running the SQANTISIM *classif* mode this file will contain the potential structural category of each transcript. The columns of the *prefix_index.tsv* file are described below:

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


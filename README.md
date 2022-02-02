# SQANTI-SIM

**SQANTI-SIM** is a simulator of controlled novelty and degradation of transcripts sequenced by long reads. It is a plugin for the SQANTI3 tool ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI3)).
SQANTI-SIM is a wrapper tool for RNA-Seq long-reads simulators such as [IsoSeqSim]() or [NanoSim]()(formerly Trans-NanoSim).
Currently this simulators don't focus in simulate the structural categories of the transcripts sequenced, neither aiming to simulate the novelty and degradation with a reliable ground-truth.
The aim of SQANTI-SIM is to simulate novel and degradated transcripts in a controlled way using as ground-truth real transcripts in the reference that fit the description of the different SQANTI3 structural categories.

## Requirements and Installation

- Python 3
- Python packages: argparse, copy, subprocess and tqdm

## Running SQANTI-SIM

SQANTI-SIM can be used in two different ways: (i) to classify all transcripts from a GTF file and then modifiy the original annotation GTF file or (ii) give as input the intermediary file with the structural categories and just modify the original GTF for different porpuses.

### Required input

You have to give either one of the next options. You can only give one of both deppending on your interest:

- **--gtf** expects a GTF reference annotation file to analyze and then generate a modified file from this one.
- **--cat** expects the intermediate file with the structural categories of the transcripts from the reference GTF. This argument is used when you have already classified the transcripts from a GTF and you just want to generate a new modified GTF from the same reference.

### Optional input

- **-o** is the prefix for the output files (default = sqanti_sim)
- **-d** is the path to the directory to save the output files (default = .)
- You can choose how many transcripts from each structural category you want to delete from the reference using **--ISM**, **--NIC** and **--NNC** (default = 0).
- TODO: **-k** to choose the number of cores to run simultaneusly in the structural categorization
- **-v** to show up the version of the tool 

```
#!bash
$ python sqanti3_sim.py -h

usage: sqanti3_sim.py [-h] [--gtf GTF | --cat CAT] [-o OUTPUT] [-d DIR] [--ISM ISM] [--NIC NIC] [--NNC NNC] [-k CORES] [-v]

SQANTI-SIM: a simulator of controlled novelty and degradation of transcripts sequence by long-reads

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Reference annotation in GTF format
  --cat CAT             File with transcripts structural categories generated with SQANTI-SIM
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files. Default: Directory where the script was run
  --ISM ISM             Number of incomplete-splice-matches to delete
  --NIC NIC             Number of novel-in-catalog to delete
  --NNC NNC             Number of novel-not-in-catalog to delete
  -k CORES, --cores CORES
                        Number of cores to run in parallel
  -v, --version         Display program version number
```

## Example run

Example data can be found in `data/example_data.tar.gz`. Unpack it using `tar -xzf example_data.tar.gz` in `data` folder.

### Using a new reference GTF as input

`python sqanti3_sim.py --gtf data/example_ref.gtf -o test -d data --ISM 100 --NIC 200 --NNC 100`

### Using an existing SC classification file
`python sqanti3_sim.py --cat data/test_categories.txt -o test -d data --ISM 200 --NIC 100 --NNC 200`

## Output explanation

This tool generates mainly 2 output files:

- The modified GTF file used to "simulate" using it as reference annotation file, with the suffix **_modified.gtf**
- An intermediate file with the suffix **_categories.txt** that contains the classification in SQANTI3 structural categories for each transcript o the reference GTF (this is the ground-truth)

```
TransID	GeneID	SC	RefGene	RefTrans
ENST00000661675.1	ENSG00000223587.2	NIC	ENSG00000223587.2	novel
ENST00000656877.1	ENSG00000223587.2	NNC	ENSG00000223587.2	novel
ENST00000660204.1	ENSG00000223587.2	NIC	ENSG00000223587.2	novel
ENST00000660491.1	ENSG00000223587.2	NNC	ENSG00000223587.2	novel
ENST00000659825.1	ENSG00000223587.2	NIC	ENSG00000223587.2	novel
ENST00000440867.1	ENSG00000223587.2	ISM	ENSG00000223587.2	ENST00000660204.1
ENST00000426697.1	ENSG00000224918.1	Intergenic	NA	NA
ENST00000663345.1	ENSG00000224318.6	Fusion	ENSG00000224318.6_ENSG00000252017.1_ENSG00000231660.1	NA
ENST00000657108.1	ENSG00000224318.6	Fusion	ENSG00000224318.6_ENSG00000252017.1_ENSG00000231660.1	NA
```

Moreover, in the terminal it shows a sumary table with all the transcripts classified for each structural category:

```
_______________________________________________________________________________
S Q A N T I - S I M 📊

Summary Table 🔎
_______________________________________________________________________________
| FSM: 278
| ISM: 2473
| NIC: 1587
| NNC: 907
| Fusion: 7112
| Antisense: 413
| Genic-genomic: 128
| Genic-intron: 350
| Intergenic: 1585
| Unclassified: 0
```


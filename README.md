![SQANTI-SIM logo](https://github.com/jorgemt98/SQANTI-SIM/blob/main/sqantisim_logo.png)

# SQANTI-SIM

**SQANTI-SIM** is a simulator of controlled novelty and degradation of transcripts sequenced by long reads. It is a plugin for the SQANTI3 tool ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI3)).
SQANTI-SIM is a wrapper tool for RNA-Seq long-reads simulators such as [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) or [NanoSim](https://github.com/bcgsc/NanoSim)(formerly Trans-NanoSim).
Currently this simulators don't focus in simulating the structural categories of the transcripts sequenced, neither aiming to simulate the novelty and degradation with a reliable ground-truth.
The aim of SQANTI-SIM is to simulate novel and degradated transcripts in a controlled way using as ground-truth real transcripts in the reference annotation that fit the description of the different SQANTI3 structural categories.

## Requirements and Installation

As SQANTI-SIM is a plugin for the real SQANTI3 programm, in order to use it is necessary to install the SQANTI3 dependencies (or anaconda enviroment) as described in its [wiki](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-dependencies-and-installation). Moreover, it is needed to intall the python package tqdm. You can installit in your SQANTI3 conda enviroment as follows:

```
conda activate SQANTI3.env

conda install -c conda-forge tqdm 
```

## Running SQANTI-SIM

SQANTI-SIM can be used in two different ways: (i) to classify all transcripts from a GTF file and then modifiy the original annotation GTF file or (ii) give as input the intermediary file with the structural categories and just modify the original GTF for different porpuses.

![SQANTI-SIM workflow](https://github.com/jorgemt98/SQANTI-SIM/blob/main/sqantisim_workflow.png)

### Required input

- **--gtf** expects a GTF reference annotation file to analyze and then generate a modified file from this one.

### Optional input

- **--cat** expects the intermediate file with the structural categories of the transcripts from the reference GTF. This argument is used when you have already classified the transcripts from a GTF and you just want to generate a new modified GTF from the same reference.
- **-o** is the prefix for the output files (default = sqanti_sim)
- **-d** is the path to the directory to save the output files (default = .)
- You can choose how many transcripts from each structural category you want to delete from the reference using **--ISM**, **--NIC** and **--NNC** (default = 0).
- **-k** to choose the number of cores to run simultaneusly in the structural categorization
- **-v** to show up the version of the tool 

```
#!bash
$ python sqanti3_sim.py -h

usage: sqanti3_sim.py [-h] --gtf GTF [--cat CAT] [-o OUTPUT] [-d DIR]
                      [--ISM ISM] [--NIC NIC] [--NNC NNC] [--Fusion FUSION]
                      [--Antisense ANTISENSE] [--GG GG] [--GI GI]
                      [--Intergenic INTERGENIC] [--read_only] [-k CORES] [-v]

SQANTI-SIM: a simulator of controlled novelty and degradation of transcripts
sequence by long-reads

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Reference annotation in GTF format
  --cat CAT             File with transcripts structural categories generated
                        with SQANTI-SIM
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files. Default: Directory where
                        the script was run
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
  --read_only           If used the program will only categorize the GTF file
                        but skipping writing a new modified GTF
  -k CORES, --cores CORES
                        Number of cores to run in parallel
  -v, --version         Display program version number
```

## Classification criteria

![classification workflow](https://github.com/jorgemt98/SQANTI-SIM/blob/main/sqantisim_class_decision_tree.png)

## Example run

Example data can be found in `data/example_data.tar.gz`. Unpack it using `tar -xzf example_data.tar.gz` in `data` folder.

### Using a new reference GTF as input

`python sqanti3_sim.py --gtf data/example_ref.gtf -o test -d data --ISM 100 --NIC 200 --NNC 100`

### Using an existing SC classification file
`python sqanti3_sim.py --gtf data/example_ref.gtf --cat data/test_categories.txt -o test -d data --ISM 200 --NIC 100 --NNC 200`

## Output explanation

This tool generates mainly 2 output files:

- The modified GTF file used to "simulate" using it as reference annotation file, with the suffix **_modified.gtf**
- An intermediate file with the suffix **_categories.txt** that contains the classification in SQANTI3 structural categories for each transcript o the reference GTF (this is the ground-truth)

```
TransID	GeneID	SC	RefGene	RefTrans
ENST00000456328.2	ENSG00000223972.5	novel_not_in_catalog	ENSG00000223972.5	novel
ENST00000450305.2	ENSG00000223972.5	novel_not_in_catalog	ENSG00000223972.5	novel
ENST00000488147.1	ENSG00000227232.5	antisense	novelGene_ENSG00000243485.5_AS	novel
ENST00000619216.1	ENSG00000278267.1	genic_intron		novel
ENST00000473358.1	ENSG00000243485.5	novel_not_in_catalog	ENSG00000243485.5	novel
ENST00000469289.1	ENSG00000243485.5	incomplete-splice_match	ENSG00000243485.5	ENST00000473358.1
ENST00000607096.1	ENSG00000284332.1	incomplete-splice_match	ENSG00000243485.5	ENST00000469289.1
ENST00000417324.1	ENSG00000237613.2	novel_not_in_catalog	ENSG00000237613.2	novel
ENST00000461467.1	ENSG00000237613.2	incomplete-splice_match	ENSG00000237613.2	ENST00000417324.1
```

Moreover, in the terminal it shows a sumary table with all the transcripts classified for each structural category:

```
S Q A N T I - S I M 📊

Summary Table 🔎
_______________________________________________________________________________
| full-splice_match: 5564
| incomplete-splice_match: 40929
| novel_in_catalog: 73698
| novel_not_in_catalog: 78958
| fusion: 2139
| antisense: 5462
| genic_intron: 88
| genic: 1574
| intergenic: 28600
```


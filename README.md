<p align="center">
  <img src="https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/sqantisim_logo.png">
</p>

# SQANTI-SIM

**SQANTI-SIM** is a simulator of controlled novelty and degradation of transcripts sequenced by long reads. SQANTI-SIM is a wrapper tool for RNA-Seq long-reads simulators such as [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) and [NanoSim](https://github.com/bcgsc/NanoSim)(formerly Trans-NanoSim) to simulate novel and degraded transcripts based on the SQANTI3 structural categories ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI3)).

The aim of SQANTI-SIM is to simulate novel and degradated transcripts in a controlled way using as ground-truth real transcripts in the reference annotation that fit the description of the different SQANTI3 structural categories.

## Lastest updates

Current version (10/03/2020): SQANTI-SIM version alpha

Updates, patches and releases:

**alpha:**
- In development

## Short guide

### Requirements and Installation

In order to use the SQANTI-SIM pipeline it is necesarry to install its dependencies. You can install the SQANTI-SIM conda enviroment as follows:

```
conda env create -f SQANTI_SIM.conda_env.yml
export PYTHONPATH=$PYTHONPATH:<path_to_Cupcake>
export PYTHONPATH=$PYTHONPATH:<path_to_Cupcake>/sequence/
conda activate SQANTI-SIM
```

### Running SQANTI-SIM

The SQANTI-SIM pipeline consists mainly in 3 different stages: (i) simulate data, (ii) reconstruct the transcriptome with your pipeline and (iii) evaluate the performance of the transcript reconstruction pipeline.

#### 1. Simulation step

First, we must generate our transcript index file with all the transcripts form the reference annotation GTF classified in SQANTI3 structural categories

```
usage: sqanti_sim.py classif [-h] --gtf GTF [-o OUTPUT] [-d DIR]
                             [--min_ref_len MIN_REF_LEN] [-k CORES]

sqanti_sim.py classif parse options

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Reference annotation in GTF format
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length (default: 0 bp as
                        in largasp challenge 1 evaluation)
  -k CORES, --cores CORES
                        Number of cores to run in parallel
```

Then create your desired expression profile and generate the modified GTF that you must use in your transcript reconstruction pipeline

```
usage: sqanti_sim.py preparatory [-h] -i TRANS_INDEX --gtf GTF [-o OUTPUT]
                                 [-d DIR] [--read_count READ_COUNT]
                                 [-nt TRANS_NUMBER] [--nbn_known NBN_KNOWN]
                                 [--nbp_known NBP_KNOWN]
                                 [--nbn_novel NBN_NOVEL]
                                 [--nbp_novel NBP_NOVEL] [--rt RT]
                                 [--pb_reads PB_READS | --ont_reads ONT_READS]
                                 [--ISM ISM] [--NIC NIC] [--NNC NNC]
                                 [--Fusion FUSION] [--Antisense ANTISENSE]
                                 [--GG GG] [--GI GI] [--Intergenic INTERGENIC]
                                 [-k CORES] [-s SEED]
                                 {equal,custom,sample}

sqanti_sim.py preparatory parse options

positional arguments:
  {equal,custom,sample}
                        Different modes to generate the expression matrix:
                        equal (simulate with equal coverage for all reads),
                        custom (simulate with diferent negative binomial
                        distributions for novel and known transcripts) or
                        sample (simulate using a real sample)

optional arguments:
  -h, --help            show this help message and exit
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTI-SIM
  --gtf GTF             Reference annotation in GTF format
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  --read_count READ_COUNT
                        Number of reads to simulate (required for "equal"
                        mode)
  -nt TRANS_NUMBER, --trans_number TRANS_NUMBER
                        Number of different transcripts to simulate (required
                        for "equal" or "custom" mode)
  --nbn_known NBN_KNOWN
                        Average read count per known transcript to simulate
                        (i.e., the parameter "n" of the Negative Binomial
                        distribution) (required for "custom" mode)
  --nbp_known NBP_KNOWN
                        The parameter "p" of the Negative Binomial
                        distribution for known transcripts (required for
                        "custom" mode)
  --nbn_novel NBN_NOVEL
                        Average read count per novel transcript to simulate
                        (i.e., the parameter "n" of the Negative Binomial
                        distribution) (required for "custom" mode)
  --nbp_novel NBP_NOVEL
                        The parameter "p" of the Negative Binomial
                        distribution for novel transcripts (required for
                        "custom" mode)
  --rt RT               Reference transcripts in FASTA format (required for
                        "sample" mode)
  --pb_reads PB_READS   Input PacBio reads for quantification (required for
                        "sample" mode)
  --ont_reads ONT_READS
                        Input ONT reads for quantification (required for
                        "sample" mode)
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
  -s SEED, --seed SEED  Randomizer seed [123]
```

Finally, simulate the PacBio and/or ONT data. You can also simulate Illumina data to use it as input in your transcript discovery and reconstruction pipeline

```
usage: sqanti_sim.py sim [-h] --gtf GTF --genome GENOME [--rt RT] -i
                         TRANS_INDEX [--read_type READ_TYPE] [-d DIR]
                         [-k CORES] [--pb] [--ont] [--illumina]
                         [--read_count READ_COUNT] [-s SEED]

sqanti_sim.py sim parse options

optional arguments:
  -h, --help            show this help message and exit
  --gtf GTF             Reference annotation in GTF format
  --genome GENOME       Reference genome FASTA
  --rt RT               Reference transcripts in FASTA format
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTI-SIM
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
  --read_count READ_COUNT
                        Number of reads to simulate (if not given it will use
                        the counts of the given expression file)
  -s SEED, --seed SEED  Randomizer seed [123]
```

#### 2. Run transcript reconstruction pipeline

Run your transcript discovery and reconstruction pipeline with the modified GTF generate in the step before as your reference annotation GTF.

#### 3. Generate SQANTI-SIM report

Finally, run SQANTI3 with the modified GTF and generate the SQANTI-SIM report to evaluate the performance of the pipeline.

```
usage: sqanti_sim.py eval [-h] --isoforms ISOFORMS --gtf GTF --genome GENOME
                          -i TRANS_INDEX [-o OUTPUT] [-d DIR]
                          [--min_ref_len MIN_REF_LEN] [-k CORES]

sqanti_sim.py eval parse options

optional arguments:
  -h, --help            show this help message and exit
  --isoforms ISOFORMS   GTF with trancriptome reconstructed with your pipeline
  --gtf GTF             Reference annotation in GTF format
  --genome GENOME       Reference genome FASTA
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTI-SIM
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length (use the same as
                        in the classif step)
  -k CORES, --cores CORES
                        Number of cores to run in parallel
```

## Wiki

For more information visit our wiki

- [What is SQANTI-SIM?]()
- [SQANTI-SIM dependencies and installation]() 
- [Running SQANTI-SIM simulation]() 
- [Running SQANTI-SIM evaluation]() 
- [SQANTI-SIM output explanation]()
- [SQANTI-SIM example]()  

## How to cite SQANTI-SIM

There is not SQANTI-SIM paper :disappointed:

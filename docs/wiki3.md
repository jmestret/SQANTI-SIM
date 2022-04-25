## Parameters summary

Modes:

* [classif](#classif)
* [preparatory](#prep)
* [sim](#sim)
* [eval](#eval)


## <a name="classif"></a>classif

| Parameter | Required | Type | Description | Default
| --- | :---: | :---: | --- | :---: |
| --gtf | T | str | Reference annotationg in GTF format | |
| -o/--output | F | str | Prefix for output index file | sqanti_sim |
| -d/--dir | F | str | Directory for output files | . |
| -k/--cores | F | int | Number of cores to run in parallel | 1 |

## <a name="prep"></a>preparatory

In this mode you have three different sub modes which share some common arguments but have their own unique arguments.

**Common arguments**

| Parameter | Required | Type | Description | Default
| --- | :---: | :---: | --- | :---: |
| -i/--trans_index | T | str | File with transcript information generated with SQANTI-SIM classif | |
| --gtf | T | Reference annotationg in GTF format | |
| -o/--output | F | str | Prefix for output files | Same as index file |
| -d/--dir | F | str | Directory for output files | . |
| -nt/--trans_number | F | int | Number of different transcripts to simulate | 10000 |
| --ISM | F | int | Number of incomplete-splice-matches to simulate | 0 |
| --NIC | F | int | Number of novel-in-catalog to simulate | 0 |
| --NNC | F | int | Number of novel-not-in-catalog to simulate | 0 |
| --Fusion | F | int | Number of Fusion to simulate | 0 |
| --Antisense | F | int | Number of Antisense to simulate | 0 |
| --GG | F | int | Number of Genic-genomic to simulate | 0 |
| --GI | F | int | Number of Genic-intron to simulate | 0 |
| --Intergenic | F | int | Number of Intergenic to simulate | 0 |
| -k/--cores | F | int | Number of cores to run in parallel | 1 |
| -s/--seed | F | int | Randomizer seed | 123 |

**equal arguments**

| Parameter | Required | Type | Description | Default
| --- | :---: | :---: | --- | :---: |
| --read_count | F | int | Number of reads to simulate | 50000 |

**custom arguments**

| Parameter | Required | Type | Description | Default
| --- | :---: | :---: | --- | :---: |
| --nbn_known | F | float | Average read count per known transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution) | 15 |
| --nbp_knwon | F | float | The parameter 'p' of the Negative Binomial distribution for known transcripts | 0.5 | 
| --nbn_novel | F | float | Average read count per novel transcript to simulate (i.e., the parameter 'n' of the Negative Binomial distribution) | 5 |
| --nbp_novel | F | float | The parameter 'p' of the Negative Binomial distribution for novel transcripts | 0.5 |

**sample arguments**

| Parameter | Required | Type | Description | Default
| --- | :---: | :---: | --- | :---: |
| --rt | T | str | Reference transcripts in FASTA format | |
| --pb_reads/--ont_reads | T | str | Input PacBio or ONT reads for quantification | |
| --diff_exp | F | | If used the program will simulate different expression values for novel and known transcripts | F |
| --low_prob | F | float | Low value of prob vector (used if --diff_exp) | 0.25 |
| --high_prob | F | float | High value of prob vector (used if --diff_exp) | 0.75 |

## <a name="sim"></a>sim

## <a name="eval"></a>eval

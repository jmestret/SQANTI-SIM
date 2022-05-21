### Table of Contents

- [Introduction](#intro)
- [Usage](#use)
- [Arguments](#args)
- [Output explanation](#out)

## <a name="intro"></a>Introduction

`sqantisim.py eval` is the final step of the SQANTISIM pipeline. In this step, SQANTI3 is run, and SQANTISIM generates a report to evaluate the performance of the pipeline employed to identify the transcripts. This mode takes your long read-defined transcriptome, the modified reference annotation and the reference genome used in the steps before to evaluate the performance of your pipeline.

In this step, you can provide **orthogonal** data such as CAGE Peak and short-read support (that maybe you used in your reconstruction pipeline) to generate some metrics in the report regarding this information.

## <a name="use"></a>Usage

SQANTISIM *eval* mode usage:

```
sqantisim.py eval [-h] --isoforms ISOFORMS --gtf GTF --genome GENOME
                          -i TRANS_INDEX [-o OUTPUT] [-d DIR]
                          [--short_reads SHORT_READS] [--cage_peak CAGE_PEAK]
                          [--min_support MIN_SUPPORT] [-k CORES]
```

With the `--help` option you can display a complete description of the arguments:

```
sqantisim.py eval parse options

optional arguments:
  -h, --help            show this help message and exit
  --isoforms ISOFORMS   Transcriptome reconstructed with your pipeline
  --gtf GTF             Reference annotation in GTF format
  --genome GENOME       Reference genome FASTA
  -i TRANS_INDEX, --trans_index TRANS_INDEX
                        File with transcript information generated with
                        SQANTISIM
  -o OUTPUT, --output OUTPUT
                        Prefix for output files
  -d DIR, --dir DIR     Directory for output files (default: .)
  --short_reads SHORT_READS
                        File Of File Names (fofn, space separated) with paths
                        to FASTA or FASTQ from Short-Read RNA-Seq. If
                        expression or coverage files are not provided,
                        Kallisto (just for pair-end data) and STAR,
                        respectively, will be run to calculate them.
  --cage_peak CAGE_PEAK
                        FANTOM5 Cage Peak (BED format, optional)
  --min_support MIN_SUPPORT
                        Minimum number of supporting reads for an isoform
  -k CORES, --cores CORES
                        Number of cores to run in parallel
```

Running the tool with the minimum input to simulate PacBio reads will look as follows:

```
(SQANTISIM.env)$ python sqantisim.py eval \
			--isoforms long_read_transcriptome.gtf \
			--trans_index prefix_index.tsv \
			--gtf modified_reference_annotation.gtf \
			--genome reference_genome.fasta \
```

## <a name="args"></a>Arguments detailed explanation

### Required input

These are the minimum parameters you will need to run `sqantisim.py eval`:

- **Transcript index** file (`-i`): This file is the *prefix_index.tsv* file generated in the previous *sim* step.
- **Long-read transcriptome** (`--isoforms`): The isoforms identified with your transcript reconstruction pipeline. The transcripts models can be in GTF, FASTA and GTF format. **We recommend to use GTF format**. If you are going to input it as FASTA or FASTQ, you must add the `--fasta` argument.
- **Reference annotation** in GTF format (`--gtf`): This file is the modified GTF reference annotation generated in the *design* step that you should have used in your transcript reconstruction pipeline.
- **Reference genome** in FASTA format (`--genome`): This is the reference genome in FASTA format.

### Optional input

- **Short reads support** (`--short_reads`): File Of File Names (fofn, space-separated) with paths to FASTA or FASTQ from Short-Read RNA-Seq.
- **CAGE Peak data** (`--cage_peak`): FANTOM5 Cage Peak (BED format, optional)
- **Output prefix**(`-o`): The output prefix for the index file. SQANTISIM will use "sqantisim" as the default prefix.
- **Output directory**(`-d`): Output directory for output files. SQANTISIM will use the directory where the script was run as the default output directory.
- **Supporting reads** (`--min_support`): Minimum number of supporting reads for an isoform to be identified as a new isoform by your transcript reconstruction pipeline. This parameter doesn't affect the results; it just generates some metrics in the report according to your pipeline limitations/thresholds.
- **Parallelization**(`-k`): Number of cores to run in parallel. Most SQANTISIM modes have code chunks that can be run in parallel. However the default option is to run SQANTISIM in one single thread.

## <a name="out"></a>Output explanation

#### prefix_SQANTISIM_report.html

This is the main output file of the SQANTISIM pipeline. This HTML report picks up all the performance metrics of your transcript identification pipeline.

- **Total transcripts**: All simulated transcripts.
- **True Positives (TP)**: Reconstructed transcript models that match all the splice junctions with the true reference transcript and the difference in the TSS and TTS is less than 50 bp.
- **Partial True Positives (PTP)**: Reconstructed transcript models that match all the splice junctions with the true reference transcript, but differs badly from the annotated reference TSS and TTS.
- **False Positives (FP)**: Transcripts that were detected but weren't simulated.
- **False Negative (FN)**: Transcripts that were simulated but not detected.
- **Precision**: TP / (TP + FP)
- **Sensitivity**: TP / (TP + FN)
- **F-score**: 2 * ((Precision * Sensitivity) / (Precision + Sensitivity))
- **Positive Detection Rate**: (TP + PTP) / (TP + FN)
- **False Discovery Rate**: (FP + PTP) / (TP + FP)
- **False Detection Rate**: FP / (TP + FP)

#### SQANTI3 output files

This is a directory with all the SQANTI3 output. A detailed explanation of its output can be found in the [SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC).


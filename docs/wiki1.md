## What is SQANTISIM? <img src="https://github.com/jorgemt98/SQANTISIM/blob/main/docs/red_small_logo.png" alt="" width="100" />

*SQANTISIM* is a simulator of controlled novelty and degradation of transcripts sequenced by long reads. It is a wrapper tool for RNA-Seq long-reads simulators such as [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) and [NanoSim](https://github.com/bcgsc/NanoSim) (formerly Trans-NanoSim) to simulate transcripts based on the SQANTI3 structural categories ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI3)).

SQANTISIM aims to simulate novel and degraded transcripts in a controlled way using ground-truth present transcripts in the reference annotation that fit the description of the different SQANTI3 structural categories. SQANTISIM simulates reads from **true annotated transcripts** and modifies the GTF reference annotation pretending that the original transcripts are not annotated (novel). This allows to have a robust ground-truth and also to use **orthogonal data** (Short-read RNA-seq data, CAGE Peak data, polyA motif...) that can be used to support the discovery and reconstruction of this "novel" transcripts that were simulated.

![workflow](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/sqantisim_workflow.png)

SQANTISIM's strategy is based on SQANTI3, being able to simulate in a controlled way the different structural categories of SQANTI3. A detailed explanation of the SQANTI3 structural categories can be found in the [SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories). Briefly, degraded transcripts with fewer 5' start and/or 3' end exons are classified as an *Incomplete Splice Match* (ISM). On the other hand, novel isoforms are mainly classified as *Novel In Catalog* (NIC) or *Novel Not in Catalog* (NNC). NIC isoforms use a new combination of known donor/acceptor sites, while NNC isoforms have at least one donor or acceptor site that is not annotated. With SQANTISIM we are also able to simulate all the rest structural categories defined by SQANTI3.

<p align = "center">
<img src="https://github.com/FJPardoPalacios/public_figures/blob/master/figuras_class_SQ3.png" alt="Trulli" style="width:100%">
</p>

<p align = "center">
Figure taken from the SQANTI3 wiki
</p>

### How it works?

The SQANTISIM pipeline consists mainly of 3 different steps: (i) simulate data, (ii) reconstruct the transcriptome with your pipeline and (iii) evaluate the performance of the transcript reconstruction pipeline. `sqantisim.py` is a wrapper script with modules assumed to be run in the following order: classif, design, sim and eval. The first three modules will simulate the data, and the last module will generate the evaluation report.

**1. Transcript classification in SQANTI3 structural category**

This is the first step of the workflow. In this step, the SQANTI3 structural category decision tree is used to classify all the annotated transcripts from a given GTF reference annotation in its potential SQANTI3 structural categories, assigning the query transcript an SQANTI3 structural category assuming it cannot match as a FSM with itself. This will generate an index tab-separated file with all the structural annotations of each transcript from the reference GTF and its potential SQANTI3 structural category.

Information about running the classification step can be found at [*SQANTISIM: classif*](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/wiki5.md).

**2. Simulation experiment design**

This step takes the reference annotation in GTF format and the index file generated before and modifies the GTF according to the structural categories you want to simulate. It also generates the expression distribution and the TPM values for each transcript that will be simulated.Three different modes can generate transcript TPM values: using a uniform distribution for the count values, using two custom negative binomial distributions (one for the novel transcripts and the other for the known) and simulating the TPM values from the quantification of the expression levels of a transcriptome mapping with Minimap2 against the reference transcriptome.

Information about running the design step can be found at [*SQANTISIM: design*](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/wiki6.md).

**3. Data simulation**

This step takes the index file with the requested counts and the original GTF reference annotation and simulates ONT or PacBio reads using NanoSim and IsoSeqSim. You can also simulate Illumina reads using Polyester to support the long-read data. All three simulators use the *requested_tpm* column from the index file to generate the expression values for each read. You must use these simulated reads and the modified GTF reference annotation generated in the *design* step as input for the transcript discovery and reconstruction pipeline you are using.

Information about running the simulation step can be found at [*SQANTISIM: sim*](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/wiki7.md).

**4. SQANTISIM evaluation report**

Finally, this step takes the modified GTF reference annotation and a GTF with the transcript transcripts model generated by your pipeline. This will run the SQANTI3 pipeline and generate a report with performance metrics on how well the pipeline reconstructed novel and known transcripts and metrics of sensitivity and precision for each SQANTI3 structural category.

Information about running the evaluation step can be found at [*SQANTISIM: eval*](https://github.com/jorgemt98/SQANTISIM/blob/main/docs/wiki8.md).


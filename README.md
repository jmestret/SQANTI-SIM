<p align="center">
  <img src="https://github.com/jmestret/SQANTISIM/blob/main/docs/sqantisim_logo.png" alt="" width="300">
</p>

<p align="center">
  <a href="https://github.com/jmestret/SQANTISIM/wiki/Requirements-and-installation">Installation</a>
  ·
  <a href="https://github.com/jmestret/SQANTISIM/wiki">Documentation</a>
  ·
  <a href="https://github.com/jmestret/SQANTISIM/issues">Report bug</a>
  ·
  <a href="https://github.com/ConesaLab/SQANTI3">SQANTI3</a>
</p>

***

*SQANTISIM* is a simulator of controlled novelty and degradation of transcripts sequenced by long reads. It is a wrapper tool for RNA-Seq long-reads simulators such as [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) and [NanoSim](https://github.com/bcgsc/NanoSim) (formerly Trans-NanoSim) to simulate transcripts based on the SQANTI3 structural categories ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI3)).

SQANTISIM aims to simulate novel and degraded transcripts in a controlled way using ground-truth present transcripts in the reference annotation that fit the description of the different SQANTI3 structural categories.

![small_workflow](https://github.com/jmestret/SQANTISIM/blob/main/docs/small_workflow.png)

## Documentation

Please refer to [Wiki](https://github.com/jmestret/SQANTISIM/wiki) for how to install and use SQANTISIM 

## Lastest updates

Current version (14/07/2022): 0.1.1-beta

Updates, patches and releases:

**0.1.1-beta**:
- Updates in the SQANTISIM evaluation stage:
	- New evaluation metrics and charts have been added  (gene level metrics and expression level, number of exons and transcript length of true positives and false negatives).
	- An option for using RMSE, MAPE, and Q-Q plots to assess transcript quantification (`-e/--expression`).
	- The evaluation stage produces a new output file, a tab-separated file containing the performance metrics.
	- New table export buttons in the report.
	- The SQANTI3 report is no longer created in order to prevent undesired mistakes 
- Improved the isoform complexity approximation.
- Correct multithreading in NanoSim modified simulation.

0.1.0-beta:
- Please, try to use the tool and notify any bug or suggestion

## How to cite SQANTISIM

SQANTISIM will be included in the upcoming SQANTI3 paper

## Installation guide

* [Prerequisites](#prerequisites)
* [Install SQANTISIM](#installation)

## <a name="prerequisites"></a>Prerequisites

Following dependencies are required for running *SQANTISIM*:

#### Python packages

![python](https://img.shields.io/badge/python-3.7-blue)
* bx-python (Tested with version 0.8.13)
* BioPython (Tested with version 1.79)
* numpy (Tested with version 1.21.5)
* pandas (Tested with version 1.3.5)
* scipy (Tested with version 1.7.3)
* tqdm (Tested with version 4.63.0)

#### R packages

![R](https://img.shields.io/badge/R-4.1.2-blue)
* [polyester](https://github.com/alyssafrazee/polyester) (Tested with version 1.30.0)
* tidyverse (Tested with version X.X.X)
* DT (Tested with version 0.21)
* fmsb (Tested with version 0.7.3)
* griExtra (Tested with version 2.3)
* knitr (Tested with version 1.37)
* rmarkdown (Tested with version 2.13)

#### External programs and scripts

* [Minimap2](https://github.com/lh3/minimap2) (Tested with version 2.24)
* [STAR](https://github.com/alexdobin/STAR) (Tested with version 2.7.10a)
* [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) (Tested with version 28.0.0)
* [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) (Tested with version 0.1)
* [NanoSim requirements](https://github.com/bcgsc/NanoSim) (Tested with version 3.1.0)
* [SQANTI3 requirements](https://github.com/ConesaLab/SQANTI3) (Tested with version 5.0)

SQANTISIM has been tested on python 3.7 and R 4.1.2. You can check your python and R versions by:

```
python --version
# Python 3.7.12
R --version
# R version 4.1.2
```

## <a name="installation"></a>Install SQANTISIM

The recommended way to install all SQANTISIM dependencies is via Anaconda. [Here](https://docs.anaconda.com/anaconda/install/) is an explanation of how to install Anaconda. Installation is only necessary once. Please follow the steps below to avoid errors during installation.

#### 1. Install and update Anaconda

Make sure you have Anaconda installed and updated. In case Anaconda is not installed, please refer to its [documentation](https://docs.anaconda.com/anaconda/install/).

```
export PATH=$HOME/anacondaPy37/bin:$PATH
conda -V
conda update conda
```

#### 2. Download and set up SQANTISIM

Clone the [GitHub repository](https://github.com/jorgemt98/SQANTISIM) and create the conda environment with the dependencies.

```
git clone https://github.com/jorgemt98/SQANTISIM
cd SQANTISIM
conda env create -f SQANTISIM.conda_env.yml
source activate SQANTISIM.env
```

#### 3. Install cDNA_Cupcake

SQANTISIM makes use of cDNA_Cupcake, which cannot be automatically installed with Anaconda. If you don't have [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) installed, first activate the SQANTISIM environment and then follow the following commands:

```
(SQANTISIM.env)$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
(SQANTISIM.env)$ cd cDNA_Cupcake
(SQANTISIM.env)$ python setup.py build
(SQANTISIM.env)$ python setup.py install
```

After this, SQANTISIM is ready to use. To use SQANTISIM, just run the scripts directly; it does not require installation for itself. You can check whether the setting up is complete by:

```
(SQANTISIM.env)$ python sqantisim.py

# [SQANTISIM] usage: python sqantisim.py <mode> --help
# [SQANTISIM] modes: classif, preparatory, sim, eval
```
 

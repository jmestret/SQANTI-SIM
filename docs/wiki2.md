## Installation guide

* [Prerequisites](#prerequisites)
* [Install SQANTI-SIM](#installation)

## <a name="prerequisites"></a>Prerequisites

Following dependencies are required for running *SQANTI-SIM*:

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
* [SQANTI3 requirements](https://github.com/ConesaLab/SQANTI3) (Tested with version 4.2)

SQANTI-SIM has been tested on python 3.7 and R 4.1.2. You can check your python and R version by:

```
python --version
# Python 3.7.12
R --version
# R version 4.1.2
```

## <a name="installation"></a>Install SQANTI-SIM

The recommended way to install all SQANTI-SIM dependencies is via Anaconda. [Here](https://docs.anaconda.com/anaconda/install/) is an explanation on how to install Anaconda. Installation is only necessary onece. Please follow the steps below to avoid errors during installation.

**(1)** Make sure you have Anaconda installed and updated. In case Anaconda is not installed, please refer to its [documentation](https://docs.anaconda.com/anaconda/install/).

```
export PATH=$HOME/anacondaPy37/bin:$PATH
conda -V
conda update conda
```

**(2)** Clone the [GitHub repository](https://github.com/jorgemt98/SQANTI-SIM) and create the conda enviroment with the dependencies.

```
git clone https://github.com/jorgemt98/SQANTI-SIM
cd SQANTI-SIM
conda env create -f SQANTI_SIM.conda_env.yml
source activate SQANTI-SIM.env
```

**(3)** SQANTI-SIM makes use of cDNA_Cupcake which cannot be automatically installed with Anaconda. If you don't have [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake) installed, first activate the SQANTI-SIM enviroment and then follow the next commands:

```
(SQANTI-SIM.env)$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
(SQANTI-SIM.env)$ cd cDNA_Cupcake
(SQANTI-SIM.env)$ python setup.py build
(SQANTI-SIM.env)$ python setup.py install
```

After this, SQANTI-SIM is ready to use. To use SQANTI-SIM just run the scripts directly, it does not require installation for itself. You can check wether the setting up is complete by:

```
(SQANTI-SIM.env)$ python sqanti_sim.py

# [SQANTI-SIM] usage: python sqanti_sim.py <mode> --help
# [SQANTI-SIM] modes: classif, preparatory, sim, eval
```


 

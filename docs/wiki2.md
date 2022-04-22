## Installation guide

* [Prerequisites](#prerequisites)
* [Install SQANTI-SIM](#installation)

## <a name="prerequisites"></a>Prerequisites

Following dependencies are required for running *SQANTI-SIM*:

#### Python packages

![python](https://img.shields.io/badge/python-3.7-blue)
* bx-python (Tested with version X.X.X)
* BioPython (Tested with version X.X.X)
* numpy (Tested with version X.X.X)
* pandas (Tested with version X.X.X)
* scipy (Tested with version X.X.X)
* tqdm (Tested with version X.X.X)

#### R packages

![R](https://img.shields.io/badge/R-4.1.2-blue)
* [polyester](https://github.com/alyssafrazee/polyester) (Tested with version X.X.X)
* tidyverse (Tested with version X.X.X)
* DT (Tested with version X.X.X)
* fmsb (Tested with version X.X.X)
* griExtra (Tested with version X.X.X)
* knitr (Tested with version X.X.X)
* rmarkdown (Tested with version X.X.X)

#### External programs and scripts

* [Minimap2]() (Tested with version X.X.X)
* [STAR]() (Tested with version X.X.X)
* [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake)
* [IsoSeqSim]() 
* [NanoSim requirements](https://github.com/bcgsc/NanoSim)
* [SQANTI3 requirements](https://github.com/ConesaLab/SQANTI3)

SQANTI-SIM has been tested on python 3.7 and R 4.1.2. You can check your python and R version by:

```
python --version
# Python 3.7.12
R --version
# R version 4.1.2
```

## <a name="installation"></a>Install SQANTI-SIM

The recommended way to install all SQANTI-SIM dependencies is via Anaconda. [Here](https://docs.anaconda.com/anaconda/install/) is an explanation on how to install Anaconda. SQANTI-SIM makes use of cDNA_Cupcake which cannot be automatically installed with conda (more information in step X). Installation is only necessary onece. Please follow the steps below to avoid errors during installation.

**(1)** Make sure you have Anaconda installed and updated. In case Anaconda is not installed, please refer to it [documentation](https://docs.anaconda.com/anaconda/install/).

```
export PATH=$HOME/anacondaPy37/bin:$PATH
conda -V
conda update conda
```

**(2)** Download the package from [GitHub](https://github.com/jorgemt98/SQANTI-SIM) and install it locally.

```
git clone https://github.com/jorgemt98/SQANTI-SIM
```

**(3)** Change into the SQANTI-SIM top-level directory and create an Anaconda enviroment.


 

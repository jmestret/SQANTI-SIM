## SQANTISIM usage

#### 1. Activate the conda enviroment

Before running SQANTISIM, you will need to activate the conda environment. If you are running SQANTISIM from the terminal, use:

```
$ conda activate SQANTISIM.env 
```

On the other hand, if you are trying to use SQANTISIM from a bash script use:

```
$ source activate SQANTISIM.env 
```

#### 2. Add `cDNA_cupcake/sequence` to `$PYTHONPATH`

Once the environment is loaded, add *cDNA Cupcake* to the path:

```
(SQANTISIM.env)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(SQANTISIM.env)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/
```

#### 3. Run SQANTISIM

To run SQANTISIM, just call the python script `sqantisim.py` and choose one of the four running modes (*classif*, *design*, *sim* or *eval*) and specify the required and desired inputs. 

```
(SQANTISIM.env)$ usage: python sqantisim.py <mode> --help
(SQANTISIM.env)$ modes: classif, design, sim, eval
```

Refer to [*Usage summary*](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki3.md) for a brief explanation of the arguments and parameters for running SQANTISIM or use `python sqantisim.py <mode> --help`. For a detailed description of the different running modes refer to:
- [SQANTISIM: classif](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki5.md)
- [SQANTISIM: design](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki6.md)
- [SQANTISIM: sim](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki7.md)
- [SQANTISIM: eval](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki8.md)


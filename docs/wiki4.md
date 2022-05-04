## SQANTI-SIM usage

Refer to [*Usage summary*](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki3.md) for a brief explanation of the arguments and paramenters for running SQANTI-SIM.

**1. Activate the conda enviroment**

Before running SQANTI-SIM, you will need to activate the conda enviroment. If you are running SQANTI-SIM from the terminal use:

```
$ conda activate SQANTI-SIM.env 
```

On the other hand, if you trying to use SQANTI-SIM from a bash script use:

```
$ source activate SQANTI-SIM.env 
```

**2. Add `cDNA_cupcake/sequence` to `$PYTHONPATH`**

Once the enviroment is loaded, add *cDNA Cupcake* to the path:

```
(SQANTI-SIM.env)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(SQANTI-SIM.env)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/
```

**3. Run SQANTI-SIM**

SQANTI-SIM is a python script so it doesn't not require an installation itself. To run SQANTI-SIM just call the python script `sqanti_sim.py` and choose one of the 4 running modes (*classif*, *preparatory*, *sim* or *eval*) and specify the requiered and desired inputs. Refer to [*Usage summary*](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki3.md) for a brief explanation of the arguments and paramenters for running SQANTI-SIM or use `python sqanti_sim.py <mode> --help`.

For a detailed explanation of the different running modes refer to:
- [*SQANTI-SIM: classif*](https://github.com/jorgemt98/SQANTI-SIM/blob/main/docs/wiki5.md)
- [*SQANTI-SIM: preparatory*]()
- [*SQANTI-SIM: sim*]()
- [*SQANTI-SIM: eval*]()


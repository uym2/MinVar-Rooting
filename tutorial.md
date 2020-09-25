This is a mini tutorial where we give examples on running each rooting method (Outgroup, Midpoint, MinVar, Root-to-Tip) found in the ```use_cases``` folder. 

If you installed FastRoot using PyPI (i.e. pip), download [use_cases.zip](https://github.com/uym2/MinVar-Rooting/edit/master/use_cases.zip) to your machine and unzip it before trying the examples.

## Outgroup Rooting

An example is given in the folder `use_cases/OG`. Starting from the base directory,

```bash
   cd use_cases/OG
```

Inside this folder you will find a list of input trees (`input.trees`), a list of outgroups (`outgroups.txt`), and two bash scripts to root these trees using Outgroup Rooting. One script defines the outgroups as a file (`run_OG_file.sh`) while the other defines them as a list (`run_OG_list.sh`).
To root the input trees using this method, either run one of the scripts
```bash
   ./run_OG_file.sh
```
or call as follows:
```bash
   FastRoot.py -m OG -i input.trees -g outgroups.txt -o output.trees
```

## Midpoint Rooting

An example is given in the folder `use_cases/MP`. Starting from the base directory,

```bash
   cd use_cases/MP
```

Inside this folder you will find a list of input trees (`input.trees`) and a bash script to root these trees using Midpoint Rooting (`run_MP.sh`).
To root the input trees using this method, either run the script 
```bash
   ./run_MP.sh
```
or call as follows:
```bash
   FastRoot.py -m MP -i input.trees -o output.trees
```

## MinVar Rooting

An example is given in the folder `use_cases/MV`. Starting from the base directory,

```bash
   cd use_cases/MV
```

Inside this folder you will find a list of input trees (`input.trees`) and a bash script to root these trees using MinVar Rooting (`run_MV.sh`).
To root the input trees using this method, either run the script 
```bash
   ./run_MV.sh
```
or call as follows:
```bash
   FastRoot.py -m MV -i input.trees -o output.trees
```

## Root-to-Tip Rooting

An example is given in the folder `use_cases/RTT`. Starting from the base directory,

```bash
   cd use_cases/RTT
```

Inside this folder you will find a list of input trees (`input.trees`), the sampling times (`sampling_times.txt`), and a bash script to root these trees using Root-to-Tip Rooting (`run_RTT.sh`).
To root the input trees using this method, either run the script 
```bash
   ./run_RTT.sh
```
or call as follows:
```bash
   FastRoot.py -m RTT -i input.trees -t sampling_times.txt -o output.trees
```

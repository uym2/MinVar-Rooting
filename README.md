### If you find MinVar-Rooting helpful for your research, please cite the following paper
- Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. “Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction.” Edited by Gabriel Moreno-Hagelsieb. PLOS ONE 12, no. 8 (2017): e0182238. doi:10.1371/journal.pone.0182238.


### Besides, MinVar-Rootig was named MCCV in an independent and simultaneous work by Tria et al. Please also consider citing the following paper
- Domingues Kümmel Tria, Fernando & Landan, Giddy & Dagan, Tal. (2017). Phylogenetic rooting using minimal ancestor deviation. Nature Ecology & Evolution. 1. 0193. 10.1038/s41559-017-0193.


## Implementation
This is a dendropy-based implementation of the MinVar-Rooting, which seeks to root the tree at the point that minimizes the variance of the root to tip distances. The code was designed to be easily generalized for a class of 'optimization-based' rooting methods (some of which are under development), which includes a linear-time version of the traditional midpoint (MP) rooting as well.

Complexity: all rooting methods are linear (with the number of species) in time and memory.

## Dependencies
- Python3
- treeswift (version 1.1.14)
- cvxopt (version 1.2.5)
- numpy (version 1.19.0)

## Install using Pip
```bash
python3 -m pip install FastRoot
```

## Install from source code
1. Download the source code.  
	* Either clone the repository to your machine
	```bash
	   git clone https://github.com/uym2/MinVar-Rooting.git
	```
	* or simply download [this zip file](https://github.com/uym2/MinVar-Rooting/archive/master.zip) to your machine and unzip it in your preferred destination.

2. To install, go to the MinVar-Rooting folder. 
	* If you have ```pip```, use
	```bash
	   python3 -m pip install .
	```
	* Otherwise, if you have root access, type
	``` bash
	   sudo python3 setup.py install
	```
	* If you do not have root access, type
	``` bash
	   python3 setup.py install --user
	```
	
After installation, run:

```bash
FastRoot.py -h
```
to see the commandline help of FastRoot.

## Usage


```bash
FastRoot.py [-h] [-i INPUT] [-m METHOD] [-g OUTGROUPS] [-t SMPLTIMES] [-o OUTFILE] [-s SCHEMA] [-f INFOFILE]
```

optional arguments:
```		
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input File (default is STDIN)
  -m METHOD, --method METHOD
                        Method (MP for midpoint, MV for minVAR, OG for outgroup, RTT for root-to-tip) (default is MV)
  -g OUTGROUPS, --outgroups OUTGROUPS
                        Specify the outgroups. If specifying a list of outgroups, put them between quotes (i.e. ").
                        Otherwise, specifying a file which containts all the outgroups. Can only be used with -m OG
  -t SMPLTIMES, --smplTimes SMPLTIMES
                        The file containing the sampling times at leaves; to be used with -m RTT
  -o OUTFILE, --outfile OUTFILE
                        Output File (default is STDOUT)
  -s SCHEMA, --schema SCHEMA
                        Schema of your input treefile (default is newick)
  -f INFOFILE, --infofile INFOFILE
                        Save all the logging to this file. Default: print to stderr
  -v, --version         Show FastRoot version and exit
```

NOTE: `FastRoot.py` works for a list of trees

### Outgroups

* Outgroup Rooting `-m OG -g <OUTGROUPS>` maximizes the number of outgroup to ingroup triplets in the tree.
* Outgroups may be defined with a file `-g <OUTGROUP_FILE>` or with a list surrounded by quotation marks `-g "OUTGROUP1 OUTGROUP2 ..."`.
* The outgroups must be defined when using the outgroup rooting method.

### Sampling Times

* The sampling time file (example file: `use_cases/RTT/sampling_times.txt`) is a tab-delimited file, with one pair of species-time per line.
* It must have two columns: the species names and the corresponding sampling times.
* The sampling time for every leaf must be specified.
* The sampling times must be defined when using the root-to-tip rooting method.


 For example, lines

```
000009  9.36668
000010  9.36668
000011  11.3667
000012  11.3667
```
show that leaves `000009` and `000010` are sampled at time 9.36668 while nodes `000011` and `000012` are sampled at time 11.3667. 

**Note:** These times are assumed to be forward; i.e, smaller values mean closer to the root of the tree.

### Output
`FastRoot.py` with `-o` will output to the specified destination. Without `-o`, it prints the tree to standard output (stdout). The optimal score of each tree is printed to stderr by default; you can direct it to a file using `-f`.

### Example Usage

Below we give examples on running each rooting method (Outgroup, Midpoint, MinVar, Root-to-Tip) found in the ```use_cases``` folder. 

* If you installed FastRoot using PyPI, download [use_cases.zip](https://github.com/uym2/MinVar-Rooting/edit/master/use_cases.zip) to your machine and unzip it before trying the examples.
	
Note: the examples below assume you are using Linux or MacOS. For Windows users who installed FastRoot using the installation wizard, change ```FastRoot.py``` to ```FastRoot``` in all commands in the examples below.


#### Outgroup Rooting

An example is given in the folder `use_cases/OG`. Starting from the base directory,

```bash
   cd use_cases/OG
```

Inside this folder you will find a list of input trees (`input.trees`), a list of outgroups (`outgroups.txt`), and two bash scripts to root these trees using Outgroup Rooting. One scripts defines the outgroups as a file (`run_OG_file.sh`) while the other defines them as a list (`run_OG_list.sh`).
To root the input trees using this method, either run one of the scripts
```bash
   ./run_OG_file.sh
```
or call as follows:
```bash
   FastRoot.py -m OG -i input.trees -g outgroups.txt -o output.trees
```

#### Midpoint Rooting

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

#### MinVar Rooting

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

#### Root-to-Tip Rooting

An example is given in the folder `use_cases/RTT`. Starting from the base directory,

```bash
   cd use_cases/RTT
```

Inside this folder you will find a list of input trees (`input.trees`), the sampling times (`sampling_times.txt`) and a bash script to root these trees using Root-to-Tip Rooting (`run_RTT.sh`).
To root the input trees using this method, either run the script 
```bash
   ./run_RTT.sh
```
or call as follows:
```bash
   FastRoot.py -m RTT -i input.trees -t sampling_times.txt -o output.trees
```

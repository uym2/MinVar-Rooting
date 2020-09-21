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

### Sampling Times

* The sampling time file (example file: `use_cases/RTT/sampling_times.txt`) is a tab-delimited file, with one pair of species-time per line.
* It must have two columns: the species names and the corresponding sampling times.
* The sampling time for every leaf must be specified.
* This file is necessary for the root-to-tip rooting method.


 For example, lines

```
000009  9.36668
000010  9.36668
000011  11.3667
000012  11.3667
```
show that leaves `000009` and `000010` are sampled at time 9.36668 while nodes `000011` and `000012` are sampled at time 11.3667. 

**Note:** These times are assumed to be forward; i.e, smaller values mean closer to the root of the tree. The top of the branch above the root is assumed to be 0.

## Output
`FastRoot.py` with `-o` will output to the specified destination. Without `-o`, it prints the tree to standard output (stdout). The optimal score of each tree is printed to stderr by default; you can direct it to a file using `-f`.

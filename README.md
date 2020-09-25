FastRoot is a python implementation of a class of 'optimization-based' rooting methods for phylogenetic trees: Minimum Variance (MV), Midpoint (MP), Outgroups (OG), and Root-to-tip (RTT). All rooting methods are linear in time and memory. FastRoot was tested on Linux or MacOS.

## Installation
FastRoot needs the following dependencies
- Python3
- treeswift (version 1.1.14)
- cvxopt (version 1.2.5)
- numpy (version 1.19.0)

You can either install from PyPI (i.e. pip) or from source code.

### Install using Pip
You need to have Python3 and Pip (in most cases are already installed with Python3, if not, please see https://pip.pypa.io/en/stable/installing/). All the other dependencies of FastRoot will be automatically installed together with the package.
```bash
python3 -m pip install FastRoot
```

### Install from source code
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

### Rooting methods
There are 4 rooting methods: minVAR (MV), midpoint (MP), outgroups (OG), and root-to-tip (RTT).

#### MinVAR Rooting (MV)
* Usage: `-m MV`
* Root the tree at the point that minimizes the variance of the root to tip distances.

#### Midpoint Rooting (MP)
* Usage: `-m MP`
* Root the tree at the midpoint of the longest path between any pair of leaves (i.e. midpoint of the diameter).

#### Root-to-tip Rooting (RTT)
* Usage: `-m RTT -t <SAMPLING_TIMES>`
* Optimizes the least squares regression of the root-to-tip time and substitutions.
* The sampling times MUST be defined via ```-t```: a tab-delimited file, with one pair of species-time per line.
* Example file: `use_cases/RTT/sampling_times.txt`.
For example, lines

```
000009  9.36668
000010  9.36668
000011  11.3667
000012  11.3667
```
show that leaves `000009` and `000010` are sampled at time 9.36668 while nodes `000011` and `000012` are sampled at time 11.3667. 

**Note:** 
- The sampling time for every leaf must be specified.
- These times are assumed to be forward; i.e, smaller values mean closer to the root of the tree.

#### Outgroup Rooting (OG)
* Usage: `-m OG -g <OUTGROUPS>`
* Maximizes the number of outgroup to ingroup triplets in the tree (i.e. maximizes the number of triplets of the forms (o,(i,i)) and (i,(o,o)) where i is an ingroup species and o is an outgroup species).
* The outgroups MUST be defined via ```-g```: can either by a file `-g <OUTGROUP_FILE>` or a list surrounded by quotation marks `-g "OUTGROUP1 OUTGROUP2 ..."`.

### Output
`FastRoot.py` with `-o` will output to the specified destination. Without `-o`, it prints the tree to standard output (stdout). 
The optimal score of each tree (depends on the rooting method) is printed to stderr by default; you can direct it to a file using `-f`.

### Example usage
You can checkout the [tutorial](https://github.com/uym2/MinVar-Rooting/blob/master/tutorial.md) for some examples.

## Publications
If you find MinVar-Rooting helpful for your research, please cite the following paper
- Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. “Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction.” Edited by Gabriel Moreno-Hagelsieb. PLOS ONE 12, no. 8 (2017): e0182238. doi:10.1371/journal.pone.0182238.


#### Besides, MinVar-Rootig was named MCCV in an independent and simultaneous work by Tria et al. Please also consider citing the following paper
- Domingues Kümmel Tria, Fernando & Landan, Giddy & Dagan, Tal. (2017). Phylogenetic rooting using minimal ancestor deviation. Nature Ecology & Evolution. 1. 0193. 10.1038/s41559-017-0193.

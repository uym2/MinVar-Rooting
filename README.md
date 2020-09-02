### If you find MinVar-Rooting helpful for your research, please cite the following paper
- Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. “Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction.” Edited by Gabriel Moreno-Hagelsieb. PLOS ONE 12, no. 8 (2017): e0182238. doi:10.1371/journal.pone.0182238.


### Besides, MinVar-Rootig was named MCCV in an independent and simultaneous work by Tria et al. Please also consider citing the following paper
- Domingues Kümmel Tria, Fernando & Landan, Giddy & Dagan, Tal. (2017). Phylogenetic rooting using minimal ancestor deviation. Nature Ecology & Evolution. 1. 0193. 10.1038/s41559-017-0193.


## Implementation
This is a dendropy-based implementation of the MinVar-Rooting, which seeks to root the tree at the point that minimizes the variance of the root to tip distances. The code was designed to be easily generalized for a class of 'optimization-based' rooting methods (some of which are under development), which includes a linear-time version of the traditional midpoint (MP) rooting as well.

Complexity: all rooting methods are linear (with the number of species) in time and memory.

## Dependencies
- Python (either 2.x or 3.x)
- treeswift (version 1.1.14)
- cvxopt (version 1.2.5)
- numpy (version 1.19.0)

### Install using Pip
```bash
pip install FastRoot==1.0
```

## Install from source code
1. Download the source code.  

	```bash
	   git clone https://github.com/uym2/MinVar-Rooting.git
	```

2. To install, go to the MinVar-Rooting folder. 
	* If you have ```pip```, use
	```bash
	   python3 -m pip install .
	```
	* Otherwise, type
	``` bash
	   python3 setup.py install
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
                        Method (MP for midpoint, MV for minVAR, OG for outgroup, RTTas or RTTqp for root-to-tip)
                        (default is MV)
  -g OUTGROUPS, --outgroups OUTGROUPS
                        Listing of the outgroups; to be used with -m OG
  -t SMPLTIMES, --smplTimes SMPLTIMES
                        The file containing the sampling times at leaves; to be used with -m RTT
  -o OUTFILE, --outfile OUTFILE
                        Output File (default is STDOUT)
  -s SCHEMA, --schema SCHEMA
                        Schema of your input treefile (default is newick)
  -f INFOFILE, --infofile INFOFILE
                        Report the optimization score to file
```

NOTE: `FastRoot.py` works for a list of trees

## Output
`FastRoot.py` with `-o` will output to the specified destination. Without `-o`, it prints the tree to standard output.

## Pseudocode
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/MV_alg.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/Eq4.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/Eq6.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/Eq7.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/MP_alg.png)

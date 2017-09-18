## Mathematical properties of the MinVar-Rooting method and justifications of the linear algorithm can be found in the following paper
- Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. “Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction.” Edited by Gabriel Moreno-Hagelsieb. PLOS ONE 12, no. 8 (2017): e0182238. doi:10.1371/journal.pone.0182238.
          - Best paper award at GLBIO 2017.
- Besides, MinVar-Rootig is named MCCV in an independent and simultaneous work published in Nature Ecology and Evolution.
Domingues Kümmel Tria, Fernando & Landan, Giddy & Dagan, Tal. (2017). Phylogenetic rooting using minimal ancestor deviation. Nature Ecology & Evolution. 1. 0193. 10.1038/s41559-017-0193. 
  

## Implementation
This is a dendropy-based implementation for a new idea of rooting a phylogenetic tree, called the minVar (MV) root. Our rooting method roots the tree at the point that minimizes the variance of the root to tip distances. The code was designed to be easily generalized for a class of 'optimization-based' rooting methods (some of which are under development), which includes a linear-time version of the traditional midpoint (MP) rooting as well.

Complexity: all rooting methods are linear (with the number of species) in time and memory.

## Dependencies
- python (either 2.x or 3.x)
- dendropy (version 4.2.0 recommended)

## Usage

```python
python FastRoot.py [-h] -i INPUT -m METHOD [-o OUTFILE] [-s SCHEMA]
```
```
  -h, --help            show help message and exit
  
  -i INPUT, --input INPUT:
                        input file
                        
  -m METHOD, --method METHOD
                        method: MP for midpoint and MV for minVAR. Default is MV
                        
  -o OUTFILE, --outfile OUTFILE
                        specify output file
                        
  -s SCHEMA, --schema SCHEMA
                        schema of your input treefile. Default is Newick
```

NOTE: FastRoot.py works for a list of trees

The following scripts might be more handy but only work with single-tree input

```python 
python MP_reroot.py <tree_file_in_newick>
python MV_reroot.py <tree_file_in_newick>
```

## Output
The FastRoot.py with -o will output to the specified destination. Without -o, it prints the tree to stdout. 
The MP\_reroot.py and MV\_reroot.py will place the output in the same directory as the input tree with suffix MP\_rooted or MV\_rooted accordingly.

## Pseudocode
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/MV_alg.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/Eq4.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/Eq6.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/Eq7.png)
![alt tag](https://github.com/uym2/MinVar-Rooting/blob/master/imgs/MP_alg.png)

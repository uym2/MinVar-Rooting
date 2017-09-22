#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import MPR_Tree,MVR_Tree
from dendropy import Tree,TreeList

from sys import stdin,stdout
import argparse

METHOD2FUNC = {'MP':MPR_Tree,'MV':MVR_Tree}

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',required=False,type=argparse.FileType('r'),default=stdin,help="Input File (default is STDIN)")
parser.add_argument('-m','--method',required=False,type=str,default="MV",help="Method (MP for midpoint, MV for minVAR) (default is MV)")
parser.add_argument('-o','--outfile',required=False,type=argparse.FileType('w'),default=stdout,help="Output File (default is STDOUT)")
parser.add_argument('-s','--schema',required=False,type=str,default="newick",help="Schema of your input treefile (default is newick)")
parser.add_argument('-f','--infofile',required=False,type=argparse.FileType('w'),default=None,help="Write info of the new root to file (mostly for research and debugging purposes) (default is None)")
args = parser.parse_args()
assert args.method in METHOD2FUNC, "Invalid method! Valid options: MP for midpoint, MV for minVAR"

for line in args.input:
	tree = Tree.get(data=line,schema=args.schema.lower(),preserve_underscores=True)
	a_tree = METHOD2FUNC[args.method](ddpTree=tree)
	d2currRoot,br2currRoot = a_tree.Reroot()

	if args.infofile:
		args.infofile.write("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")

	a_tree.tree_as_newick(outfile=args.outfile,append=True)

if args.infofile:
	args.infofile.close()

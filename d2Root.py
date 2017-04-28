#! /usr/bin/env python

# usage: python LabelTree.py <tree_file>

import os
from Tree_extend import Tree_extend,MVR_Tree,MPR_Tree,MBR_Tree
from dendropy import Tree,TreeList

from sys import argv
from os.path import splitext
import argparse
import optparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-m','--method',required=False,help="rooting method; default is MV")
parser.add_argument('-s','--schema',required=False,help="schema of your input treefile. Default is newick")

args = vars(parser.parse_args())

tree_file = args["input"]
base_name,ext = splitext(tree_file)
schema=args["schema"] if args["schema"] else "newick"


trees = TreeList.get_from_path(tree_file,schema)

for tree in trees:
        if args["method"] == "MV":
        	a_tree = MVR_Tree(ddpTree=tree)
        elif args["method"] == "MP":
        	a_tree = MPR_Tree(ddpTree=tree)
        elif args["method"] == "MB":
        	a_tree = MBR_Tree(ddpTree=tree)
	a_tree.Reroot()
        #a_tree.tree_as_newick(outfile=outfile,append=True,label_by_name=True)
        D = a_tree.compute_distances()
        print(D)

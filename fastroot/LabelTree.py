#! /usr/bin/env python

# usage: python LabelTree.py <tree_file>

import os
from fastroot.Tree_extend import Tree_extend
#from dendropy import Tree,TreeList
from treeswift import *

from sys import argv
from os.path import splitext
import argparse
import optparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-o','--outfile',required=False,help="specify output file")
parser.add_argument('-s','--schema',required=False,help="schema of your input treefile. Default is newick")
parser.add_argument('-l','--label',required=False,help="labeling style: 'leaves', 'all', or 'internal'. Default: all")

args = vars(parser.parse_args())

tree_file = args["input"]
base_name,ext = splitext(tree_file)
schema=args["schema"] if args["schema"] else "newick"
style = args["label"] if args["label"] else "all"

if args["outfile"]:
	outfile = args["outfile"]
else:
	#outfile = base_name + "_labeled" + ext
        outfile = None
try:
	os.remove(outfile)
except:
	pass

trees = read_tree(tree_file,schema)

for tree in trees:
    a_tree = Tree_extend(ddpTree=tree)
    a_tree.Topdown_label(label_type=style)
    a_tree.tree_as_newick(outfile=outfile,append=True,label_by_name=True)

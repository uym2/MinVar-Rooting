#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import MPR_Tree, OGR_Tree
from MinVar import *
from RTT import *
#from dendropy import Tree,TreeList
from treeswift import *

from sys import stdin, stdout
import argparse

METHOD2FUNC = {'MP': MPR_Tree, 'MV': MV00_Tree, 'OG': OGR_Tree, 'RTT': RTT_Tree}

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=False, type=argparse.FileType('r'), default=stdin,
                    help="Input File (default is STDIN)")
#parser.add_argument('-m', '--method', required=False, type=str, default="MV",
                   # help="Method (MP for midpoint, MV for minVAR, OG for outgroup) (default is MV)")
parser.add_argument('-m', '--method', required=False, type=str, default="RTT",
                    help="Method (MP for midpoint, MV for minVAR, OG for outgroup) (default is MV)")
# temporarily changed default to RTT
parser.add_argument('-g', '--outgroups', required=False, type=str,
                    help="Listing of the outgroups; to be used with -m OG")
parser.add_argument('-t', '--smplTimes', required=False, type=argparse.FileType('r'),
                    help="The file containing the sampling Times at leaves; to be used with -m RTT")
parser.add_argument('-o', '--outfile', required=False, type=argparse.FileType('w'), default=stdout,
                    help="Output File (default is STDOUT)")
parser.add_argument('-s', '--schema', required=False, type=str, default="newick",
                    help="Schema of your input treefile (default is newick)")
parser.add_argument('-f', '--infofile', required=False, type=argparse.FileType('w'), default=None,
                    help="Report the optimization score to file")

args = parser.parse_args()

assert args.method in METHOD2FUNC, "Invalid method! Valid options: MP for midpoint, MV for minVAR"

OGs = args.outgroups.split() if args.outgroups else None
smplTimes = {}
for line in args.smplTimes:
    sp,t = line.strip().split()
    smplTimes[sp] = float(t)

for line in args.input:
    tree = read_tree(line, schema=args.schema.lower())
    if args.method == 'OG':
        a_tree = OGR_Tree(OGs, ddpTree=tree)
    elif args.method == 'RTT':
        a_tree = RTT_Tree(smplTimes, ddpTree=tree)
    else:
        a_tree = METHOD2FUNC[args.method](ddpTree=tree)

    #a_tree.Bottomup_update()
    #a_tree.prepare_root() ########
    a_tree.Reroot()
    #print(a_tree.report_score())

    if args.infofile:
        args.infofile.write(a_tree.report_score() + "\n")

    # args.outfile.write(a_tree.ddpTree.as_string("newick"))

    a_tree.tree_as_newick(outstream=args.outfile)

# if args.infofile:
# args.infofile.close()

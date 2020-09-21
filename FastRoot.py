#! /usr/bin/env python3

from fastroot.Tree_extend import MPR_Tree, OGR_Tree
from fastroot.MinVar import *
from fastroot.RTT import *
from treeswift import *
import fastroot
from sys import stdin, stdout, argv, exit, stderr
import argparse
from os import path
import re

def main():

    # parse arguments
    parser = argparse.ArgumentParser()
    
    
    parser.add_argument('-i', '--input', required=False, type=argparse.FileType('r'), default=stdin,
                        help="Input File (default is STDIN)")
    parser.add_argument('-m', '--method', required=False, type=str, default="MV",
                        help="Method (MP for midpoint, MV for minVAR, OG for outgroup, RTT for root-to-tip) (default is MV)")
    parser.add_argument('-g', '--outgroups', required=False, type=str,
                        help="Specify the outgroups. If specifying a list of outgroups, put them between quotes (i.e. \"). Otherwise, specifying a file which containts all the outgroups. Can only be used with -m OG")
    parser.add_argument('-t', '--smplTimes', required=False, type=argparse.FileType('r'),
                        help="The file containing the sampling times at leaves; to be used with -m RTT")
    parser.add_argument('-o', '--outfile', required=False, type=argparse.FileType('w'), default=stdout,
                        help="Output File (default is STDOUT)")
    parser.add_argument('-s', '--schema', required=False, type=str, default="newick",
                        help="Schema of your input treefile (default is newick)")
    parser.add_argument('-f', '--infofile', required=False, type=argparse.FileType('w'), default=None,
                        help="Save all the logging to this file. Default: print to stderr")
    parser.add_argument("-v", "--version", action='version',
                        version=fastroot.PROGRAM_NAME + " " + fastroot.PROGRAM_VERSION,
                        help="Show FastRoot version and exit")
    
    # print help message if no argument is given
    if len(argv) == 1:
        logger = fastroot.new_logger(__name__)
        logger.info("Running " +  fastroot.PROGRAM_NAME +  " version " + fastroot.PROGRAM_VERSION) 
        parser.print_help()
        exit(0) 

    args = parser.parse_args()
    stream = args.infofile if args.infofile else stderr
    logger = fastroot.new_logger(__name__,myStream=stream)
    logger.info("Running " +  fastroot.PROGRAM_NAME +  " version " + fastroot.PROGRAM_VERSION)
    new_argv = re.sub(' +', ' ', " ".join(argv).replace("\n", " "))
    logger.info(fastroot.PROGRAM_NAME + " was called as follows: " + new_argv)

    METHOD2FUNC = {'MP': MPR_Tree, 'MV': MV00_Tree, 'OG': OGR_Tree, 'RTT': RTT_Tree}
    METHOD2DESC = {'MP': "Midpoint", 'MV': "MinVar", 'OG': "Outgroup", 'RTT': "Root-to-Tip"}
    method = args.method.upper()

    # reading outgroups
    if args.outgroups:
        if method != 'OG':
            method = 'OG'
            logger.warning("The rooting method is set to outgroup rooting (OG) due to the presence of outgroups")
        if path.exists(args.outgroups):
            OGs = []
            for line in open(args.outgroups, 'r'):
                OGs.append(line.strip())
        else:
            OGs = args.outgroups.split()
    else:
        OGs = None

    # reading sampling times
    if args.smplTimes:
        smplTimes = {}
        for line in args.smplTimes:
            sp,t = line.strip().split()
            smplTimes[sp] = float(t)
        if method != 'RTT':
            method = 'RTT'
            logger.warning("The rooting method is set to root-to-tip rooting (RTT) due to the presence of sampling times")

    if method == 'RTT' and args.smplTimes is None:
        logger.error("Need sampling times for root-to-tip rooting")
        exit()
    elif method == 'OG' and args.outgroups is None:
        logger.error("Need outgroups for outgroup rooting")
        exit()

    assert method in METHOD2FUNC, "Invalid method! Valid options: MP for midpoint, MV for minVAR, OG for outgroups, RTT for root-to-tip"
    logger.info("Rooting Method: " + METHOD2DESC[method] + " Rooting")

    # read and root each tree
    for i,line in enumerate(args.input):
        tree = read_tree(line, schema=args.schema.lower())
        if method == 'OG':
            a_tree = OGR_Tree(OGs, ddpTree=tree,logger_id=i+1,logger_stream=stream)
        elif method == 'RTT':
            a_tree = RTT_Tree(smplTimes, ddpTree=tree,logger_id=i+1,logger_stream=stream)
        else:
            a_tree = METHOD2FUNC[method](ddpTree=tree,logger_id=i+1,logger_stream=stream)

        a_tree.Reroot()
        logger.info("Tree " + str(i+1) + " " + a_tree.report_score())
        a_tree.tree_as_newick(outstream=args.outfile)

if __name__ == "__main__":
    main()

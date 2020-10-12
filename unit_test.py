import os
import unittest
from fastroot.Tree_extend import MPR_Tree, OGR_Tree
from treeswift import *
from fastroot.MinVar import MV00_Tree, minVAR_Base_Tree
from fastroot.RTT import RTT_Tree
from random import gauss
#from sys import stdin, stdout, argv, exit, stderr
#import argparse
#import re


# check if two newick strings are equal (same clades):
def compute_clades(tree):
    clades = []
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.clade = [node.label]
        else:
            node.clade = []
            for child in node.children:
                node.clade += child.clade
        clades.append(node.clade)
    return clades

def check_clades_include(cl1, cl2):
    # check if cl1 is a subset of cl2
    for c1 in cl1:
        found = False
        for c2 in cl2:
            if sorted(c1) == sorted(c2):
                found = True
                break
        if not found:
            return False
    return True

def check_two_nwk_str(nwk1, nwk2):
    tree1 = read_tree_newick(nwk1)
    tree2 = read_tree_newick(nwk2)

    L1 = sorted([x.label for x in tree1.traverse_leaves()])
    L2 = sorted([x.label for x in tree2.traverse_leaves()])
    if L1 != L2:
        return False

    cld1 = compute_clades(tree1)
    cld2 = compute_clades(tree2)

    return check_clades_include(cld1, cld2)


def root_trees(inputFile,method='MV',timeFile=None,OGFile=None,t0=0.0):
    # define OGs in file, not list
    # reading outgroups
    if OGFile:
        method = 'OG'
        OGs = []
        f = open(OGFile, 'r')
        for line in f:
                OGs.append(line.strip())
        f.close()

    # reading sampling times
    if timeFile:
        smplTimes = {}
        f = open(timeFile, 'r')
        for line in f:
            sp, t = line.strip().split()
            smplTimes[sp] = float(t) + t0
        method = 'RTT'
        f.close()

    METHOD2FUNC = {'MP': MPR_Tree, 'MV': MV00_Tree, 'OG': OGR_Tree, 'RTT': RTT_Tree}

    # read and root each tree
    score = {}  # dictionary of tree number and its score
    tree_list = []
    branches = {}
    inputTrees = open(inputFile)
    for i, line in enumerate(inputTrees.readlines()):
        tree = read_tree(line, 'newick')
        if method == 'OG':
            a_tree = OGR_Tree(OGs, ddpTree=tree, logger_id=i + 1)
        elif method == 'RTT':
            a_tree = RTT_Tree(smplTimes, ddpTree=tree, logger_id=i + 1)
        else:
            a_tree = METHOD2FUNC[method](ddpTree=tree, logger_id=i + 1)

        a_tree.Reroot()

        if method == 'RTT':
            score[i+1] = a_tree.return_values()  # tree number: RTT score, mu, t0
        else:
            score[i+1]= a_tree.opt_score()

        tree_list.append(a_tree.ddpTree.newick())

        branches[i + 1] = sorted([child.edge_length for child in a_tree.ddpTree.root.child_nodes()])
    inputTrees.close()
    return score,tree_list, branches


def branch_lengths(file):
    branch_lengths = {}
    f = open(file,'r')
    for i, line in enumerate(f):
        l1,l2 = line.strip().split()
        branch_lengths[i+1] = sorted([float(l1), float(l2)])
    f.close()
    return branch_lengths

def score_from_file(file, RTT=False):
    score = {}
    f = open(file, 'r')
    if RTT:
        for i,line in enumerate(f):
            RTT_score, mu, t0 = line.strip().split()
            score[i+1] = [float(RTT_score), float(mu), float(t0)]
    else:
        for i, line in enumerate(f):
            score[i+1] = float(line.strip())
    f.close()
    return score


class RootingTestCase(unittest.TestCase):
    """Tests for `FastRoot.py`."""

    def test_OG(self):
        """Does OG root properly?"""
        print("Testing OG")
        EPSILON_SCORE = 0.0001

        OG = open("unit_test/OG/output.trees")
        score = score_from_file("unit_test/OG/score.txt")
        #branches = branch_lengths("unit_test/OG/branches.txt")
        score_test, OG_test, branches_test = root_trees("unit_test/OG/input.trees", method='OG', OGFile="unit_test/OG/outgroups.txt")
        for line, line_test in zip(OG.readlines(), OG_test):
            self.assertTrue(check_two_nwk_str(line, line_test), msg="OG Rooting failed: incorrect newick string.")
        for i in score:
            self.assertTrue(score[i] - score_test[i] < EPSILON_SCORE, msg="OG Rooting failed: incorrect triplet distance score.")
        #for i in branches:
        #    self.assertTrue(branches[i][0] - branches_test[i][0] < 0.001)
        #    self.assertTrue(branches[i][1] - branches_test[i][1] < 0.001)
        OG.close()

    def test_MP(self):
        """Does MP root properly?"""
        print("Testing MP")  #
        EPSILON_SCORE = 0.0001
        EPSILON_BRANCH = 0.0001

        MP = open("unit_test/MP/output.trees")
        score = score_from_file("unit_test/MP/score.txt")
        branches = branch_lengths("unit_test/MP/branches.txt")
        score_test, MP_test, branches_test = root_trees("unit_test/MP/input.trees", method='MP')
        for line, line_test in zip(MP.readlines(), MP_test):
            self.assertTrue(check_two_nwk_str(line, line_test), msg="MP Rooting failed: incorrect newick string.")
        for i in score:
            self.assertTrue(score[i]-score_test[i] < EPSILON_SCORE, msg="MP Rooting failed: incorrect tree height.")
        for i in branches:
            self.assertTrue(branches[i][0]-branches_test[i][0] < EPSILON_BRANCH, msg="MP Rooting failed: incorrect branch length.")
            self.assertTrue(branches[i][1] - branches_test[i][1] < EPSILON_BRANCH, msg="MP Rooting failed: incorrect branch length.")
        MP.close()


        MP = open("unit_test/MP/special/output.trees")
        score = score_from_file("unit_test/MP/special/score.txt")
        branches = branch_lengths("unit_test/MP/special/branches.txt")
        score_test, MP_test, branches_test = root_trees("unit_test/MP/special/input.trees", method='MP')
        for line, line_test in zip(MP.readlines(), MP_test):
            self.assertTrue(check_two_nwk_str(line, line_test), msg="MP Rooting failed: incorrect newick string.")
        for i in score:
            self.assertTrue(score[i] - score_test[i] < EPSILON_SCORE, msg="MP Rooting failed: incorrect tree height.")
        for i in branches:
            self.assertTrue(branches[i][0] - branches_test[i][0] < EPSILON_BRANCH,
                            msg="MP Rooting failed: incorrect branch length.")
            self.assertTrue(branches[i][1] - branches_test[i][1] < EPSILON_BRANCH,
                            msg="MP Rooting failed: incorrect branch length.")
        MP.close()

    def test_MV(self):
        """Does MV root properly?"""
        print("Testing MV")
        EPSILON_SCORE = 0.0001
        EPSILON_BRANCH = 0.0001

        MV = open("unit_test/MV/output.trees")
        score = score_from_file("unit_test/MV/score.txt")
        branches = branch_lengths("unit_test/MV/branches.txt")
        score_test, MV_test, branches_test = root_trees("unit_test/MV/input.trees", method='MV')
        for line, line_test in zip(MV.readlines(), MV_test):
            self.assertTrue(check_two_nwk_str(line, line_test), msg="MV Rooting failed: incorrect newick string.")
        for i in score:
            self.assertTrue(score[i] - score_test[i] < EPSILON_SCORE, msg="MV Rooting failed: incorrect MV score.")
        for i in branches:
            self.assertTrue(branches[i][0] - branches_test[i][0] < EPSILON_BRANCH, msg="MV Rooting failed: incorrect branch length.")
            self.assertTrue(branches[i][1] - branches_test[i][1] < EPSILON_BRANCH, msg="MV Rooting failed: incorrect branch length.")
        MV.close()

        MV = open("unit_test/MV/special/output.trees")
        score = score_from_file("unit_test/MV/special/score.txt")
        branches = branch_lengths("unit_test/MV/special/branches.txt")
        score_test, MV_test, branches_test = root_trees("unit_test/MV/special/input.trees", method='MV')
        for line, line_test in zip(MV.readlines(), MV_test):
            self.assertTrue(check_two_nwk_str(line, line_test), msg="MV Rooting failed: incorrect newick string.")
        for i in score:
            self.assertTrue(score[i] - score_test[i] < EPSILON_SCORE, msg="MV Rooting failed: incorrect MV score.")
        for i in branches:
            self.assertTrue(branches[i][0] - branches_test[i][0] < EPSILON_BRANCH,
                            msg="MV Rooting failed: incorrect branch length.")
            self.assertTrue(branches[i][1] - branches_test[i][1] < EPSILON_BRANCH,
                            msg="MV Rooting failed: incorrect branch length.")
        MV.close()

    def test_RTT(self):
        """Does RTT root properly?"""
        print("Testing RTT")
        EPSILON_SCORE = 0.0001
        EPSILON_MU = 0.0001
        EPSILON_TIME = 0.001
        EPSILON_BRANCH = 0.0001

        # True Trees
        RTT = open("unit_test/RTT/true_trees/output.trees")
        score = score_from_file("unit_test/RTT/true_trees/score.txt", RTT=True)
        branches = branch_lengths("unit_test/RTT/true_trees/branches.txt")
        score_test, RTT_test, branches_test = root_trees("unit_test/RTT/true_trees/input.trees", method='RTT', timeFile="unit_test/RTT/true_trees/sampling_times.txt")
        for line, line_test in zip(RTT.readlines(), RTT_test):
            self.assertTrue(check_two_nwk_str(line, line_test), msg="RTT Rooting for true trees failed: incorrect newick string.")
        for i in score:
            self.assertTrue(0 - score_test[i][0] < EPSILON_SCORE, msg="RTT Rooting for true trees failed: incorrect RTT score.")  # RTT score
            self.assertTrue(1 - score_test[i][1] < EPSILON_MU, msg="RTT Rooting for true trees failed: incorrect mutation rate.")  # mu
            self.assertTrue(0 - score_test[i][2] < EPSILON_TIME, msg="RTT Rooting for true trees failed: incorrect initial time.")  # t0
        for i in branches:
            self.assertTrue(branches[i][0] - branches_test[i][0] < EPSILON_BRANCH, msg="RTT Rooting for true trees failed: incorrect branch length.")
            self.assertTrue(branches[i][1] - branches_test[i][1] < EPSILON_BRANCH, msg="RTT Rooting for true trees failed: incorrect branch length.")


        # True Trees with Random Initial Time
        t0 = gauss(1900, 120)
        score_test, RTT_test, branches_test = root_trees("unit_test/RTT/true_trees/input.trees", method='RTT',
                                                         timeFile="unit_test/RTT/true_trees/sampling_times.txt",t0=t0)
        for i in score:
            self.assertTrue(t0 - score_test[i][2] < EPSILON_TIME, msg="RTT Rooting for true trees with random initial time failed: incorrect initial time.")  # t0

        RTT.close()

        # Random Trees
        score = score_from_file("unit_test/RTT/random_trees/score.txt", RTT=True)
        score_test, RTT_test, branches_test = root_trees("unit_test/RTT/random_trees/input.trees", method='RTT',
                                                         timeFile="unit_test/RTT/random_trees/sampling_times.txt")
        for i in score:
            self.assertTrue(score[i][0] - score_test[i][0] < EPSILON_SCORE, msg="RTT Rooting for random trees failed: incorrect RTT score.")  # RTT score
            self.assertTrue(score[i][1] - score_test[i][1] < EPSILON_MU, msg="RTT Rooting for random trees failed: incorrect mutation rate.")  # mu
            self.assertTrue(score[i][2] - score_test[i][2] < EPSILON_TIME, msg="RTT Rooting for random trees failed: incorrect initial time.")  # t0

if __name__ == '__main__':
    unittest.main()
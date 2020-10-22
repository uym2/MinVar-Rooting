import unittest
from fastroot.Tree_extend import MPR_Tree, OGR_Tree
from treeswift import *
from fastroot.MinVar import MV00_Tree, minVAR_Base_Tree
from fastroot.RTT import RTT_Tree
from random import gauss
from os.path import dirname, realpath
from test_suite.utils import *

path=dirname(realpath(__file__))

class RootingTestCase(unittest.TestCase):
    """Tests for `FastRoot.py`."""

    def test_OG(self):
        """Does OG root properly?"""
        print("Testing OG")
        EPSILON_SCORE = 0.0001
        correct_nwk = True
        correct_score = True

        OG = open(path+"/unit_test/OG/output.trees")
        score = score_from_file(path+"/unit_test/OG/score.txt")
        #branches = branch_lengths(path+"/unit_test/OG/branches.txt")
        score_test, OG_test, branches_test = root_trees(path+"/unit_test/OG/input.trees", method='OG', OGFile=path+"/unit_test/OG/outgroups.txt")

        for line, line_test in zip(OG.readlines(), OG_test):
            correct_nwk &= check_two_nwk_str(line,line_test)
        for i in score:
            correct_score &= score[i] - score_test[i] < EPSILON_SCORE
        #for i in branches:
        #    self.assertTrue(branches[i][0] - branches_test[i][0] < 0.001)
        #    self.assertTrue(branches[i][1] - branches_test[i][1] < 0.001)

        OG.close()

        if not correct_nwk:
            print("OG Rooting failed: incorrect newick string.")
        if not correct_score:
            print("OG Rooting failed: incorrect triplet distance score.")

        self.assertTrue(correct_nwk and correct_score, msg="OG rooting failed.")

    def test_MP(self):
        """Does MP root properly?"""
        print("Testing MP")
        EPSILON_SCORE = 0.0001
        EPSILON_BRANCH = 0.0001
        correct_nwk = True
        correct_score = True
        correct_branches = True

        MP = open(path+"/unit_test/MP/output.trees")
        score = score_from_file(path+"/unit_test/MP/score.txt")
        branches = branch_lengths(path+"/unit_test/MP/branches.txt")
        score_test, MP_test, branches_test = root_trees(path+"/unit_test/MP/input.trees", method='MP')
        for line, line_test in zip(MP.readlines(), MP_test):
            correct_nwk &= check_two_nwk_str(line, line_test)
        for i in score:
            correct_score &= score[i]-score_test[i] < EPSILON_SCORE
        for i in branches:
            correct_branches &= branches[i][0]-branches_test[i][0] < EPSILON_BRANCH
            correct_branches &= branches[i][1] - branches_test[i][1] < EPSILON_BRANCH
        MP.close()

        if not correct_nwk:
            print("MP Rooting failed: incorrect newick string.")
        if not correct_score:
            print("MP Rooting failed: incorrect tree height.")
        if not correct_branches:
            print("MP Rooting failed: incorrect branch lengths.")

        self.assertTrue(correct_nwk and correct_score and correct_branches, msg="MP rooting failed.")

    def test_MP_special_trees(self):
        """Does MP root special trees properly?"""
        print("Testing MP for special trees")
        EPSILON_SCORE = 0.0001
        EPSILON_BRANCH = 0.0001
        correct_nwk = True
        correct_score = True
        correct_branches = True

        types = {1:["(balanced & ultrametric)",True,True,True], 2:["(caterpillar & ultrametric)",True,True,True], 3:["(random & not ultrametric)",True,True,True], 4:["(balanced & ultrametric)",True,True,True], 5:["(caterpillar & not ultrametric)",True,True,True]}

        MP = open(path+"/unit_test/MP/special/output.trees")
        score = score_from_file(path+"/unit_test/MP/special/score.txt")
        branches = branch_lengths(path+"/unit_test/MP/special/branches.txt")
        score_test, MP_test, branches_test = root_trees(path+"/unit_test/MP/special/input.trees", method='MP')

        for i, line in enumerate(zip(MP.readlines(), MP_test)):
            types[i+1][1] &= check_two_nwk_str(line[0], line[1])
            correct_nwk &= types[i+1][1]
        for i,j in enumerate(score):
            types[i+1][2] &= score[j] - score_test[j] < EPSILON_SCORE
            correct_score &= types[i+1][2]
        for i,j in enumerate(branches):
            types[i + 1][3] &= branches[j][0] - branches_test[j][0] < EPSILON_BRANCH
            types[i + 1][3] &= branches[j][1] - branches_test[j][1] < EPSILON_BRANCH
            correct_branches &= types[i + 1][3]

        MP.close()

        if not correct_nwk:
            errors = []
            for type in types:
                if not types[type][1]:
                    errors.append(types[type][0])
            print("MP Rooting for special trees failed: incorrect newick string for the following trees:")
            print(*[type for type in errors],sep=", ")
        if not correct_score:
            errors = []
            for type in types:
                if not types[type][2]:
                    errors.append(types[type][0])
            print("MP Rooting for special trees failed: incorrect tree height for the following trees:")
            print(*[type for type in errors], sep=", ")
        if not correct_branches:
            errors = []
            for type in types:
                if not types[type][3]:
                    errors.append(types[type][0])
            print("MP Rooting for special trees failed: incorrect branch lengths for the following trees:")
            print(*[type for type in errors], sep=", ")

        self.assertTrue(correct_nwk and correct_score and correct_branches, msg="MP rooting for special trees failed.")

    def test_MV(self):
        """Does MV root properly?"""
        print("Testing MV")
        EPSILON_SCORE = 0.0001
        EPSILON_BRANCH = 0.0001
        correct_nwk = True
        correct_score = True
        correct_branches = True

        MV = open(path+"/unit_test/MV/output.trees")
        score = score_from_file(path+"/unit_test/MV/score.txt")
        branches = branch_lengths(path+"/unit_test/MV/branches.txt")
        score_test, MV_test, branches_test = root_trees(path+"/unit_test/MV/input.trees", method='MV')

        for line, line_test in zip(MV.readlines(), MV_test):
            correct_nwk &= check_two_nwk_str(line, line_test)
        for i in score:
            correct_score &= score[i]-score_test[i] < EPSILON_SCORE
        for i in branches:
            correct_branches &= branches[i][0]-branches_test[i][0] < EPSILON_BRANCH
            correct_branches &= branches[i][1] - branches_test[i][1] < EPSILON_BRANCH

        MV.close()

        if not correct_nwk:
            print("MV Rooting failed: incorrect newick string.")
        if not correct_score:
            print("MV Rooting failed: incorrect MV score.")
        if not correct_branches:
            print("MV Rooting failed: incorrect branch lengths.")

        self.assertTrue(correct_nwk and correct_score and correct_branches, msg="MV rooting failed.")

    def test_MV_special_trees(self):
        """Does MV root special trees properly?"""
        print("Testing MV for special trees")
        EPSILON_SCORE = 0.0001
        EPSILON_BRANCH = 0.0001
        correct_nwk = True
        correct_score = True
        correct_branches = True

        types = {1:["(balanced & ultrametric)",True,True,True], 2:["(caterpillar & ultrametric)",True,True,True], 3:["(random & not ultrametric)",True,True,True], 4:["(balanced & ultrametric)",True,True,True], 5:["(caterpillar & not ultrametric)",True,True,True]}

        MV = open(path+"/unit_test/MV/special/output.trees")
        score = score_from_file(path+"/unit_test/MV/special/score.txt")
        branches = branch_lengths(path+"/unit_test/MV/special/branches.txt")
        score_test, MV_test, branches_test = root_trees(path+"/unit_test/MV/special/input.trees", method='MV')

        for i, line in enumerate(zip(MV.readlines(), MV_test)):
            types[i+1][1] &= check_two_nwk_str(line[0], line[1])
            correct_nwk &= types[i+1][1]
        for i,j in enumerate(score):
            types[i+1][2] &= score[j] - score_test[j] < EPSILON_SCORE
            correct_score &= types[i+1][2]
        for i,j in enumerate(branches):
            types[i + 1][3] &= branches[j][0] - branches_test[j][0] < EPSILON_BRANCH
            types[i + 1][3] &= branches[j][1] - branches_test[j][1] < EPSILON_BRANCH
            correct_branches &= types[i + 1][3]

        MV.close()

        if not correct_nwk:
            errors = []
            for type in types:
                if not types[type][1]:
                    errors.append(types[type][0])
            print("MV Rooting for special trees failed: incorrect newick string for the following trees:")
            print(*[type for type in errors],sep=", ")
        if not correct_score:
            errors = []
            for type in types:
                if not types[type][2]:
                    errors.append(types[type][0])
            print("MV Rooting for special trees failed: incorrect MV score for the following trees:")
            print(*[type for type in errors], sep=", ")
        if not correct_branches:
            errors = []
            for type in types:
                if not types[type][3]:
                    errors.append(types[type][0])
            print("MV Rooting for special trees failed: incorrect branch lengths for the following trees:")
            print(*[type for type in errors], sep=", ")

        self.assertTrue(correct_nwk and correct_score and correct_branches, msg="MV rooting for special trees failed.")

if __name__ == '__main__':
    unittest.main()

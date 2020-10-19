import unittest
from fastroot.Tree_extend import MPR_Tree, OGR_Tree
from treeswift import *
from fastroot.MinVar import MV00_Tree, minVAR_Base_Tree
from fastroot.RTT import RTT_Tree
from random import gauss


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


# Root input trees
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
        correct_nwk = True
        correct_score = True

        OG = open("unit_test/OG/output.trees")
        score = score_from_file("unit_test/OG/score.txt")
        #branches = branch_lengths("unit_test/OG/branches.txt")
        score_test, OG_test, branches_test = root_trees("unit_test/OG/input.trees", method='OG', OGFile="unit_test/OG/outgroups.txt")

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

        MP = open("unit_test/MP/output.trees")
        score = score_from_file("unit_test/MP/score.txt")
        branches = branch_lengths("unit_test/MP/branches.txt")
        score_test, MP_test, branches_test = root_trees("unit_test/MP/input.trees", method='MP')
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

        MP = open("unit_test/MP/special/output.trees")
        score = score_from_file("unit_test/MP/special/score.txt")
        branches = branch_lengths("unit_test/MP/special/branches.txt")
        score_test, MP_test, branches_test = root_trees("unit_test/MP/special/input.trees", method='MP')

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

        MV = open("unit_test/MV/output.trees")
        score = score_from_file("unit_test/MV/score.txt")
        branches = branch_lengths("unit_test/MV/branches.txt")
        score_test, MV_test, branches_test = root_trees("unit_test/MV/input.trees", method='MV')

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

        MV = open("unit_test/MV/special/output.trees")
        score = score_from_file("unit_test/MV/special/score.txt")
        branches = branch_lengths("unit_test/MV/special/branches.txt")
        score_test, MV_test, branches_test = root_trees("unit_test/MV/special/input.trees", method='MV')

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

    def test_RTT_true_trees(self):
        """Does RTT root true trees properly?"""
        print("Testing RTT for true trees")
        EPSILON_SCORE = 0.0001
        EPSILON_MU = 0.0001
        EPSILON_TIME = 0.001
        EPSILON_BRANCH = 0.0001
        correct_nwk = True
        correct_score = True
        correct_mu = True
        correct_t0 = True
        correct_branches = True
        correct_random_t0 = True

        # True Trees
        RTT = open("unit_test/RTT/true_trees/output.trees")
        score = score_from_file("unit_test/RTT/true_trees/score.txt", RTT=True)
        branches = branch_lengths("unit_test/RTT/true_trees/branches.txt")
        score_test, RTT_test, branches_test = root_trees("unit_test/RTT/true_trees/input.trees", method='RTT', timeFile="unit_test/RTT/true_trees/sampling_times.txt")
        for line, line_test in zip(RTT.readlines(), RTT_test):
            correct_nwk &= check_two_nwk_str(line, line_test)
        for i in score:
            correct_score &= (0 - score_test[i][0] < EPSILON_SCORE)
            correct_mu &= (1 - score_test[i][1] < EPSILON_MU)
            correct_t0 &= (0 - score_test[i][2] < EPSILON_TIME)
        for i in branches:
            correct_branches &= (branches[i][0] - branches_test[i][0] < EPSILON_BRANCH)
            correct_branches &= (branches[i][1] - branches_test[i][1] < EPSILON_BRANCH)

        # True Trees with Random Initial Time
        t0 = gauss(1900, 120)
        score_test2, RTT_test2, branches_test2 = root_trees("unit_test/RTT/true_trees/input.trees", method='RTT',
                                                         timeFile="unit_test/RTT/true_trees/sampling_times.txt",t0=t0)
        for i in score:
            correct_random_t0 &= (t0 - score_test2[i][2] < EPSILON_TIME)

        RTT.close()

        if not correct_nwk:
            print("RTT Rooting for true trees failed: incorrect newick string.")
        if not correct_score:
            print("RTT Rooting for true trees failed: incorrect RTT score.")
        if not correct_mu:
            print("RTT Rooting for true trees failed: incorrect mutation rate.")
        if not correct_t0:
            print("RTT Rooting for true trees failed: incorrect initial time.")
        if not correct_branches:
            print("RTT Rooting for true trees failed: incorrect branch lengths.")
        if not correct_random_t0:
            print("RTT Rooting for true trees with random initial time failed: incorrect initial time.")

        self.assertTrue(correct_nwk and correct_score and correct_t0 and correct_mu and correct_branches and correct_random_t0, msg="RTT rooting for true trees failed.")

        '''
        # Random Trees
        score = score_from_file("unit_test/RTT/random_trees/score.txt", RTT=True)
        score_test, RTT_test, branches_test = root_trees("unit_test/RTT/random_trees/input.trees", method='RTT',
                                                         timeFile="unit_test/RTT/random_trees/sampling_times.txt")
        for i in score:
            self.assertTrue(score[i][0] - score_test[i][0] < EPSILON_SCORE, msg="RTT Rooting for random trees failed: incorrect RTT score.")  # RTT score
            self.assertTrue(score[i][1] - score_test[i][1] < EPSILON_MU, msg="RTT Rooting for random trees failed: incorrect mutation rate.")  # mu
            self.assertTrue(score[i][2] - score_test[i][2] < EPSILON_TIME, msg="RTT Rooting for random trees failed: incorrect initial time.")  # t0
        '''

    def test_RTT_random_trees(self):
        """Does RTT root random trees properly?"""
        print("Testing RTT for random trees")
        EPSILON_SCORE = 0.0001
        EPSILON_MU = 0.0001
        EPSILON_TIME = 0.001
        correct_score = True
        correct_mu = True
        correct_t0 = True

        # Random Trees
        scores = {100: (1.654308800928493, -9.615687825811976e-06, 484861.32610241993), 500: (2.9545620705410256, 0.00022718310490997106, -27427.154466620195), 1000: (2.297664689022771, 0.017003747287065468, -365.4412497523147), 5000: (3.3270020784333463, 0.000911530692283927, -7636.059486474005)}  # RTT, mu, t0

        score_test_100, RTT_test, branches_test = root_trees("unit_test/RTT/random_trees/100/input.tre", method='RTT',
                                                         timeFile="unit_test/RTT/random_trees/100/sampling_times.txt")
        score_test_500, RTT_test, branches_test = root_trees("unit_test/RTT/random_trees/500/input.tre", method='RTT',
                                                             timeFile="unit_test/RTT/random_trees/500/sampling_times.txt")
        score_test_1000, RTT_test, branches_test = root_trees("unit_test/RTT/random_trees/1000/input.tre", method='RTT',
                                                             timeFile="unit_test/RTT/random_trees/1000/sampling_times.txt")
        score_test_5000, RTT_test, branches_test = root_trees("unit_test/RTT/random_trees/5000/input.tre", method='RTT',
                                                             timeFile="unit_test/RTT/random_trees/5000/sampling_times.txt")
        score_test = score_test_100,score_test_500,score_test_1000,score_test_5000
        for i,n in enumerate(scores):
            for j in score_test[i]:
                correct_score &= (scores[n][0] - score_test[i][j][0] < EPSILON_SCORE)
                correct_mu &= (scores[n][1] - score_test[i][j][1] < EPSILON_MU)
                correct_t0 &= (scores[n][2] - score_test[i][j][2] < EPSILON_TIME)

        if not correct_score:
            print("RTT Rooting for random trees failed: incorrect RTT score.")
        if not correct_mu:
            print("RTT Rooting for random trees failed: incorrect mutation rate.")
        if not correct_t0:
            print("RTT Rooting for random trees failed: incorrect initial time.")

        self.assertTrue(correct_score and correct_t0 and correct_mu, msg="RTT rooting for random trees failed.")


if __name__ == '__main__':
    unittest.main()
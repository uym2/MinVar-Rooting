import unittest
from random import gauss
from os.path import dirname, realpath
from fastroot_tests.utils import *

path=dirname(realpath(__file__))

class RTTTestCase(unittest.TestCase):
    """Tests for `FastRoot.py`."""

    def test_RTT_true_trees(self):
        """Does RTT root true trees properly?"""
        #print("Testing RTT for true trees")
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
        RTT = open(path+"/unit_test/RTT/true_trees/output.trees")
        score = score_from_file(path+"/unit_test/RTT/true_trees/score.txt", RTT=True)
        branches = branch_lengths(path+"/unit_test/RTT/true_trees/branches.txt")
        score_test, RTT_test, branches_test = root_trees(path+"/unit_test/RTT/true_trees/input.trees", method='RTT', timeFile=path+"/unit_test/RTT/true_trees/sampling_times.txt")
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
        score_test2, RTT_test2, branches_test2 = root_trees(path+"/unit_test/RTT/true_trees/input.trees", method='RTT',
                                                         timeFile=path+"/unit_test/RTT/true_trees/sampling_times.txt",t0=t0)
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
        score = score_from_file(path+"/unit_test/RTT/random_trees/score.txt", RTT=True)
        score_test, RTT_test, branches_test = root_trees(path+"/unit_test/RTT/random_trees/input.trees", method='RTT',
                                                         timeFile=path+"/unit_test/RTT/random_trees/sampling_times.txt")
        for i in score:
            self.assertTrue(score[i][0] - score_test[i][0] < EPSILON_SCORE, msg="RTT Rooting for random trees failed: incorrect RTT score.")  # RTT score
            self.assertTrue(score[i][1] - score_test[i][1] < EPSILON_MU, msg="RTT Rooting for random trees failed: incorrect mutation rate.")  # mu
            self.assertTrue(score[i][2] - score_test[i][2] < EPSILON_TIME, msg="RTT Rooting for random trees failed: incorrect initial time.")  # t0
        '''

    def test_RTT_random_trees(self):
        """Does RTT root random trees properly?"""
        #print("Testing RTT for random trees")
        EPSILON_SCORE = 0.0001
        EPSILON_MU = 0.0001
        EPSILON_TIME = 0.001
        correct_score = True
        correct_mu = True
        correct_t0 = True

        # Random Trees
        scores = {100: (1.654308800928493, -9.615687825811976e-06, 484861.32610241993), 500: (2.9545620705410256, 0.00022718310490997106, -27427.154466620195), 1000: (2.297664689022771, 0.017003747287065468, -365.4412497523147), 5000: (3.3270020784333463, 0.000911530692283927, -7636.059486474005)}  # RTT, mu, t0

        score_test_100, RTT_test, branches_test = root_trees(path+"/unit_test/RTT/random_trees/100/input.tre", method='RTT',
                                                         timeFile=path+"/unit_test/RTT/random_trees/100/sampling_times.txt")
        score_test_500, RTT_test, branches_test = root_trees(path+"/unit_test/RTT/random_trees/500/input.tre", method='RTT',
                                                             timeFile=path+"/unit_test/RTT/random_trees/500/sampling_times.txt")
        score_test_1000, RTT_test, branches_test = root_trees(path+"/unit_test/RTT/random_trees/1000/input.tre", method='RTT',
                                                             timeFile=path+"/unit_test/RTT/random_trees/1000/sampling_times.txt")
        score_test_5000, RTT_test, branches_test = root_trees(path+"/unit_test/RTT/random_trees/5000/input.tre", method='RTT',
                                                             timeFile=path+"/unit_test/RTT/random_trees/5000/sampling_times.txt")
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


#if __name__ == '__main__':
#    unittest.main()

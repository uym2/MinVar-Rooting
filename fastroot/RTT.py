from fastroot.Tree_extend import *
from fastroot.quadprog_solvers import *
from sys import stderr

EPSILON = 1e-5

class RTT_Tree(Tree_extend):
    # supportive base class to implement VAR-reroot, hence the name
    def __init__(self,smplTimes, ddpTree=None, tree_file=None,schema="newick",logger_id=1,logger_stream=stderr):
        super(RTT_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("RTT_Tree_" + str(logger_id),myStream=logger_stream)
        self.smplTimes = smplTimes
        self.reset()

    def reset(self):
        self.RTT = None
        self.opt_root = self.ddpTree.root
        self.opt_y = 0
        self.opt_x = 0
        self.opt_mu = 0
        self.tmin = min(self.smplTimes.values())

    def Node_init(self, node, nleaf=1, SDI=0, SD=0, ST=0, SDT=0, SSD=0):
        node.SDI = SDI
        node.SD = SD
        node.nleaf = nleaf
        node.ST = ST
        node.SDT = SDT
        node.SSD = SSD

    def Opt_function(self, node, SST, deltaT, deltaD, SDT, SSD, ST, SD):
        n = self.total_leaves
        a, b, c, d, e, f = n, SST, (-2 * deltaT), (2 * deltaD), (-2 * SDT), SSD
        k, m, r = 2*(n-2*node.nleaf), -2*ST, 2*SD

        tmin = self.tmin

        # use quadprog to compute mu_star, y_star, and x_star
        P = array([[a,k/2,c/2.],[k/2,n,m/2],[c/2,m/2,b]])
        q = array([d/2.,r/2,e/2])
        G = array([[-1.,0.,0.], [0.,0.,-1.], [1.,0.,0.],[0.,1.,-tmin]])
        h = array([0., EPSILON, node.edge_length,0]).reshape((4,))
        solution = cvxopt_solve_qp(P,q,G,h)
        x_star = solution[0]
        y_star = solution[1]
        mu_star = solution[2]
        curr_RTT = a*x_star*x_star + b*mu_star*mu_star + c*x_star*mu_star + d*x_star + e*mu_star + f + n*y_star*y_star + k*x_star*y_star + m*mu_star*y_star + r*y_star
        
        if self.RTT is None or (curr_RTT - self.RTT < -EPSILON):
            self.RTT = curr_RTT
            self.opt_root = node
            self.opt_x = node.edge_length - x_star
            self.opt_y = y_star
            self.opt_mu = mu_star

    def bUp_update(self, node):
        if node.is_leaf():
            node.nleaf = 1
            node.SDI = 0
            node.ST = self.smplTimes[node.label]
        else:
            node.nleaf = 0
            node.SDI = 0
            node.ST = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.SDI += child.SDI + child.nleaf * child.edge_length
                node.ST += child.ST

    def Update_var(self, child, node, edge_length):
        SST = self.SST
        deltaT = self.ddpTree.root.ST - 2 * child.ST
        deltaD = -2 * child.nleaf * edge_length - 2 * child.SDI + node.SD
        SDT = node.SDT
        SSD = node.SSD
        ST = self.ST
        SD = node.SD
        return SST, deltaT, deltaD, SDT, SSD, ST, SD

    def tDown_update(self, node, opt_function):
        for child in node.child_nodes():
            child.SD = node.SD + (self.total_leaves - 2 * child.nleaf) * child.edge_length
            child.SDT = node.SDT + child.edge_length * (self.ddpTree.root.ST - 2 * child.ST)
            child.SSD = node.SSD + (self.total_leaves - 4 * child.nleaf) * (child.edge_length ** 2) + 2 * (node.SD - 2 * child.SDI) * child.edge_length
            SST, deltaT, deltaD, SDT, SSD, ST, SD = self.Update_var(child, node, child.edge_length)
            opt_function(child, SST, deltaT, deltaD, SDT, SSD, ST, SD)

    def prepare_root(self):
        root = self.get_root()
        root.SD = root.SDI
        self.total_leaves = root.nleaf
        self.ST = root.ST
        self.ddpTree.root.droot = 0
        self.ddpTree.root.troot = 0
        root.SD, root.SSD, root.SDT, self.SST = 0, 0, 0, 0
        for v in self.ddpTree.traverse_preorder():
            if not v.is_root():
                # must have defined edge lengths
                v.droot = v.parent.droot + v.edge_length
                if v.is_leaf():
                    root.SSD += (v.droot ** 2)
                    self.SST += (self.smplTimes[v.label] ** 2)
                    root.SD += v.droot
                    root.SDT += (v.droot * self.smplTimes[v.label])

    def opt_score(self):
        return self.RTT

    def report_score(self):
        return "RTT=" + str(self.opt_score()/self.total_leaves) + "\tmu=" + str(self.opt_mu) +  "\tt0=" + str(self.opt_y/self.opt_mu)

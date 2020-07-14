from Tree_extend import *

class RTT_Tree(Tree_extend):
    # supportive base class to implement VAR-reroot, hence the name
    def __init__(self,smplTimes, ddpTree=None, tree_file=None,schema="newick"):
        super(RTT_Tree, self).__init__(ddpTree, tree_file, schema)
        self.smplTimes = smplTimes
        self.reset()

    def reset(self):
        self.RTT = None
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, nleaf=1, SDI=0, SD=0, var=-1, ST=0, SDT=0, SSD=0):
        node.SDI = SDI
        node.SD = SD
        node.nleaf = nleaf
        #node.var = var # is this necessary? should I take out?
        node.ST = ST
        node.SDT = SDT
        node.SSD = SSD

    '''
        def Opt_function(self, node, a, b, c):
            print("Abstract method! Should never be called")
    '''

    def Opt_function(self, node, SST, deltaT, deltaD, SDT, SSD):
        n = self.total_leaves
        m = (SDT - deltaT * deltaD / n) / (SST - deltaT * deltaT / n)
        x = (deltaT * m - deltaD) / n
        if x >= 0 and x <= node.edge_length:
            curr_RTT = n * x * x + SST * m * m - 2 * deltaT * x * m + 2 * deltaD * x - 2 * SDT * m + SSD
            if self.RTT is None or curr_RTT < self.RTT:
                self.RTT = curr_RTT
                self.opt_root = node
                self.opt_x = node.edge_length - x

    '''
    def compute_dRoot_VAR(self):################
        cumm = {'ssq': 0, 'sum': 0}

        def compute_dRoot(node, cumm_l):
            if node.is_leaf():
                cumm['ssq'] += cumm_l ** 2
                cumm['sum'] += cumm_l
            else:
                for child in node.child_nodes():
                    compute_dRoot(child, cumm_l + child.edge_length)

        compute_dRoot(self.get_root(), 0)
        N = self.get_root().nleaf
        root_var = cumm['ssq'] / N - (cumm['sum'] / N) ** 2
        self.get_root().var = root_var
    '''

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
        return SST, deltaT, deltaD, SDT, SSD

    def tDown_update(self, node, opt_function):
        for child in node.child_nodes():
            child.SD = node.SD + (self.total_leaves - 2 * child.nleaf) * child.edge_length
            child.SDT = node.SDT + child.edge_length * (self.ddpTree.root.ST - 2 * child.ST)
            child.SSD = node.SSD + (self.total_leaves - 4 * child.nleaf) * (child.edge_length ** 2) + 2 * (node.SD - 2 * child.SDI) * child.edge_length
            SST, deltaT, deltaD, SDT, SSD = self.Update_var(child, node, child.edge_length)
            opt_function(child, SST, deltaT, deltaD, SDT, SSD)

    def prepare_root(self):
        root = self.get_root()
        root.SD = root.SDI
        #self.compute_dRoot_VAR() ########
        self.total_leaves = root.nleaf
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
        print("SD:",root.SD,"  SSD:",root.SSD,"  SDT:", root.SDT, "  SST:", self.SST, "  ST:", root.ST)
        # function works for sample1 and sample2

    def opt_score(self):
        return self.RTT

    def report_score(self):
        return "RTT score: " + str(self.opt_score())
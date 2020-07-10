from Tree_extend import *

class RTT_Tree(Tree_extend):
    # supportive base class to implement VAR-reroot, hence the name
    def __init__(self,smplTimes, ddpTree=None, tree_file=None,schema="newick"):
        super(RTT_Tree, self).__init__(ddpTree, tree_file, schema)
        self.smplTimes = smplTimes
        self.reset()

    def reset(self):
        self.RTT = None
        self.opt_root = self.ddpTree.root #########
        self.opt_x = 0

    def Node_init(self, node, nleaf=1, sum_in=0, sum_total=0, var=-1):
        node.sum_in = sum_in
        node.sum_total = sum_total
        node.nleaf = nleaf
        node.var = var

    def Opt_function(self, node, a, b, c):
        print("Abstract method! Should never be called")

    def compute_dRoot_VAR(self):
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

    def bUp_update(self, node):
        if node.is_leaf():
            node.nleaf = 1
            node.sum_in = 0
            node.ST = self.smplTimes[node.label]
        else:
            node.nleaf = 0
            node.sum_in = 0
            node.ST = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.sum_in += child.sum_in + child.nleaf * child.edge_length
                node.ST += child.ST

    def Update_var(self, child, node, edge_length):
        alpha = 2 * (node.sum_total - 2 * (child.sum_in + child.nleaf * edge_length)) / self.total_leaves
        beta = 1 - 2 * float(child.nleaf) / self.total_leaves
        a = 1 - beta * beta
        b = alpha - 2 * node.sum_total * beta / self.total_leaves
        c = node.var
        child.var = a * edge_length * edge_length + b * edge_length + c
        return a, b, c

    def tDown_update(self, node, opt_function):
        for child in node.child_nodes():
            child.sum_total = node.sum_total + (self.total_leaves - 2 * child.nleaf) * child.edge_length
            a, b, c = self.Update_var(child, node, child.edge_length)
            opt_function(child, a, b, c)

    def prepare_root(self):
        root = self.get_root()
        root.sum_total = root.sum_in
        self.compute_dRoot_VAR()
        self.total_leaves = root.nleaf

    def opt_score(self):
        return self.RTT

    def report_score(self):
        return "RTT score: " + str(self.opt_score())
try:
	from dendropy4 import Tree,Node
except:
	from dendropy import Tree,Node

class Tree_extend(object):
	def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
		if tree_file:
			self.ddpTree = Tree.get_from_path(tree_file,schema)
		else:
			self.ddpTree = ddpTree
		self.Tree_records = Tree_records

	def New_record(self,old_label=None):
		print "Abstract method!"
		return Node_record(old_label=old_label)

	def Bottomup_update(self):
		i = 0
		for node in self.ddpTree.postorder_node_iter():
		# store old label
			#node_record = Node_record(old_label=node.label)
			node_record = self.New_record(old_label=node.label)
		# assign new label
			node.label = i
			node_record.Bottomup_update(node,self.Tree_records)
			self.Tree_records.append(node_record)
			i = i+1
	
	def Topdown_update(self):
		for node in self.ddpTree.preorder_node_iter():
			self.Tree_records[node.label].Topdown_update(node,self.Tree_records,self.Opt_function)


	def Reroot(self):
		print "Abstract method!"

	def Opt_function(self,node):
		print "Abstract method!"


	def tree_as_newick(self,outfile=None,restore_label=False):
	# dendropy's method to write newick seems to have problem ...
		if outfile:
			outstream = open(outfile,'w')
		else:
			outstream = stdout
	
		self.__write_newick(self.ddpTree.seed_node,outstream,restore_label)
		outstream.write(";")

		if outfile:
			outstream.close()	

	def __write_newick(self,node,outstream,restore_label=False):
		if node.is_leaf():
				outstream.write(node.taxon.label)
		else:
			outstream.write('(')
			is_first_child = True
			for child in node.child_node_iter():
				if is_first_child:
					is_first_child = False
				else:
					outstream.write(',')
				self.__write_newick(child,outstream,restore_label=restore_label)
			outstream.write(')')
		if not node.is_leaf():
			if restore_label: 
				not node.label is None and not self.Tree_records[node.label].old_label is None and outstream.write(str(self.Tree_records[node.label].old_label))
			else:
				not node.label is None and outstream.write(str(node.label))
		if node.edge_length:
			outstream.write(":"+str(node.edge_length))

	def reroot_at_edge(self,edge,length1,length2,new_root=None):
	# the method provided by dendropy DOESN'T seem to work ...
		head = edge.head_node
		tail = edge.tail_node
	
		if not new_root:
			new_root = Node()		

		tail.remove_child(head)
		
		new_root.add_child(head)
		head.edge_length=length2
		head.edge.length = length2

		p = tail.parent_node
		l = tail.edge_length

		new_root.add_child(tail)
		tail.edge_length=length1
		tail.edge.length=length1

		if tail.label == self.ddpTree.seed_node.label:
			head = new_root

		while tail.label != self.ddpTree.seed_node.label:
			head = tail
			tail = p
			p = tail.parent_node

			l1 = tail.edge_length
			tail.remove_child(head)

			head.add_child(tail)
			tail.edge.length=l
			tail.edge_length=l
			l = l1
			
		# out of while loop: tail IS now tree.seed_node
		if tail.num_child_nodes() < 2:
			# merge the 2 branches of the old root and adjust the branch length
			sis = [child for child in tail.child_node_iter()][0]
			l = sis.edge_length
			tail.remove_child(sis)	
			head.add_child(sis)
			sis.edge_length=l+tail.edge_length
			head.remove_child(tail)
	
		self.ddpTree.seed_node = new_root

class MPR_Tree(Tree_extend):
# supportive class to implement midpoint-reroot (mpr = mid point reroot, hence the name)
	def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
		if tree_file:
			self.ddpTree = Tree.get_from_path(tree_file,schema)
		else:
			self.ddpTree = ddpTree
		self.Tree_records = Tree_records
		self.max_distance = -1
		self.opt_root = self.ddpTree.seed_node
		self.opt_x = -1

	def New_record(self,old_label=None):
		return MPR_Node_record(old_label=old_label)

	def Opt_function(self,node):
		m = max(self.Tree_records[node.label].max_in) 
		curr_max_distance = m + self.Tree_records[node.label].max_out
		x = (self.Tree_records[node.label].max_out - m)/2
		if curr_max_distance > self.max_distance and x >= 0 and x <= node.edge_length:
			self.max_distance = curr_max_distance
			self.opt_x = x
			self.opt_root = node

	def Reroot(self):
		self.Bottomup_update()
		self.Topdown_update()
		self.reroot_at_edge(self.opt_root.edge,self.opt_root.edge_length-self.opt_x,self.opt_x)

class Node_record(object):
	def __init__(self,old_label=None):
		self.old_label = old_label
	
	def Bottomup_update(self,node,Tree_records):
		print "Just an abstract method! You should never see this message. Otherwise please check your code!"

class MPR_Node_record(Node_record):
# supportive class to implement midpoint-reroot (mpr = mid point reroot, hence the name)
	def __init__(self,old_label=None,max_in=[0,0],max_out=-1):
		self.old_label=old_label
		self.max_in = max_in
		self.max_out = max_out

	def Bottomup_update(self,node,Tree_records):
		if not node.is_leaf():
			self.max_in=[]
			for child in node.child_node_iter():
				self.max_in.append(max(Tree_records[child.label].max_in) + child.edge_length)	
	
	def Topdown_update(self,node,Tree_records,opt_function):
		child_idx = 0
		for child in node.child_node_iter():	
			Tree_records[child.label].max_out = max([self.max_out]+[self.max_in[k] for k in range(len(self.max_in)) if k != child_idx])+child.edge_length
			opt_function(child)
			child_idx = child_idx+1


class STD_Node_record(Node_record):
# supportive class to implement std-reroot, hence the name
	def __init__(self,sum_in=0,sum_out=0):
		self.sum_in = sum_in
		self.sum_out = sum_out1

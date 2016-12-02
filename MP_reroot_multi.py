from dendropy4 import Tree,Node
from sys import argv,stdout

tree_file = argv[1]

class record(object):
	def __init__(self,max_in=[0,0],max_out=0):
		self.max_in = max_in
		self.max_out = max_out

def tree_as_newick(tree,outfile=None):
	if outfile:
		outstream = open(outfile,'w')
	else:
		outstream = stdout
	
	_write_newick(tree.seed_node,outstream)
	outstream.write(";")

	if outfile:
		outstream.close()	

def _write_newick(node,outstream):
	if node.is_leaf():
			outstream.write(str(node.taxon))
	else:
		outstream.write('(')
		is_first_child = True
		for child in node.child_node_iter():
			if is_first_child:
				is_first_child = False
			else:
				outstream.write(',')
			_write_newick(child,outstream)
		outstream.write(')')
	if not node.is_leaf() and node.label:
		outstream.write(str(node.label))
	if node.edge_length:
		outstream.write(":"+str(node.edge_length))

def deroot_tree(tree,parent_side):
	root = tree.seed_node
	left_child,right_child = [child for child in node.child_node_iter()]
	combined_edge = left_child.edge_length + right_child.edge_length
	root.remove_child(left_child)
	root.remove_child(right_child)
	if parent_side == 'r':
		right_child.add_child(left_child,edge_length=combined_edge)
	else:
		left_child.add_child(right_child,edge_length=combined_edge)
	
def reroot_at_edge(tree,edge,length1,length2,new_root=None):
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
	is_root_edge = False
	if tail.label == tree.seed_node.label:
		head = new_root
        	is_root_edge = True
		l = tail.edge_length

	while tail.label != tree.seed_node.label:
		head = tail
		tail = p
		p = tail.parent_node
		
		if tail.label == tree.seed_node.label:
			break
		l1 = tail.edge_length

		tail.remove_child(head)

		p = tail.parent_node
		l1 = tail.edge_length

		head.add_child(tail)
		tail.edge.length=l
		tail.edge_length=l
		
		l = l1
		
	# out of while loop: tail IS now tree.seed_node
	sis = [child for child in tail.child_node_iter() if child.label != head.label][0]
	l1 = sis.edge_length
	if not is_root_edge:
		tail.remove_child(head)
	tail.remove_child(sis)	
	head.add_child(sis)
	sis.edge.length=l+l1
	sis.edge_length=l+l1	

	tree.seed_node = new_root

a_tree = Tree.get_from_path(tree_file,"newick")

i = 0
record_arr = []


for node in a_tree.postorder_node_iter():
	node.label = i
	if node.is_leaf():
		node_record = record()
		record_arr.append(node_record)
	else:
		max_in = []
		for child in node.child_node_iter():
			max_in.append(max(record_arr[child.label].max_in) + child.edge_length)	
		node_record = record(max_in=max_in)
		record_arr.append(node_record)
	i = i+1

record_arr[a_tree.seed_node.label].max_outside = 0
max_distance = 0
opt_root = a_tree.seed_node
opt_x = 0

for node in a_tree.preorder_node_iter():
	child_idx = 0
	node_record = record_arr[node.label]

	for child in node.child_node_iter():
		record_arr[child.label].max_out = max([node_record.max_out]+[node_record.max_in[k] for k in range(len(node_record.max_in)) if k != child_idx])+child.edge_length
		m = max(record_arr[child.label].max_in) 
		curr_max_distance = m + record_arr[child.label].max_out
		x = (record_arr[child.label].max_out - m)/2
		if curr_max_distance > max_distance and x >= 0 and x <= child.edge_length:
			max_distance = curr_max_distance
			opt_x = x
			opt_root = child
		child_idx = child_idx+1

opt_idx = opt_root.label
reroot_at_edge(a_tree,opt_root.edge,opt_root.edge_length-opt_x,opt_x)

tree_as_newick(a_tree,outfile=tree_file.split(".tre")[0]+"_rooted.tre")

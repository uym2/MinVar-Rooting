from dendropy4 import Tree

a_tree = Tree.get_from_path("test.tre","newick")

i = 0
max_left = []
max_right = []

for node in a_tree.postorder_node_iter():
	node.label = i
	if node.is_leaf():
		max_left.append(0)
		max_right.append(0)
	else:
		left_child,right_child = [child for child in node.child_node_iter()]
		max_left.append(max(max_left[left_child.label],max_right[left_child.label])+left_child.edge_length)
		max_right.append(max(max_left[right_child.label],max_right[right_child.label])+right_child.edge_length)
		
	i = i+1

for node in a_tree.postorder_node_iter():
	print str(node.label) + " " + str(node.edge_length) + " " + str(max_left[node.label]) + " " + str(max_right[node.label])
	#print max_left[edge.label]
	#print max_right[edge.label]


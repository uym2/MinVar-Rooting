# Usage: python d2NewRoot.py <d_head> <d_tail> <edge_length> <x>

from sys import argv

d_head=float(argv[1])
d_tail=float(argv[2])
l=float(argv[3])
x=float(argv[4])

epsilon = 10**-5

if ( abs(d_head + d_tail - l) < epsilon):
	print "I am in between!"
	print "d = " + str(abs(d_head-x))
elif ( abs(d_head - d_tail - l) < epsilon):
	print "I am near tail!"
	print "d = " + str(abs(d_head-x))
elif ( abs(d_tail - d_head - l) < epsilon):
	print "I am near head!"
	print "d = " + str(d_head+x)
else:
	print "I am at nowhere! :("


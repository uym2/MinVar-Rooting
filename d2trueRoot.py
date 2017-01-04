#! /usr/bin/env python

# Usage: python d2NewRoot.py <d_head> <d_tail> <edge_length> <x>

from sys import argv

d_head=float(argv[1])
l=float(argv[3])
x=float(argv[4])

if argv[2] == "sum":
	d_tail=d_head+l
elif argv[2] == "sub":
	d_tail=d_head-l
else:
	d_tail=float(argv[2])

epsilon = 10**-5



if (l == 0):
	print("I am right there! :)")
	print("d = 0")
elif ( abs(d_head + d_tail - l) < epsilon):
	print("I am in between!")
	print("d = " + str(abs(d_head-x)))
elif ( abs(d_head - d_tail - l) < epsilon):
	print("I am near tail!")
	print("d = " + str(abs(d_head-x)))
elif ( abs(d_tail - d_head - l) < epsilon):
	print("I am near head!")
	print("d = " + str(d_head+x))
else:
	print("I am at nowhere! :(")


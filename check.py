class A:
	def check(self):
		print "check A!"

class B(A):
	def check_inherit(self):
		B.check()

b = B()
b.check()

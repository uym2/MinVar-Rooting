import numpy
import quadprog
from numpy import *

def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
# minimize x^T*P*x + q^T*x s.t. Gx <= h and Ax = b
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -numpy.vstack([A, G]).T
        qp_b = -numpy.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]

P = array([[2.,0.],[0,2.]])
q = array([-8.,-6.])
G = array([[-1.,0.],[0.,-1.],[1.,1.]])
h = array([0.,0.,5.]).reshape((3,))

print(quadprog_solve_qp(P, q, G, h))    

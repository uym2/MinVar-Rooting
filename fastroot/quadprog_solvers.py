import numpy
#import quadprog
import cvxopt
from numpy import *
import logging
from sys import stdout

logger = logging.getLogger("quadprog_solvers")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False

'''
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
'''

def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [cvxopt.matrix(P), cvxopt.matrix(q)]
    if G is not None:
        args.extend([cvxopt.matrix(G), cvxopt.matrix(h)])
        if A is not None:
            args.extend([cvxopt.matrix(A), cvxopt.matrix(b)])
    sol = cvxopt.solvers.qp(*args,options={'show_progress':False})
    if 'optimal' not in sol['status']:
        return None
    return numpy.array(sol['x']).reshape((P.shape[1],))

if __name__ == "__main__":
    P = array([[2.,0.],[0,2.]])
    q = array([-8.,-6.])
    G = array([[-1.,0.],[0.,-1.],[1.,1.]])
    h = array([0.,0.,5.]).reshape((3,))
    #logger.info(quadprog_solve_qp(P,q,G,h))
    logger.info(cvxopt_solve_qp(P,q,G,h))

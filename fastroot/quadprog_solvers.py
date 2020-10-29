import numpy
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

def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, maxIter=1000):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [cvxopt.matrix(P), cvxopt.matrix(q)]
    if G is not None:
        args.extend([cvxopt.matrix(G), cvxopt.matrix(h)])
        if A is not None:
            args.extend([cvxopt.matrix(A), cvxopt.matrix(b)])
    sol = cvxopt.solvers.qp(*args,options={'show_progress':False,'maxiters':maxIter})
    if 'optimal' not in sol['status']:
        if "unknown" in sol['status']:
            logger.warning("Couldn't find optimal solution on one branch. Perhaps due to maximum iterations exceeded. Consider increasing the maximum iterations via -x.")
        else:
            logger.warning("Couldn't find optimal solution on one branch. Solution status: " + sol['status'])
	#return None
    return numpy.array(sol['x']).reshape((P.shape[1],))

if __name__ == "__main__":
    P = array([[2.,0.],[0,2.]])
    q = array([-8.,-6.])
    G = array([[-1.,0.],[0.,-1.],[1.,1.]])
    h = array([0.,0.,5.]).reshape((3,))
    #logger.info(quadprog_solve_qp(P,q,G,h))
    logger.info(cvxopt_solve_qp(P,q,G,h))

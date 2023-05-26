##
##  PhD Works
##  Quite a simple simulation
##
##  Carlos Persiani
##  Henrique Garcia
##

import numpy as np
from math import pi
from numpy.linalg import matrix_rank, eig
import matplotlib.pyplot as plt
import mRPI
import scipy
from scipy.spatial import ConvexHull
from scipy.signal import place_poles
from control import ctrb
from sim_2d import Dynamics, disturbances


def dlqr(A,B,Q,R):

    """Solve the discrete time lqr controller.
     
    x[k+1] = A x[k] + B u[k]
     
    cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
    """
    #ref Bertsekas, p.151
     
    #first, try to solve the ricatti equation
    X = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
     
    #compute the LQR gain
    K = np.matrix(scipy.linalg.inv(B.T*X*B+R)*(B.T*X*A))
     
    #eigVals, eigVecs = scipy.linalg.eig(A-B*K)
     
    return K

def define_bench() -> Dynamics:

    h = 0.01

    Ac = np.array([[0, 1],
                   [0, -0.010/4.595*1e-3]])
                
    Bc = np.array([[0, 0],
                   [0.25/4.595*1e-3, -0.25/4.595*1e-3]])

    Ad = h*Ac + np.eye(2)
    Bd = h*Bc

    bench = Dynamics()

    bench.Aol = Ad
    bench.B = Bd
    bench.K = dlqr(Ad, Bd, 1e11*np.diag([1,0.00001]), 1e-5*np.eye(2))
    bench.Acl = bench.Aol - np.dot(bench.B, bench.K)

    return bench

def main():

    bench = define_bench()

    # Disturbances ranges (a n-uple for each system dimension)
    w = np.array([[-5*pi/180,5*pi/180],[-0.2,0.2]])

    # First RPI guess set vertices
    I = np.array([[-3.14,3.14],[-2,2]])

    ## initial state
    x = np.array([0,0])
    u = np.array([0,0])

    phi = mRPI.mRPI(np.array(bench.Acl), w, I, 100)

    plt.figure()
    for simplex in phi.simplices:
        plt.plot(phi.points[simplex, 0], phi.points[simplex, 1], 'k-')
    plt.grid(True)

    X1 = [x[0]]
    X2 = [x[1]]

    for i in range(2000):

        x = bench.update(x, u, w)
        u = bench.control(x)

        X1.append(x[0])
        X2.append(x[1])

    plt.plot(X1,X2,'*b-')
    plt.xlabel(r'$x_{1}$ State')
    plt.ylabel(r'$x_{2}$ State')

    plt.show()

if __name__ == '__main__':
    main()

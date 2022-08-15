##
##  PhD Works
##  Quite a simple simulation
##
##  Carlos Persiani
##  Henrique Garcia
##

import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt
import mRPI
from scipy.spatial import ConvexHull

class Dynamics(object):

    def __init__(self):
        
        self.Aol = np.array([[1.000, 0.650],
                             [0.000, 1.000]])

        self.B = np.array([0.211,
                           0.650])

        self.K = - np.array([0.575, 1.217])

        self.Acl = np.array([[0.879, 0.393],
                             [-0.374, 0.209]])

    def update(self,x:np.array, u:np.array, w:np.array) -> np.array:

        return np.dot(self.Aol, x) + np.dot(self.B,u) + disturbances(*w)

    def control(self, x:np.array) -> np.array:

        return np.dot(self.K, x)

def disturbances(bound1:np.array, bound2:np.array) -> np.array:

    sigma1 = abs(bound1[0])/5           ## p[Out of domain] < 0.01%
    sigma2 = abs(bound2[0])/5           ## p[Out of domain] < 0.01%

    rand1 = sigma1*np.random.randn()
    rand2 = sigma2*np.random.randn()

    if rand1 > bound1[0]: rand1 = bound1[0]
    if rand1 < bound1[1]: rand1 = bound1[1]
    if rand2 > bound2[0]: rand2 = bound2[0]
    if rand2 < bound2[1]: rand2 = bound2[1]

    return np.array([rand1, rand2])

def main():

    # Dynamic controled matrix (eig must lie inside the unit circle)
    Acl = np.array([[0.879, 0.393],
                    [-0.374, 0.209]])

    #print(abs(eig(Acl)[0][0]))
    #print(abs(eig(Acl)[0][1]))

    # Disturbances ranges (a n-uple for each system dimension)
    w = np.array([[0.211,-0.211],[0.65,-0.65]])
    w = np.array([[0.1,-0.1],[0.3,-0.3]])
    w1 = np.array([[0.07,-0.07],[0.35,-0.35]])

    # First RPI guess set vertices
    I = np.array([[5,-5],[5,-5]])

    # Just to see if the selected point lies inside the calculated mRPI
    # test_point = np.array([0.5,1.075])

    # Initial position
    x = np.array([0,0])
    u = 0

    phi = mRPI.mRPI(Acl, w, I, 50)

    ES = mRPI.Event_Set(np.array(phi.equations), Acl, w)

    plt.figure()

    for simplex in phi.simplices:
        plt.plot(phi.points[simplex, 0], phi.points[simplex, 1], 'k-')
    plt.grid(True)

    for simplex in ES.simplices:
        plt.plot(ES.points[simplex, 0], ES.points[simplex, 1], 'r--')
    plt.grid(True)

    Dyn = Dynamics()

    X1 = [x[0]]
    X2 = [x[1]]
    U = [u]

    count= 0

    for i in range(1000):

        x = Dyn.update(x, u, w)
        
        if not mRPI.whether_inside(np.array(ES.equations),x): 
            u = Dyn.control(x)
            count += 1
        else: 
            u = 0

        X1.append(x[0])
        X2.append(x[1])
        U.append(u)

    print(count)

    plt.plot(X1,X2,'*b-')
    plt.xlabel(r'$x_{1}$ State')
    plt.ylabel(r'$x_{2}$ State')
    plt.axis("Equal")

    plt.figure()
    plt.plot(X1, 'k')
    plt.xlabel('Time Scale [ ]')
    plt.ylabel(r'$x_{1}$ State')
    plt.grid(True)

    plt.figure()
    plt.plot(X2, 'k')
    plt.xlabel('Time Scale [ ]')
    plt.ylabel(r'$x_{2}$ State')
    plt.grid(True)

    plt.figure()
    plt.plot(U, 'k-*')
    plt.xlabel('Time Scale [ ]')
    plt.ylabel(r'$u$ Input')
    plt.grid(True)

    plt.show()

if __name__ == '__main__':

    main()



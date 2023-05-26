##
## São Carlos
## Aug, 2022
##
## Carlos Persiani
## Henrique Garcia
##

import scipy
from xmlrpc.client import boolean
from scipy.spatial import ConvexHull, convex_hull_plot_2d, HalfspaceIntersection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import itertools
from typing import List, Text, Union

SetOfPoints = np.array
Matrix = np.array

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

def E(w:np.array) -> ConvexHull:

    return ConvexHull(list(itertools.product(*w)))

def Event_Set(HP:np.array, Aol:np.array, w:np.array) -> ConvexHull:
    
    w_plus = np.array([[w[0,0]],[w[1,0]]])
    w_minus = np.array([[w[0,1]],[w[1,1]]])

    H = HP[:,:-1]
    P = np.reshape(HP[:,-1],(H.shape[0],1))

    P_plus = P - np.dot(H, w_plus)
    P_minus = P - np.dot(H, w_minus)

    HA = np.dot(H, Aol)

    P_new = np.concatenate((P_plus,P_minus),axis=0)
    H_new = np.concatenate((HA,HA),axis=0)

    halfspaces = np.concatenate((H_new,P_new), axis=1)

    return h2v_representation(halfspaces, np.array([0.0,0.0]))

def Event_Set_2(HP:np.array, Aol:np.array, w:np.array) -> ConvexHull:
    
    w_plus = np.array([[w[0,0]],[w[1,0]]])
    w_minus = np.array([[w[0,1]],[w[1,1]]])

    H = HP[:,:-1]
    P = np.reshape(HP[:,-1],(H.shape[0],1))

    E_star = Aol + np.eye(Aol.shape[0])

    P_plus = P - np.dot(np.dot(H,E_star), w_plus)
    P_minus = P - np.dot(np.dot(H,E_star), w_minus)

    HAA = np.dot(np.dot(H, Aol), Aol)

    P_new = np.concatenate((P_plus,P_minus),axis=0)
    H_new = np.concatenate((HAA,HAA),axis=0)

    halfspaces = np.concatenate((H_new,P_new), axis=1)

    return h2v_representation(halfspaces, np.array([0,0]))

<<<<<<< HEAD
def Event_Set_3(HP:np.array, Aol:np.array, w:np.array) -> ConvexHull:
    
    w_plus = np.array([[w[0,0]],[w[1,0]]])
    w_minus = np.array([[w[0,1]],[w[1,1]]])

    H = HP[:,:-1]
    P = np.reshape(HP[:,-1],(H.shape[0],1))

    E_star = np.dot(Aol,Aol) + Aol + np.eye(Aol.shape[0])

    P_plus = P - np.dot(np.dot(H,E_star), w_plus)
    P_minus = P - np.dot(np.dot(H,E_star), w_minus)

    HAA = np.dot(np.dot(H, np.dot(Aol,Aol)), Aol)

    P_new = np.concatenate((P_plus,P_minus),axis=0)
    H_new = np.concatenate((HAA,HAA),axis=0)

    halfspaces = np.concatenate((H_new,P_new), axis=1)

    return h2v_representation(halfspaces, np.array([0,0]))
=======
>>>>>>> 81574d53993bf5def78f180de44c8dd2734b427b

def h2v_representation(HP:np.array, inner_point:np.array) -> ConvexHull:
    
    #halfspaces = np.concatenate(H,P,axis=1)

    hs = HalfspaceIntersection(HP, inner_point)

    return ConvexHull(np.asarray(list(zip(*hs.intersections))).transpose())

def getBoundary(convex_hull:ConvexHull) -> SetOfPoints:
    return np.array(list(map(lambda i: list(convex_hull.points[i]), convex_hull.vertices)))
    
def matrix_dot_set(matrix: Matrix, convex_hull:ConvexHull) -> ConvexHull:
    new_set:SetOfPoints = list(map(lambda i : np.dot(matrix, convex_hull.points[i]), convex_hull.vertices))
    return ConvexHull(new_set)
    
def Minkowski(convex_hull_a:ConvexHull, convex_hull_b:ConvexHull) -> ConvexHull:

    xablau = list(itertools.product(convex_hull_a.vertices, convex_hull_b.vertices))

    return ConvexHull(np.asarray(list(map(lambda i : list(convex_hull_a.points[i[0]] + convex_hull_b.points[i[1]]), xablau))))
    
def mRPI(Acl:Matrix, w:SetOfPoints, I:SetOfPoints, max_iteration:int = 10) -> ConvexHull:

    EW:SetOfPoints = E(w)
    PHI:Union[SetOfPoints,ConvexHull] = E(I)

    for i in range(max_iteration):

        PHI = matrix_dot_set(Acl, PHI)
        PHI = Minkowski(PHI,EW)
                
        #yield PHI

    return PHI

def whether_inside(ES:ConvexHull, x:np.array) -> boolean:
    '''
    Returns whether or not the point is inside the envelope
    True: It is inside
    False: It is NOT inside
    '''

    hull_path = Path(ES.points)
    
    return hull_path.contains_point(x)


def main():

    ## Inputs ##
    
    # Open Loop matrix 
    Aol = np.array([[1.000, 0.650],
                    [0.000, 1.000]])

    # Dynamic controled matrix A - BK 
    Acl = np.array([[0.879, 0.393],
                    [-0.374, 0.209]])

    # Disturbances ranges (here E is eye(dim(x)))
    w = np.array([[0.211,-0.211],[0.65,-0.65]])

    # First RPI guess set vertices
    I = np.array([[5,-5],[5,-5]])

    # See if the selected point lies inside the calculated mRPI
    test_point = np.array([0.5,1.075])

    ## Code per se ##

    phi = mRPI(Acl, w, I, 50)

    #print(whether_inside(phi.equations, test_point))

    ES = Event_Set(np.array(phi.equations), Aol, w)

    plt.figure()
    for simplex in phi.simplices:
        plt.plot(phi.points[simplex, 0], phi.points[simplex, 1], 'k-')
    plt.grid(True)

    for simplex in ES.simplices:
        plt.plot(ES.points[simplex, 0], ES.points[simplex, 1], 'r--')
    plt.grid(True)

    plt.show()

if __name__ == '__main__':

    main()
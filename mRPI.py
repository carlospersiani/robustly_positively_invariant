##
## São Carlos
## Aug, 2022
##
## Carlos Persiani
## Henrique Garcia
##

import scipy
# reference: PhD thesis from Bruno BOISSEAU entitled "Event-based control: application to robotic systems" (https://www.theses.fr/2017GREAT037.pdf), advisors Nicolas MARCHAND and John-Jairo MARTINEZ-MOLINA

from enum import Enum, auto
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
    # HP is a [A ; B] matrix where Ax<B  is the collection of equations that defines the halfspace expected to be the mRPI for a dynamic system

    amount_of_equations = HP.shape[0]
    amount_of_dimensions = HP.shape[1]-1

    #w_plus = np.array([[w[0,0]],[w[1,0]]])
    w_plus = np.reshape(w[:,0],(amount_of_dimensions,1))
    #w_minus = np.array([[w[0,1]],[w[1,1]]])
    w_minus = np.reshape(w[:,1],(amount_of_dimensions,1))

    H = HP[:,:-1]
    P = np.reshape(HP[:,-1],(amount_of_equations,1))

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

    return h2v_representation(halfspaces, np.array([0.0,0.0]))

# Attempt to use 3x3 system failed: too much RAM needed (around 60,000 equations)
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

def h2v_representation(HP:np.array, inner_point:np.array) -> ConvexHull:
    
    #halfspaces = np.concatenate(H,P,axis=1)
    print(inner_point)
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
    
def mRPI(Acl:Matrix, w:SetOfPoints, I:SetOfPoints, max_iteration:int = 10, plot_RPI=True) -> ConvexHull:

    # set that contains the disturbances
    EW:SetOfPoints = E(w)

    # initial guess for the Robust Positively Invariant (RPI) set
    PHI:Union[SetOfPoints,ConvexHull] = E(I)


    # iteration to improve the RPI staring from its initial guess. as the number of interactions approach inifity, RPI will actually be the exact Minimal Robust Positively Invariant (mPRI) set
    for i in range(max_iteration):
        # next RPI set will be PHI(next) = MinkowskiSum(Acl*PHI(current), EW)

        # performs the Acl*Phi(current)
        PHI = matrix_dot_set(Acl, PHI)
        # apply minkowski sum on EW and the result of the line above
        PHI = Minkowski(PHI,EW)

        #yield PHI

    if plot_RPI:
        fig, ax = plt.subplots(1)
        convex_hull_plot_2d(PHI,ax)
        fig.canvas.draw()
        fig.show()
        plt.pause(.05)


    return PHI

def whether_inside(ES:ConvexHull, x:np.array) -> boolean:
    '''
    Returns whether or not the point is inside the envelope
    True: It is inside
    False: It is NOT inside
    '''

    hull_path = Path(ES.points)
    
    return hull_path.contains_point(x)


class SystemType(Enum):
    TwoByTwo = auto()
    ThreeByThree = auto()


def main():

    ## Inputs ##

    # Open Loop matrix 
    Aol = np.array([[1.000, 0.650],
                    [0.000, 1.000]])

	B=np.array([[.211],[.65]])
	K=np.array([[.575, 1.217]])
    # Dynamic controled matrix A - BK (eig must lie inside the unit circle)
	Acl = Aol-B*K
	#Acl = np.array([[0.879, 0.393],
	#                [-0.374, 0.209]])

    # Disturbances ranges (here E is eye(dim(x)))
    w = np.array([[0.211,-0.211],[0.65,-0.65]])

	# Disturbances ranges (a n-uple for each system dimension)
	w = np.array([[0.211,-0.211],[0.65,-0.65]])

    # See if the selected point lies inside the calculated mRPI
    test_point = np.array([0.5,1.075])

    elif system == SystemType.ThreeByThree:
        # Dynamic controled matrix (eig must lie inside the unit circle)
        Acl = np.array([[0.879, 0.393, .01],
                        [-0.374, 0.209, .01],
                        [0, 0.2, .51]])

        # Disturbances ranges (an n-uple for each system dimension)
        w = np.array([[0.211,-0.211],[0.65,-0.65], [-.1, .1]])

        # First RPI guess set vertices
        I = np.array([[5,-5],[5,-5], [0.1, -.1]])

        # Just to see if the selected point lies inside the calculated mRPI
        test_point = np.array([0.5,1.075, .05])
    else:
        print(f"system should be \"2x2\" or \"3x3\". Currently system is: {system}")

    ## Code per se
    #phi is the set of states as infinity
    phi = mRPI(Acl, w, I, 50, plot_RPI=False)

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
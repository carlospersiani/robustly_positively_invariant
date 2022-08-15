##
## São Carlos
## Aug, 2022
##
## Carlos Persiani
## Henrique Garcia
##

from xmlrpc.client import boolean
from scipy.spatial import ConvexHull, convex_hull_plot_2d, HalfspaceIntersection
import numpy as np
import matplotlib.pyplot as plt
import itertools
from typing import List, Text, Union

SetOfPoints = np.array
Matrix = np.array

def E(w:np.array) -> ConvexHull:

    return ConvexHull(list(itertools.product(*w)))

def Event_Set(HP:np.array, Acl:np.array, w:np.array) -> ConvexHull:
    
    w_plus = np.array([[w[0,0]],[w[1,0]]])
    w_minus = np.array([[w[0,1]],[w[1,1]]])

    H = HP[:,:-1]
    P = np.reshape(HP[:,-1],(H.shape[0],1))

    P_plus = P - np.dot(H, w_plus)
    P_minus = P - np.dot(H, w_minus)

    P_new = np.concatenate((P_plus,P_minus),axis=0)
    H_new = np.concatenate((H,H),axis=0)

    halfspaces = np.concatenate((H_new,P_new), axis=1)

    return h2v_representation(halfspaces, np.array([0.0,0.0]))

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

def whether_inside(matrix_inequality:Matrix, vector:np.array) -> boolean:

    '''
    Returns whether or not the point is inside the envelope
    True: It is inside
    False: It is NOT inside
    '''

    statement = True
    
    for lines in matrix_inequality:

        result = np.inner(lines[:-1],vector) - lines[-1]
        if result < 0: statement = False

    return statement

def main():

    ## Inputs ##
    
    # Dynamic controled matrix (eig must lie inside the unit circle)
    Acl = np.array([[0.879, 0.393],
                    [-0.374, 0.209]])

    # Disturbances ranges (a n-uple for each system dimension)
    w = np.array([[0.211,-0.211],[0.65,-0.65]])

    # First RPI guess set vertices
    I = np.array([[5,-5],[5,-5]])

    # Just to see if the selected point lies inside the calculated mRPI
    test_point = np.array([0.5,1.075])

    ## Code per se

    phi = mRPI(Acl, w, I, 50)

    #print(whether_inside(phi.equations, test_point))

    ES = Event_Set(np.array(phi.equations), Acl, w)

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





##
## São Carlos
## Aug, 2022
##
## Carlos Persiani
## Henrique Garcia
##

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
import matplotlib.pyplot as plt
import itertools
from typing import List, Text, Union

SetOfPoints = np.array
Matrix = np.array

def E(w:np.array) -> ConvexHull:

    return ConvexHull(list(itertools.product(*w)))

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

def main():

    Acl = np.array([[-2, 0],
                    [0, -1]])

    w = np.array([[1,-1],[5,4]])

    points = np.array([[-5.5, 1.1],
                       [-5.2, -0.6],
                       [-6.5, -1],
                       [-7,-2],
                       [-9,12]])

    points_x = np.array([[5, 1],
                         [5.3, -0.5],
                         [6.6, -1.2]])

    points_hull = ConvexHull(points)

    plt.plot(points[:,0], points[:,1], 'o')
    for simplex in points_hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

    points_hull_x = ConvexHull(points_x)

    plt.plot(points_x[:,0], points_x[:,1], 'o')
    for simplex in points_hull_x.simplices:
        plt.plot(points_x[simplex, 0], points_x[simplex, 1], 'k-')


    #print(hull.equations)
    #print(hull.vertices)
    #print(hull.points)

    #ACL_phi = matrix_dot_set(hull, Acl)

    #points2 = getBoundary(ACL_phi)

    E(w)

    mink = Minkowski(points_hull, points_hull_x)

    points2 = getBoundary(mink)

    plt.plot(points2[:,0], points2[:,1], 'o')
    for simplex in mink.simplices:
        plt.plot(mink.points[simplex, 0], mink.points[simplex, 1], 'k-')

    plt.show()

if __name__ == '__main__':
    
    Acl = np.array([[0.879, 0.393],
                    [-0.374, 0.209]])

    w = np.array([[0.211,-0.211],[0.65,-0.65]])

    I = np.array([[5,-5],[5,-5]])

    plt.figure()

    #for phi in mRPI(Acl, w, I,100):
    #    for simplex in phi.simplices:
    #        plt.plot(phi.points[simplex, 0], phi.points[simplex, 1], 'k-')
    
    phi = mRPI(Acl, w, I, 5000)
    for simplex in phi.simplices:
        plt.plot(phi.points[simplex, 0], phi.points[simplex, 1], 'k-')
    plt.grid(True)
    plt.show()


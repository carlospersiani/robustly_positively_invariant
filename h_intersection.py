
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np

def Event_Set(HP:np.array, Acl:np.array) -> ConvexHull:
    
    H = HP[:,:-1]
    P = np.reshape(HP[:,-1],(H.shape[0],1))

    halfspaces = np.concatenate((np.dot(H, Acl),P), axis=1)

    return h2v_representation(halfspaces, np.array([0,0]))

def h2v_representation(HP:np.array, inner_point:np.array) -> ConvexHull:

    #halfspaces = np.concatenate(H,P,axis=1)

    hs = HalfspaceIntersection(HP, inner_point)

    return ConvexHull(np.asarray(list(zip(*hs.intersections))).transpose())

halfspaces = np.array([[-1, 0., 0.],
                       [0., -1., 0.],
                       [2., 1., -4.],
                       [-0.5, 1., -2.]])

feasible_point = np.array([0.01, 0.01])

Acl = np.array([[0.879, 0.393],
                [-0.374, 0.209]])

ES = Event_Set(halfspaces, Acl)



'''

phi = h2v_representation(halfspaces, feasible_point)

plt.figure()
for simplex in phi.simplices:
    plt.plot(phi.points[simplex, 0], phi.points[simplex, 1], 'k-')
plt.grid(True)

plt.show()
'''
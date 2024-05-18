'''
    Calculation of energy-distances from points to hull
'''

import numpy as np
import scipy.spatial
import matplotlib.path as mpltPath

np.set_printoptions(precision=3)

def dist(hull, debug=False):

    NP = hull.points.shape[0]  # total points
    NT = hull.simplices.shape[0] # no of triangles
    distances = np.zeros(NP)  # distances array

    for P in range(NP):
        D = np.inf
        for T in range(NT):
            trio = hull.simplices[T]
            eq = hull.equations[T]
            Dn = (eq[:-1] @ hull.points[P] + eq[-1]) / eq[2]
            if Dn < D and Dn > 0:
                D = Dn
        distances[P] = D

    return distances


def dist_old(hull, TOL=1e-6, radius=1e-6, distance_messages=False, debug=False):

    points = np.copy(hull.points)
    #####################################################################
    # HACK: bringing copy of binary points to within the triangle
    mask = np.where((hull.points[:, 0] < TOL) & (hull.points[:, 1] > TOL))
    points[mask, 0] = radius
    mask = np.where((hull.points[:, 0] > TOL) & (hull.points[:, 1] < TOL))
    points[mask, 1] = radius
    mask = np.where(points[:, 0] + hull.points[:, 1] > 1 - TOL)
    points[mask, 0] = points[mask, 0] * (1  - radius)
    points[mask, 1] = points[mask, 1] * (1  - radius)
    ####################################################################

    NP = hull.points.shape[0]  # total points
    NT = hull.simplices.shape[0] # no of triangles
    distances = np.zeros(NP)  # distances array
    
    for P in range(NP):  # all points
        if P not in hull.vertices:  # only for points not in hull
            for T in range(NT):  # each triangle (simplex)
                trio = hull.simplices[T]
                eq = hull.equations[T]
                if np.abs(eq[-1]) < TOL and np.abs(eq[0]) < TOL and np.abs(eq[1]) < TOL and np.abs(eq[2]) > 1 - TOL:
                    message = f'Point #{P}={hull.points[P, :2]}: simplex {trio} with eq={eq} is the lid of the hull; distance not calculated'
                else:
                    path = mpltPath.Path(hull.points[trio, :2])
                    if path.contains_point(points[P, :2], radius=radius):  # HACK: use points here, not hull.points!
                        message = f'Point #{P}={hull.points[P, :2]} in simplex {trio}, with eq: {eq}'
                        if np.abs(eq[2]) < TOL:
                            message += ', which is horizontal; distance not calculated'
                        else:
                            D = (eq[:-1] @ hull.points[P] + eq[-1]) / eq[2]
                            distances[P] = D
                            message += f'; distance: {D}'
                        if distance_messages:
                            print(message)
        else:
            if distance_messages:
                print(f'Point #{P}={hull.points[P, :2]} already in hull; distance calculation unnecessary')
        
    if debug:
        nprob = 0
        print('PROBLEMATIC POINTS:')
        for P in range(len(hull.points)):
            if np.abs(distances[P]) < TOL and P not in hull.vertices:
                print(f'#{P}: {hull.points[P, :2]} {distances[P]}')
                nprob += 1
                for T in range(NT):  # each triangle (simplex)
                    trio = hull.simplices[T]
                    eq = hull.equations[T]

                    path = mpltPath.Path(hull.points[trio, :2])
                    if path.contains_point(points[P, :2], radius=radius):  # HACK: use points here, not hull.points!
                        print(f'\ttrio: {T}: eq={eq}')
        print(f'# TOTAL: {nprob}')
    return distances
            
if __name__ == '__main__':
    
    points = np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [.2, .3, -1], [.25, .4, -.5]])
          
    hull = scipy.spatial.ConvexHull(points)
    print(hull.vertices)
    print(hull.simplices)
    print(hull.equations)

    print(dist(hull))

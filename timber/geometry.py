# timber

import math
import numpy as np

##############################################################################

## x,y,z cartesian coordinate object ##
class Coord(object):
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z

## distance ##
def cart_distance(coord1,coord2):
    x1=coord2.x-coord1.x
    y1=coord2.y-coord1.y
    z1=coord2.z-coord1.z

    tot=(x1*x1)+(y1*y1)+(z1*z1)
    return math.sqrt(tot)

## angle ##
def get_angle(point1,point2,point3):

    a = np.array([point1.x,point1.y,point1.z])
    b = np.array([point2.x,point2.y,point2.z])
    c = np.array([point3.x,point3.y,point3.z])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

## dihedral ##
def get_dihedral(point1,point2,point3,point4):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array([point1.x,point1.y,point1.z])
    p1 = np.array([point2.x,point2.y,point2.z])
    p2 = np.array([point3.x,point3.y,point3.z])
    p3 = np.array([point4.x,point4.y,point4.z])

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


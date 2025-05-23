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

## centroid of points ##
def find_centroid(coord_list):
    assert isinstance(coord_list,list)

    x_tot=0.0
    y_tot=0.0
    z_tot=0.0

    for point in coord_list:
        x_tot=x_tot+point.x
        y_tot=y_tot+point.y
        z_tot=z_tot+point.z

    return Coord(x_tot/float(len(coord_list)),y_tot/float(len(coord_list)),z_tot/float(len(coord_list)))

## closest to anchor ##
def rank_closest(anchor,point_list):

    dist_val=[]
    for point in point_list:
        dist_val.append(cart_distance(anchor,point))

    return int(dist_val.index(min(dist_val)))

## further from anchor ##
def rank_further(anchor,point_list):

    dist_val=[]
    for point in point_list:
        dist_val.append(cart_distance(anchor,point))

    return int(dist_val.index(max(dist_val)))



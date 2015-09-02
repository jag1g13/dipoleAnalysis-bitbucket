
import numpy as np
from math import sqrt
from math import acos
from math import pi
import sys

def normalToPlane(a, b, c):
    """
    Calculate the vector normal to a plane containing 3 atoms, by taking the cross product of the two vectors
    forming the "triangle". e.g. for atoms a, b & c, using the cross product: (a->b)x(b->c).
    :param a: Atom coordinates as numpy array
    :param b: Atom coordinates as numpy array
    :param c: Atom coordinates as numpy array
    :return: Vector normal to plane containing 3 atoms
    """
    x1 = b[0] - a[0]
    y1 = b[1] - a[1]
    z1 = b[2] - a[2]
    bond1 = [x1, y1, z1]
    # bond1 = b - a

    x2 = c[0] - b[0]
    y2 = c[1] - b[1]
    z2 = c[2] - b[2]
    bond2 = [x2, y2, z2]
    # bond2 = c - b

    norm = np.cross(bond1, bond2)
    norm_mag = sqrt((norm[0]*norm[0])+(norm[1]*norm[1])+(norm[2]*norm[2]))
    norm = [(norm[0]/norm_mag), (norm[1]/norm_mag), (norm[2]/norm_mag)]
    # norm /= norm_mag

    return norm

def angleBetweenVectors(a, b):
    """
    Calculate the angle between two vectors using the following formula:

                                           a . b
                                cos(x) = ---------
                                          |a| |b|

    For now, this function only returns an unsigned angle...
    :param a: Vector
    :param b: Vector
    :return: Angle between the two vectors (in degrees)
    """
    a_mag = sqrt((a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]))
    b_mag = sqrt((b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]))
    # a_mag = np.sqrt(np.dot(a, a))
    # b_mag = np.sqrt(np.dot(b, b))

    dotProd = (a[0]*b[0]) + (a[1]*b[1]) + (a[2]*b[2])
    # dotProd = np.dot(a, b)

    cosX = dotProd / (a_mag*b_mag)
    X_rad = acos(cosX)
    # X_deg = X_rad*180 / pi

    # return X_deg
    return X_rad

if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)

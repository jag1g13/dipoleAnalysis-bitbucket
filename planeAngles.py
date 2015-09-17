
import numpy as np
from math import sqrt
from math import cos
from math import sin
from math import pi
from math import atan2
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
    Calculate the signed angle between two vectors.
    :param a: Vector
    :param b: Vector
    :return: Angle between the two vectors
    """
    a_mag = sqrt((a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]))
    b_mag = sqrt((b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]))
    a = [(a[0]/a_mag), (a[1]/a_mag), (a[2]/a_mag)]
    b = [(b[0]/b_mag), (b[1]/b_mag), (b[2]/b_mag)]
    # a_mag = np.sqrt(np.dot(a, a))
    # b_mag = np.sqrt(np.dot(b, b))

    dotProd = (a[0]*b[0]) + (a[1]*b[1]) + (a[2]*b[2])
    crossProd = np.cross(a, b)
    crossMag = sqrt(crossProd[0]*crossProd[0] + crossProd[1]*crossProd[1] + crossProd[2]*crossProd[2])

    # dotProd = np.dot(a, b)

    # cosX = dotProd
    # X_rad = acos(cosX)
    X_rad = atan2(crossMag, dotProd)
    # X_deg = X_rad*180 / pi

    # return X_deg
    return X_rad


def coneCenter(atom1, atom2, atom3, improp, tolerance):
    """
    Function to find the center of the cone around which dipoles rotate
    :param atom1: Atom (containing dipole)
    :param atom2: Atom
    :param atom3: Atom
    :param improp: Improper dihedral angle (using the cone center, not the dipole)
    :param tolerance: Tolerance for the improper match
    :return: Vector at the cone center
    """
    a_vec = np.zeros(3)
    a_vec[0] = atom1.coords[0]-atom2.coords[0]
    a_vec[1] = atom1.coords[1]-atom2.coords[1]
    a_vec[2] = atom1.coords[2]-atom2.coords[2]
    b_vec = np.zeros(3)
    b_vec[0] = atom1.coords[0]-atom3.coords[0]
    b_vec[1] = atom1.coords[1]-atom3.coords[1]
    b_vec[2] = atom1.coords[2]-atom3.coords[2]

    n_vec = np.cross(a_vec, b_vec)
    n_mag = sqrt((n_vec[0]*n_vec[0]) + (n_vec[1]*n_vec[1]) + (n_vec[2]*n_vec[2]))
    n_vec /= n_mag

    c_vec = a_vec + b_vec
    c_mag = sqrt((c_vec[0]*c_vec[0]) + (c_vec[1]*c_vec[1]) + (c_vec[2]*c_vec[2]))
    c_vec /= c_mag

    u_vec = np.cross(n_vec, c_vec)
    u_mag = sqrt((u_vec[0]*u_vec[0]) + (u_vec[1]*u_vec[1]) + (u_vec[2]*u_vec[2]))
    u_vec /= u_mag

    theta = (54*pi/180) # Start testing half of 109.5 deg
    improp *= (pi/180)
    tolerance *= (pi/180)

    # Defining elements of rotation matrix
    R_00 = cos(theta) + (u_vec[0]*u_vec[0])*(1 - cos(theta))
    R_01 = (u_vec[0]*u_vec[1])*(1 - cos(theta)) - u_vec[2]*sin(theta)
    R_02 = (u_vec[0]*u_vec[2])*(1 - cos(theta)) + u_vec[1]*sin(theta)
    R_10 = (u_vec[1]*u_vec[0])*(1 - cos(theta)) + u_vec[2]*sin(theta)
    R_11 = cos(theta) + (u_vec[1]*u_vec[1])*(1 - cos(theta))
    R_12 = (u_vec[1]*u_vec[2])*(1 - cos(theta)) - u_vec[0]*sin(theta)
    R_20 = (u_vec[2]*u_vec[0])*(1 - cos(theta)) - u_vec[1]*sin(theta)
    R_21 = (u_vec[2]*u_vec[1])*(1 - cos(theta)) + u_vec[0]*sin(theta)
    R_22 = cos(theta) + (u_vec[2]*u_vec[2])*(1 - cos(theta))

    # print(" / " + str(R_00) + " " + str(R_01) + " " + str(R_02) + " \\")
    # print("| " + str(R_10) + " " + str(R_11) + " " + str(R_12) + "  |")
    # print(" \\ " + str(R_20) + " " + str(R_21) + " " + str(R_22) + " /")

    d_vec = np.zeros(3)
    d_vec[0] = R_00*c_vec[0] + R_01*c_vec[1] + R_02*c_vec[2]
    d_vec[1] = R_10*c_vec[0] + R_11*c_vec[1] + R_12*c_vec[2]
    d_vec[2] = R_20*c_vec[0] + R_21*c_vec[1] + R_22*c_vec[2]

    if abs((float(improp) - improperVector(d_vec, atom1, atom2, atom3))) < float(tolerance): # Need to create or adapt a function for this improper
        print(c_vec)
        print(d_vec)
        print(" ")
        return d_vec  # Successful find of vector
    else:
        print("\nIncorrect first guess. Difference of:")
        print(str(abs(float(improp) - improperVector(d_vec, atom1, atom2, atom3))*(180/pi)) + " degrees")

        theta *= (-1)

        R_00 = cos(theta) + (u_vec[0]*u_vec[0])*(1 - cos(theta))
        R_01 = (u_vec[0]*u_vec[1])*(1 - cos(theta)) - u_vec[2]*sin(theta)
        R_02 = (u_vec[0]*u_vec[2])*(1 - cos(theta)) + u_vec[1]*sin(theta)
        R_10 = (u_vec[1]*u_vec[0])*(1 - cos(theta)) + u_vec[2]*sin(theta)
        R_11 = cos(theta) + (u_vec[1]*u_vec[1])*(1 - cos(theta))
        R_12 = (u_vec[1]*u_vec[2])*(1 - cos(theta)) - u_vec[0]*sin(theta)
        R_20 = (u_vec[2]*u_vec[0])*(1 - cos(theta)) - u_vec[1]*sin(theta)
        R_21 = (u_vec[2]*u_vec[1])*(1 - cos(theta)) + u_vec[0]*sin(theta)
        R_22 = cos(theta) + (u_vec[2]*u_vec[2])*(1 - cos(theta))

        d_vec = np.zeros(3)
        d_vec[0] = R_00*c_vec[0] + R_01*c_vec[1] + R_02*c_vec[2]
        d_vec[1] = R_10*c_vec[0] + R_11*c_vec[1] + R_12*c_vec[2]
        d_vec[2] = R_20*c_vec[0] + R_21*c_vec[1] + R_22*c_vec[2]

    if abs((float(improp) - improperVector(d_vec, atom1, atom2, atom3))) < float(tolerance):
        print(c_vec)
        print(d_vec)
        print(" ")
        return d_vec # Have now found the cone center
    else:
        print("Could not find suitable vector of center of dipole cone. Difference of:")
        print(str(abs(float(improp) - improperVector(d_vec, atom1, atom2, atom3))*(180/pi)) + " degrees")
        return [0, 0, 0]
        # sys.exit()


def improperVector(vector, atom1, atom2, atom3):
    """
    Function to test the improper dihedral angle of the center of a dipole cone.
    :param vector: Vector when not using dipoles
    :param atom1: Atom
    :param atom2: Atom
    :param atom3: Atom
    :return: improper angle
    """

    crossProd1 = np.cross(vector, (atom1.coords-atom2.coords))  # consistency with dihedral_dipole.cpp
    mag1 = sqrt((crossProd1[0]*crossProd1[0])+(crossProd1[1]*crossProd1[1])+(crossProd1[2]*crossProd1[2]))
    crossProd1 /= mag1
    crossProd2 = np.cross((atom3.coords-atom2.coords), (atom1.coords-atom2.coords))
    mag2 = sqrt((crossProd2[0]*crossProd2[0])+(crossProd2[1]*crossProd2[1])+(crossProd2[2]*crossProd2[2]))
    crossProd2 /= mag2

    middle = atom2.coords-atom1.coords
    mag3 = sqrt(middle[0]*middle[0] + middle[1]*middle[1] + middle[2]*middle[2])
    middle /= mag3

    dot = np.dot(crossProd1, crossProd2)
    det = (crossProd1[0]*crossProd2[1]*middle[2] - crossProd1[0]*crossProd2[2]*middle[1] - crossProd1[1]*crossProd2[0]*middle[2])
    det += (crossProd1[1]*crossProd2[2]*middle[0] + crossProd1[2]*crossProd2[0]*middle[1] - crossProd1[2]*crossProd2[1]*middle[0])
    torsion = (-1)*atan2(det, dot)

    return torsion


if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)



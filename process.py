#!/usr/bin/env python

import sys
import numpy as np
import time
import matplotlib.pyplot as plt
from math import sqrt
from optparse import OptionParser

from FrameReader import FrameReader
import planeAngles


class Atom:
    def __init__(self, atom_type=None, coords=None, dipole=None):
        self.atom_id = atom_type

        if coords is None:
            self.coords = np.zeros(3)
        else:
            self.coords = coords

        if dipole is None:
            self.dipole = np.zeros(3)
        else:
            self.dipole = dipole

    def __repr__(self):
        return "<Atom {0} @ {1:8.3f}, {2:8.3f}, {3:8.3f} with dipole {4:8.3f}, {5:8.3f}, {6:8.3f}>"\
               .format(self.atom_id, self.coords[0], self.coords[1], self.coords[2],
                       self.dipole[0], self.dipole[1], self.dipole[2])


class Frame:
    def __init__(self, natoms):
        self.natoms = natoms
        self.atoms = []
        for i in xrange(natoms):
            self.atoms.append(Atom())

    def __repr__(self):
        return "<Frame containing {1} atoms>".format(natoms)

    def angle_norm_bisect(self, num1, num2, num3):
        """
        return normal vector to plane formed by 3 atoms and their bisecting vector
        """
        vec1 = (self.atoms[num2].loc - self.atoms[num1].loc)
        vec2 = (self.atoms[num3].loc - self.atoms[num2].loc)
        vec1 = vec1 / np.linalg.norm(vec1)
        vec2 = vec2 / np.linalg.norm(vec2)
        normal = np.cross(vec1, vec2)
        bisec = (vec1 + vec2) / 2.
        return polar_coords(normal), polar_coords(bisec)

    def show_atoms(self, start=0, end=-1):
        """
        print coordinates of the atoms numbered start to end
        """
        if end == -1:
            end = len(self.atoms)
        for i in xrange(start, end):
            print(self.atoms[i])

    def dipoleAngle(self, a, b):
        """
        Calculate angle of dipole on atom a with respect to bond vector a->b
        :param a: Atom number
        :param b: Atom number
        :return: Angle of dipole with respect to bond
        """
        if np.any(self.atoms[a].dipole):
            return planeAngles.angleBetweenVectors(self.atoms[a].dipole,
                                                   self.atoms[b].coords-self.atoms[a].coords)
        else:
            return 0

    def planeNormal(self, a, b, c):
        """
        Calculate normal to plane containing three atoms
        :param a: Atom number
        :param b: Atom number
        :param c: Atom number
        :return: Normal to plane containing all three atoms
        """
        return planeAngles.normalToPlane(self.atoms[a].coords,
                                         self.atoms[b].coords,
                                         self.atoms[c].coords)



def polar_coords(xyz, axis1=np.array([0, 0, 0]), axis2=np.array([0, 0, 0]), mod=True):
    """
    Convert cartesian coordinates to polar, if axes are given it will be reoriented.
    axis points to the north pole (latitude), axis2 points to 0 on equator (longitude)
    if mod, do angles properly within -pi, +pi
    """
    tpi = 2*np.pi
    polar = np.zeros(3)
    xy = xyz[0]**2 + xyz[1]**2
    polar[0] = np.sqrt(xy + xyz[2]**2)
    polar[1] = np.arctan2(np.sqrt(xy), xyz[2]) - axis1[1]
    polar[2] = np.arctan2(xyz[1], xyz[0]) - axis2[2]

    if axis2[1] < 0:
        polar[2] = polar[2] + tpi
    if mod:
        polar[1] = polar[1] % (tpi)
        polar[2] = polar[2] % (tpi)

    return polar


def print_output(output_all, output, request):
    for name, val in zip(request, output):
        print("{0}: {1:4.3f}".format("-".join(name), val))


def graph_output(output_all):
    rearrange = zip(*output_all)
    plt.figure()

    for i, item in enumerate(rearrange):
        plt.subplot(2, 3, i+1)
        data = plt.hist(item, bins=100, normed=1)


def boltzmannInvert(vals, temp=300):
    """
    Perform Boltzmann Inversion on a list of numbers assumed to be normally distributed
    :param vals: Array-like of numbers
    :param vals: Temperature of simulation in Kelvin, default 300
    :return: Tuple containing mean value and force constant
    """
    if not np.any(vals):
        return (0, 0)

    mean = 180 * np.mean(vals) / np.pi
    sdev = np.std(vals)
    fc = 1.987e-3 * temp / (sdev*sdev)

    return (mean, fc)


def calcAnglesAll(frame, offset=1, natoms=6):
    """
    Calculate dipole angle for every atom in ring
    :param frame: Frame instance
    :param offset: Calculate angle with respect to N around the ring
    :param natoms: Number of atoms in ring
    :return: Numpy array containing angles
    """
    angles = np.zeros(natoms)
    for i in xrange(natoms):
        angles[i] = frame.dipoleAngle(i, (i+offset) % natoms)
    return angles


def calcAnglesPlane(frame, offset=2, natoms=6):
    """
    Calculate angle between dipole and plane for every atom in ring
    :param frame: Frame instance
    :param offset: Calculate angle with respect to N around the ring
    :param natoms: Number of atoms in ring
    :return: Numpy array containing angles
    """
    angles = np.zeros(natoms)
    for i in xrange(natoms):
        if np.any(frame.atoms[i].dipole):
            norm = frame.planeNormal(i, (i+offset)%natoms, (i+2*offset)%natoms)
            angles[i] = planeAngles.angleBetweenVectors(frame.atoms[i].dipole, norm)
    return angles


def analyse(filename, natoms=-1):
    """
    Perform analysis of dipoles in LAMMPS trajectory
    :param lammpstrj: Filename of LAMMPS trajectory
    :return: Number of frames in trajectory
    """
    np.set_printoptions(precision=3, suppress=True)
    reader = FrameReader(filename)

    if natoms == -1:
        natoms = reader.total_atoms
    else:
        if natoms > reader.total_atoms:
            print("ERROR: Requested more atoms than are in trajectory")
            sys.exit(1)
    nframes = reader.total_frames

    angle1 = np.zeros((nframes, 6))
    angle2 = np.zeros((nframes, 6))
    angle3 = np.zeros((nframes, 6))

    frame = Frame(natoms)
    for i in xrange(nframes):
        # Read in frame from trajectory and process
        reader.readFrame(i, frame)
        angle1_tmp = calcAnglesAll(frame)
        angle2_tmp = calcAnglesAll(frame, 3)
        angle3_tmp = calcAnglesPlane(frame)
        # frame.show_atoms(0,6)

        for j in xrange(6):
            angle1[i, j] = angle1_tmp[j]
            angle2[i, j] = angle2_tmp[j]
            angle3[i, j] = angle3_tmp[j]

    for j in xrange(6):
        print("-"*5)
        # plt.hist(angle1[:, j], 50, normed=1)
        # plt.show()
        print(boltzmannInvert(angle1[:, j]))
        print(boltzmannInvert(angle2[:, j]))
        print(boltzmannInvert(angle3[:, j]))

    # analyseAngles(angle1)
    # analyseAngles(angle2)

    return nframes


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input",
                      action="store", type="string", dest="lammpstrj", default="",
                      help="Input file - LAMMPS trajectory", metavar="FILE")
    parser.add_option("-n", "--natoms",
                      action="store", type="int", dest="natoms", default="-1",
                      help="Number of atoms to calculate for")
    # parser.add_option("-v", "--verbose",
    #                   action="store_true", dest="verbose", default=False,
    #                   help="Make more verbose")
    (options, args) = parser.parse_args()
    print("="*25)
    if not options.lammpstrj:
        print("Must provide LAMMPS trajectory to run")
        sys.exit(1)
    t_start = time.clock()
    nframes = analyse(options.lammpstrj)
    t_end = time.clock()
    print("-"*25 + "\nCalculated {0} frames in {1}s\n".format(nframes, (t_end - t_start)) + "="*25)

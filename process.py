#!/usr/bin/env python

import sys
import numpy as np
import time
import matplotlib.pyplot as plt
from optparse import OptionParser

from FrameReader import FrameReader
from Frame import Frame
import planeAngles


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

    mdat = np.ma.masked_array(vals, np.logical_or(np.isnan(vals), np.equal(vals, np.zeros_like(vals))))
    mean = 180 * np.mean(mdat) / np.pi
    sdev = np.std(mdat)
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
        if np.any(frame.atoms[i].dipole):
            angles[i] = frame.dipoleAngle(i, (i+offset) % natoms)
    return angles


def calcImpropersAll(frame, natoms=6):
    """
    Calculate dipole improper angle for every atom in ring
    :param frame: Frame instance
    :param natoms: Number of atoms in ring
    :return: Numpy array containing impropers
    """
    impropers = np.zeros(natoms)
    for i in xrange(natoms):
        if np.any(frame.atoms[i].dipole):
            impropers[i] = frame.dipoleImproper(i, (i+1) % natoms, (i+5) % natoms)
    return impropers


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
            norm = frame.planeNormal(i, (i+offset) % natoms, (i+2*offset) % natoms)
            angles[i] = planeAngles.angleBetweenVectors(frame.atoms[i].dipole, norm)
    return angles


def main(filename, nframes=-1, natoms=-1):
    """
    Perform analysis of dipoles in LAMMPS trajectory
    :param lammpstrj: Filename of LAMMPS trajectory
    :return: Number of frames in trajectory
    """
    np.set_printoptions(precision=3, suppress=True)
    reader = FrameReader(filename)

    if natoms < 0:
        natoms = reader.total_atoms
    else:
        natoms = min(natoms, reader.total_atoms)

    if nframes < 0:
        nframes = reader.total_frames
    else:
        nframes = min(nframes, reader.total_frames)

    angle1 = np.zeros((nframes, 6))
    angle2 = np.zeros((nframes, 6))
    improper = np.zeros((nframes, 6))
    angle3 = np.zeros((nframes, 6))

    frame = Frame(natoms)
    print(nframes)
    for i in xrange(nframes):
        # Read in frame from trajectory and process
        if i % 1000 == 0:
            print(i)
        reader.readFrame(i, frame)
        angle1_tmp = calcAnglesAll(frame)
        angle2_tmp = calcAnglesAll(frame, 3)
        angle3_tmp = calcAnglesPlane(frame)
        improper_tmp = calcImpropersAll(frame)
        # frame.show_atoms(0,6)

        if not np.any(angle1_tmp):
            print("AAAAAGH!!!!")

        for j in xrange(6):
            angle1[i, j] = angle1_tmp[j]
            angle2[i, j] = angle2_tmp[j]
            improper[i, j] = improper_tmp[j]
            angle3[i, j] = angle3_tmp[j]

    for j in xrange(6):
        print("-"*5)
        print(boltzmannInvert(angle1[:, j]))
        print(boltzmannInvert(angle2[:, j]))
        print(boltzmannInvert(angle3[:, j]))
        print(boltzmannInvert(improper[:, j]))
        plotHistogram(reduceArrays(improper[:, j]))

    # analyseAngles(angle1)
    # analyseAngles(angle2)

    return nframes


def plotHistogram(array):
    plt.hist(array, 100, normed=1)
    plt.show()


def reduceArrays(array):
    m_array = np.ma.masked_array(array, np.logical_or(np.isnan(array), np.equal(array, np.zeros_like(array))))
    return m_array.compressed()


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input",
                      action="store", type="string", dest="lammpstrj", default="",
                      help="Input file - LAMMPS trajectory", metavar="FILE")
    parser.add_option("-n", "--natoms",
                      action="store", type="int", dest="natoms", default="-1",
                      help="Number of atoms to calculate for")
    parser.add_option("-f", "--nframes",
                      action="store", type="int", dest="nframes", default="-1",
                      help="Number of frames to calculate")
    # parser.add_option("-v", "--verbose",
    #                   action="store_true", dest="verbose", default=False,
    #                   help="Make more verbose")
    (options, args) = parser.parse_args()
    print("="*25)
    if not options.lammpstrj:
        print("Must provide LAMMPS trajectory to run")
        sys.exit(1)
    t_start = time.clock()
    nframes = main(options.lammpstrj, options.nframes)
    t_end = time.clock()
    print("-"*25 + "\nCalculated {0} frames in {1}s\n".format(nframes, (t_end - t_start)) + "="*25)

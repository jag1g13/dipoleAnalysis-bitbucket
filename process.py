#!/usr/bin/env python

import sys
import numpy as np
import time
import matplotlib.pyplot as plt
from math import sqrt
from optparse import OptionParser

from FrameReader import FrameReader

class Atom:
    def __init__(self, atom_type=None, coords=None, dipole=None):
        self.atom_type = atom_type

        if coords is None:
            self.coords = np.zeros(3)
        else:
            self.coords = coords

        if dipole is None:
            self.coords = np.zeros(3)
        else:
            self.dipole = dipole

    def __repr__(self):
        return "<Atom {0} @ {1}, {2}, {2} with dipole {4}, {5}, {6}>"\
               .format(self.atom_type, self.coords[0], self.coords[1], self.coords[2],
                       self.dipole[0], self.dipole[1], self.dipole[2])

class Frame:
    def __init__(self, natoms):
        self.natoms = natoms
        self.atoms = []
        for i in xrange(natoms):
            self.atoms.append(Atom())

    def __repr__(self):
        return "<Frame {0} containing {1} atoms>\n{2}".format(self.num, len(self.atoms), self.title)
        
    def bond_length(self, num1, num2):
        """
        return Euclidean distance between two atoms
        """
        dist = np.linalg.norm(self.atoms[num1].loc - self.atoms[num2].loc)
        #print("Distance between atoms {0} and {1} is {2}".format(num1, num2, dist))
        return dist

    def bond_length_atoms(self, atom1, atom2):
        """
        return Euclidean distance between two atoms
        """
        dist = sqrt((atom1.loc[0]-atom2.loc[0])**2 +  (atom1.loc[1]-atom2.loc[1])**2 + (atom1.loc[2]-atom2.loc[2])**2)
        #dist = sqrt(sum((atom1.loc-atom2.loc)**2))
        #dist = np.linalg.norm(atom1.loc - atom2.loc)
        #print("Distance between atoms {0} and {1} is {2}".format(num1, num2, dist))
        return dist
    
    def bond_angle(self, num1, num2, num3, num4):
        """
        angle at atom2 formed by bonds: 1-2 and 2-3
        """
        vec1 = (self.atoms[num2].loc - self.atoms[num1].loc)
        vec2 = (self.atoms[num4].loc - self.atoms[num3].loc)
        #vec1 = vec1 / np.linalg.norm(vec1)
        vec1 = vec1 / sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2)
        #vec2 = vec2 / np.linalg.norm(vec2)
        vec2 = vec2 / sqrt(vec2[0]**2 + vec2[1]**2 + vec2[2]**2)
        angle = np.arccos(np.dot(vec1, vec2))
        return 180 - (angle * 180 / np.pi)
    
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
        print(self.title)
        if end == -1:
            end = len(self.atoms)
        for i in xrange(start, end):
            print(self.atoms[i])

def calc_measures(frames, req, request, export=True):
    print("Calculating bond "+req+"s")
    if export:
        f = open("bond_"+req+"s.csv", "a")
    t_start = time.clock()
    measures = []
    if export:
        measure_names = ["-".join(name) for name in request]
        f.write(",".join(measure_names) + "\n")
    for i, frame in enumerate(frames):
        perc = i * 100. / len(frames)
        if(i%100 == 0):
            sys.stdout.write("\r{:2.0f}% ".format(perc) + "X" * int(0.2*perc) + "-" * int(0.2*(100-perc)) )
            sys.stdout.flush()
        measures.append(frame.calc_measure[req](request))
        if export:
            measures_text = [str(num) for num in measures[-1]]
            f.write(",".join(measures_text) + "\n")
    avg = np.mean(measures, axis=0)
    t_end = time.clock()
    if export:
        f.truncate()
        f.close()
    print("\rCalculated {0} frames in {1}s\n".format(len(frames), (t_end - t_start)) + "-"*20)
    return measures, avg


def polar_coords(xyz, axis1=np.array([0,0,0]), axis2=np.array([0,0,0]), mod=True):
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

    if axis2[1]<0:
        polar[2] = polar[2] + tpi
    if mod:
        polar[1] = polar[1]%(tpi)
        polar[2] = polar[2]%(tpi)

    return polar


def print_output(output_all, output, request):
    for name, val in zip(request, output):
        print("{0}: {1:4.3f}".format("-".join(name), val))


def graph_output(output_all):
    rearrange = zip(*output_all)
    plt.figure()

    for i, item in enumerate(rearrange):
        plt.subplot(2,3, i+1)
        data = plt.hist(item, bins=100, normed=1)


def analyse(filename, natoms):
    """
    Perform analysis of dipoles in LAMMPS trajectory
    :param lammpstrj: Filename of LAMMPS trajectory
    :return: Nothing
    """
    t_start = time.clock()
    np.set_printoptions(precision=3, suppress=True)
    reader = FrameReader(filename)

    if natoms == -1:
        natoms = reader.total_atoms
    else:
        if natoms > reader.total_atoms:
            print("ERROR: Requested more atoms than are in trajectory")
            sys.exit(1)
    nframes = reader.total_frames

    angle1 = np.zeros((nframes, natoms))
    angle2 = np.zeros((nframes, natoms))

    frame = Frame(natoms)
    for i in xrange(nframes):
        # Read in frame from trajectory and process
        reader.getFrame(frame)
        angle1_tmp = calcAngles1(frame)
        angle2_tmp = calcAngles2(frame)

        for j in xrange(3):
            angle1[i, j] = angle1_tmp[j]
            angle2[i, j] = angle2_tmp[j]

    analyseAngles(angle1)
    analyseAngles(angle2)

    t_end = time.clock()
    print("\rCalculated {0} frames in {1}s\n".format(len(cg_frames), (t_end - t_start)) + "-"*20)
    return len(cg_frames)


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
    if not options.lammpstrj:
        print("Must provide LAMMPS trajectory to run")
        sys.exit(1)
    analyse(options.lammpstrj, options.natoms)

import numpy as np
import planeAngles
import math
import sys


class Atom:
    def __init__(self, atom_type=None, coords=None, dipole=None):
        self.id = 0
        self.type = 0
        self.mol = 0

        if coords is None:
            self.coords = np.zeros(3)
        else:
            self.coords = coords

        if dipole is None:
            self.dipole = np.zeros(3)
        else:
            self.dipole = dipole

        self.mass = 0
        self.diameter = 0

    def __repr__(self):
        return "<Atom {0} @ {1:8.3f}, {2:8.3f}, {3:8.3f} with dipole {4:8.3f}, {5:8.3f}, {6:8.3f}>" \
            .format(self.atom_id, self.coords[0], self.coords[1], self.coords[2],
                    self.dipole[0], self.dipole[1], self.dipole[2])


class Frame:
    def __init__(self, natoms):
        self.natoms = natoms
        self.box = np.zeros((3, 2))
        self.timestep = 0
        self.atomformat = ""
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
        if not np.any(self.atoms[a].dipole):
            return 0

        return planeAngles.angleBetweenVectors(self.atoms[a].dipole,
                                               self.atoms[b].coords-self.atoms[a].coords)

    def dipoleImproper(self, a, b, c):
        """
        Calculate improper of dipole on atom a with respect to bond vector b->c (or c->b?)
        :param a: Atom number (with dipole)
        :param b: Atom number
        :param c: Atom number
        :return: Improper dihedral angle
        """
        if not np.any(self.atoms[a].dipole):
            return 0

        crossProd1 = np.cross(self.atoms[a].dipole, (self.atoms[a].coords-self.atoms[b].coords))  # consistency with dihedral_dipole.cpp
        mag1 = math.sqrt((crossProd1[0]*crossProd1[0])+(crossProd1[1]*crossProd1[1])+(crossProd1[2]*crossProd1[2]))
        crossProd1 /= mag1
        crossProd2 = np.cross((self.atoms[c].coords-self.atoms[b].coords), (self.atoms[a].coords-self.atoms[b].coords))
        mag2 = math.sqrt((crossProd2[0]*crossProd2[0])+(crossProd2[1]*crossProd2[1])+(crossProd2[2]*crossProd2[2]))
        crossProd2 /= mag2

        middle = self.atoms[b].coords-self.atoms[a].coords
        mag3 = math.sqrt(middle[0]*middle[0] + middle[1]*middle[1] + middle[2]*middle[2])
        middle /= mag3

        dot = np.dot(crossProd1, crossProd2)
        det = (crossProd1[0]*crossProd2[1]*middle[2] - crossProd1[0]*crossProd2[2]*middle[1] - crossProd1[1]*crossProd2[0]*middle[2])
        det += (crossProd1[1]*crossProd2[2]*middle[0] + crossProd1[2]*crossProd2[0]*middle[1] - crossProd1[2]*crossProd2[1]*middle[0])
        torsion = (-1)*math.atan2(det, dot)
        return torsion

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

    def centreOnMolecule(self, mol):
        mol_com = np.zeros(3)
        mol_mass = 0
        mol_natoms = 0
        for atom in self.atoms:
            if atom.mol == mol:
                mol_natoms += 1
                mol_mass += atom.mass
                mol_com += atom.mass * atom.coords
        mol_com /= mol_mass

        box = self.box[:, 1] - self.box[:, 0]
        print(box)

        min = np.zeros(3)
        max = np.zeros(3)
        for atom in self.atoms:
            atom.coords -= mol_com
            if atom.mol == mol:
                min = np.minimum(min, atom.coords)
                max = np.maximum(max, atom.coords)
        diag = max - min
        print(diag)

        pbc = np.greater(diag, box/2)
        print(pbc)
        for atom in self.atoms:
            if pbc[0]:
                pass

if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)

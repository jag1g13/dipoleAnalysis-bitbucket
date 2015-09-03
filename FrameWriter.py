from __future__ import print_function
import sys


class FrameWriter(object):
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.filename = filename
        print("Opened file: {0}".format(self.filename))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    def writeFrame(self, num, frame):
        print("ITEM: TIMESTEP", file=self.file)
        print(frame.timestep, file=self.file)
        print("ITEM: NUMBER OF ATOMS", file=self.file)
        print(frame.natoms, file=self.file)
        print("ITEM: BOX BOUNDS pp pp pp", file=self.file)
        print("{0:12.4f}{1:12.4f}".format(frame.box[0, 0], frame.box[0, 1]), file=self.file)
        print("{0:12.4f}{1:12.4f}".format(frame.box[1, 0], frame.box[1, 1]), file=self.file)
        print("{0:12.4f}{1:12.4f}".format(frame.box[2, 0], frame.box[2, 1]), file=self.file)
        print("ITEM: ATOMS id type mol x y z mux muy muz mass diameter", file=self.file)
        for atom in frame.atoms:
            print("{0:5d}{1:5d}{2:5d}{3:10.4f}{4:10.4f}{5:10.4f}{6:10.4f}{7:10.4f}{8:10.4f}{9:8.3f}{10:8.3f}"
                  .format(atom.id, atom.type, atom.mol, atom.coords[0], atom.coords[1], atom.coords[2],
                          atom.dipole[0], atom.dipole[1], atom.dipole[2], atom.mass, atom.diameter), file=self.file)



if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)

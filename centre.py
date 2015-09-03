#!/usr/bin/env python

import sys
import time
from optparse import OptionParser
from FrameReader import FrameReader
from FrameWriter import FrameWriter
from Frame import Frame


def main(filename, outname, nframes=-1, natoms=-1):
    """
    Perform analysis of dipoles in LAMMPS trajectory
    :param lammpstrj: Filename of LAMMPS trajectory
    :return: Number of frames in trajectory
    """
    reader = FrameReader(filename)

    if natoms < 0:
        natoms = reader.total_atoms
    else:
        natoms = min(natoms, reader.total_atoms)

    if nframes < 0:
        nframes = reader.total_frames
    else:
        nframes = min(nframes, reader.total_frames)

    frame = Frame(natoms)
    print(nframes)
    with FrameWriter(outname) as writer:
        for i in xrange(nframes):
            # Read in frame from trajectory and process
            if i % 1000 == 0:
                print(i)
            reader.readFrame(i, frame)
            writer.writeFrame(i, frame)


    return nframes


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input",
                      action="store", type="string", dest="intrj", default="",
                      help="Input file - LAMMPS trajectory", metavar="FILE")
    parser.add_option("-o", "--output",
                      action="store", type="string", dest="outtrj", default="",
                      help="Output file - LAMMPS trajectory", metavar="FILE")
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
    if not options.intrj:
        print("Must provide LAMMPS trajectory to run")
        sys.exit(1)
    t_start = time.clock()
    nframes = main(options.intrj, options.outtrj, options.nframes)
    t_end = time.clock()
    print("-"*25 + "\nCalculated {0} frames in {1}s\n".format(nframes, (t_end - t_start)) + "="*25)

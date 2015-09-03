#!/usr/bin/env python

from __future__ import print_function

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
    with FrameWriter(outname) as writer:
        for i in xrange(nframes):
            # Read in frame from trajectory and process
            progressBar(i, nframes)
            reader.readFrame(i, frame)
            frame.centreOnMolecule(1)
            writer.writeFrame(i, frame)
    return nframes

def progressBar(num, total, length=50, char_done="+", char_remain="-"):
    """
    Print a progress bar
    :param num: Current number of items processed
    :param total: Total number of items to process
    :param length: Length of progress bar - default 50
    :param char_done: Character to use for left of bar
    :param char_remain: Character to use for right of bar
    :return: Nothing
    """
    prog = length * (num+1) / total
    remain = length - prog
    print("\r" + char_done*prog + char_remain*remain, end="")


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
    print("\nCalculated {0} frames in {1}s\n".format(nframes, (t_end - t_start)) + "="*25)

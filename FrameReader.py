import sys
import itertools


class FrameReader(object):
    def __init__(self, filename):
        with open(filename) as file:
            self.all_lines = file.read().splitlines()
        self.total_frames = self.getNumFrames()
        self.total_atoms = self.getNumAtoms()
        print("Read file: {0}".format(filename))
        # self.atom_lines = extractAtomLines_all(self)

    def getNumFrames(self):
        """
        Get the number of frames present in the open trajectory
        :return: Number of frames in trajectory
        """
        timestep_values = self.readTimesteps()    # Takes the values of all timesteps
        return len(timestep_values)

    def getNumAtoms(self):
        """
        Get the number of atoms present in the open trajectory
        :return: Number of atoms in trajectory
        """
        count = 0
        for line in self.all_lines:
            if line[:21] == "ITEM: NUMBER OF ATOMS":
                numOfAtoms = int(self.all_lines[count+1])
                break
            else:
                count += 1
        return numOfAtoms

    def readTimesteps(self):
        """ Function to read the timesteps from the input
        trajectory, so that these can be written into the
        output file.
        lammpstrj - the raw, unedited trajectory. """
        count = 0       # Used to keep track of the lines
        timesteps = []  # Used to store the timestep values
        for line in self.all_lines:
            if line[:14] == "ITEM: TIMESTEP":
                timesteps.append(self.all_lines[count + 1])   # Line after the declaration contains the timestep
                count += 1
                continue
            else:
                count += 1
                continue
        return timesteps    # Length of this list will be equal to the number of frames

    def readAtomCoords(self, frame_lines, frame):
        """ Function to obtain the x,y,z coordinates for
        a particular atom at each frame in the trajectory.
        lammpstrj - raw, unedited lammps trajectory.
        raw_atom_lines - raw, unedited lines, specific to
          one atom. """
        # Function works in the same way as readAtomDipoles(), but written as separate functions for clarity
        start_line = 0

        atom_line = ""
        for i, line in enumerate(frame_lines):
            if line[:11] == "ITEM: ATOMS":
                atom_line = line    # Obtains atom declaration line, as above
                start_line = i+1
                break

        atom_columns = atom_line.split()
        atom_columns = atom_columns[2:]   # Breaks the line into columns and removes the unnecessary first two columns
        column_count = 0
        for header in atom_columns:
            if header == "x":             # Looks for the column containing the x-component of the atom position
                break
            else:
                column_count += 1

        position_x_column = column_count
        # for atom in frame.atoms:
        for i, line in enumerate(frame_lines[start_line:]):
            data_columns = line.split()
            frame.atoms[i].coords[0] = data_columns[position_x_column]
            frame.atoms[i].coords[1] = data_columns[position_x_column+1]
            frame.atoms[i].coords[2] = data_columns[position_x_column+2]
    # No return for the function

    def parseFrame(self, frame_lines, frame):
        """
        Read in frame and atomic data from a single frame of a LAMMPS trajectory
        :param frame_lines: Lines corresponding to a single frame
        :param frame: Frame instance to read data into
        :return: Nothing
        """
        start_line = 0
        for i, line in enumerate(frame_lines):
            if line.startswith("ITEM: TIMESTEP"):
                frame.timestep = int(frame_lines[i+1])
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                if frame.natoms > int(frame_lines[i+1]):
                    print("WARNING: Input file does not fully populate atoms in Frame")
            if line.startswith("ITEM: BOX BOUNDS"):
                for j in xrange(3):
                    vals = frame_lines[i+1+j].split()
                    frame.box[j, 0] = vals[0]
                    frame.box[j, 1] = vals[1]
            if line.startswith("ITEM: ATOMS"):
                start_line = i+1
                break

        atom_cols = frame_lines[start_line-1].split()
        cols = dict(itertools.izip(atom_cols[2:], xrange(len(atom_cols)-2)))

        for i, line in enumerate(frame_lines[start_line:]):
            if i >= frame.natoms:
                break
            data = line.split()
            frame.atoms[i].id = int(data[cols["id"]])
            frame.atoms[i].type = int(data[cols["type"]])
            frame.atoms[i].mol = int(data[cols["mol"]])
            frame.atoms[i].coords[0] = data[cols["x"]]
            frame.atoms[i].coords[1] = data[cols["y"]]
            frame.atoms[i].coords[2] = data[cols["z"]]
            frame.atoms[i].dipole[0] = data[cols["mux"]]
            frame.atoms[i].dipole[1] = data[cols["muy"]]
            frame.atoms[i].dipole[2] = data[cols["muz"]]
            frame.atoms[i].mass = float(data[cols["mass"]])
            frame.atoms[i].diameter = float(data[cols["diameter"]])

    def readAtomDipoles(self, frame_lines, frame):
        """ Function to extract the dipole vectors for each
        atom, so that the positions of the dummy atoms can
        be calculated.
        lammpstrj - the raw, unedited trajectory.
        raw_atom_lines - the raw, unedited lines, specific
          one atom number. """
        start_line = 0

        for i, line in enumerate(frame_lines):
            if line[:11] == "ITEM: ATOMS":
                atom_line = line    # Searches for the line that details the atom columns
                start_line = i+1
                break               # in order to identify the dipole columns

        atom_columns = atom_line.split()  # Breaks up the column headings
        atom_columns = atom_columns[2:]   # Removes the first two elements as they do not correspond to columns of data
        column_count = 0                  # Used to keep track of the column numbers

        for header in atom_columns:
            if header == "mux":           # Looks for the column containing the x-component of the dipole
                break
            else:
                column_count += 1

        dipole_x_column = column_count    # Sets the column number for the dipole x-component
        atom_dipoles = []                 # Defines a list to contain the dipole vectors of the atom at each frame

        for i, line in enumerate(frame_lines[start_line:]):
            data_columns = line.split()
            frame.atoms[i].dipole[0] = data_columns[dipole_x_column]
            frame.atoms[i].dipole[1] = data_columns[dipole_x_column+1]
            frame.atoms[i].dipole[2] = data_columns[dipole_x_column+2]
        # Function doesn't return anything...

    def readFrame(self, frame_num, frame):
        frame_lines = self.extractFrame(frame_num)
        # self.readAtomCoords(frame_lines, frame)
        # self.readAtomDipoles(frame_lines, frame)
        self.parseFrame(frame_lines, frame)

    def filterRemove(raw_atom_lines, atom_type_remove):
        """ Function to remove atoms of a certain atom type
        (e.g. CG water) from the output trajectory.
        raw_atom_lines - raw, unedited lines of all atoms.
        atom_type_remove - atom type to be removed from the
          trajectory. """
        filtered_lines = []   # Defines a list to store the lines that pass through the filter
        removed_IDs = []      # Stores the ID numbers of atoms that have been removed
        remove_num = str(atom_type_remove)    # Converts the atom type to string format

        for line in raw_atom_lines:
            data_columns = line.split()

            if data_columns[1] != remove_num:  # The line passes through the filter if the atom type
                filtered_lines.append(line)    # does not match that which is to be removed
            else:
                removed_IDs.append(data_columns[0])  # Adds the removed ID to the list

        return filtered_lines, removed_IDs   # Returns the set of lines containing the desired atom types

    def extractFrame(self, number):
        """
        Function to extract a specific frame from a trajectory.
        :param number: Frame number to be extracted
        :return: Lines pertaining to the selected frame.
        """
        start = (self.total_atoms + 9) * number
        end = start + self.total_atoms + 9
        return self.all_lines[start:end]

if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)

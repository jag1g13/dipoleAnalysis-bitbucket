import sys

class FrameReader(object):
    def __init__(self, filename):
        file_ = open(filename, 'r')
        print("Read file: {0}".format(filename))
        self.all_lines = file_.readlines()
        file_.close()
        self.total_frames = self.getNumFrames()
        self.total_atoms = self.getNumAtoms()
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

    def readAtomCoords(lammpstrj, raw_atom_lines):
        """ Function to obtain the x,y,z coordinates for
        a particular atom at each frame in the trajectory.
        lammpstrj - raw, unedited lammps trajectory.
        raw_atom_lines - raw, unedited lines, specific to
          one atom. """
        # Function works in the same way as readAtomDipoles(), but written as separate functions for clarity
        for line in lammpstrj:
            if line[:11] == "ITEM: ATOMS":
                atom_line = line    # Obtains atom declaration line, as above
                break
            else:
                continue
        atom_columns = atom_line.split()
        atom_columns = atom_columns[2:]   # Breaks the line into columns and removes the unnecessary first two columns
        column_count = 0
        for header in atom_columns:
            if header == "x":             # Looks for the column containing the x-component of the atom position
                break
            else:
                column_count += 1
                continue
        position_x_column = column_count
        atom_positions = []           # Creates a list to store the position vector of the atom at each frame
        line_count = 0
        for line in raw_atom_lines:
            position_vector = []    # Creates a list to store the vector for this frame
            data_columns = line.split()
            position_vector.append(data_columns[position_x_column])       # Stores the x-component into the position vector
            position_vector.append(data_columns[position_x_column + 1])   # Stores the y-component into the position vector
            position_vector.append(data_columns[position_x_column + 2])   # Stores the z-component into the position vector
            atom_positions.append(position_vector)    # Adds the position vector into the list of position vectors
        return atom_positions    # Returns the list of atom positions


    def readAtomDipoles(lammpstrj, raw_atom_lines):
        """ Function to extract the dipole vectors for each
        atom, so that the positions of the dummy atoms can
        be calculated.
        lammpstrj - the raw, unedited trajectory.
        raw_atom_lines - the raw, unedited lines, specific
          one atom number. """
        for line in frame_lines:
            if line[:11] == "ITEM: ATOMS":
                atom_line = line    # Searches for the line that details the atom columns
                break               # in order to identify the dipole columns
            else:
                continue
        atom_columns = atom_line.split()  # Breaks up the column headings
        atom_columns = atom_columns[2:]   # Removes the first two elements as they do not correspond to columns of data
        column_count = 0                  # Used to keep track of the column numbers
        for header in atom_columns:
            if header == "mux":           # Looks for the column containing the x-component of the dipole
                break
            else:
                column_count += 1
                continue
        dipole_x_column = column_count    # Sets the column number for the dipole x-component
        atom_dipoles = []                 # Defines a list to contain the dipole vectors of the atom at each frame
        for atom in frame.atoms:
            for line in frame_lines:
                data_columns = line.split()
                if data_columns[0] == str(atom.atom_id):
                    atom.dipole[0] = data_columns[dipole_x_column]
                    atom.dipole[1] = data_columns[dipole_x_column+1]
                    atom.dipole[2] = data_columns[dipole_x_column+2]
                else:
                    break
        # Function doesn't return anything...

    
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
                continue
            else:
                removed_IDs.append(data_columns[0])  # Adds the removed ID to the list
                continue
        return filtered_lines, removed_IDs   # Returns the set of lines containing the desired atom types

    def extractFrame(self, number):
        """
        Function to extract a specific frame from a trajectory.
        :param number: Frame number to be extracted
        :return: Lines pertaining to the selected frame.
        """
        # Assuming that the number is counted from 1
        count_line = 0
        count_frame = 0
        for line in self.all_lines:
            count_line += 1
            if line[:14] == "ITEM: TIMESTEP":
                count_frame += 1
                if count == number:
                    break
                continue
            else:
                continue
        count = count_line
        frame_lines = []
        for line in range(count_line, len(self.all_lines)):
            if line[:11] == "ITEM: ATOMS":
                for x in range(self.total_atoms):                        # Upon reaching the atom data section, the next (total_atoms)
                    frame_lines.append(self.all_lines[count + x + 1])     # lines are extracted and stored into a list
        return frame_lines

if __name__ == "__main__":
    print("This file is intended to be imported as a module, not run from the command line")
    sys.exit(0)

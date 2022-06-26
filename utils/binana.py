#######################################################################
# This script was written by Professor Jacob Durrant  -               #
# https://durrantlab.pitt.edu/ and minor modifications have been made #
# by @sammoneykyrle to reduce verbosity and integrate the script into #
# our scoring function. Below is the original credit and license      #
# included in the original script.                                    #
#######################################################################

# BINANA is released under the GNU General Public License (see
# http://www.gnu.org/licenses/gpl.html). If you have any questions, comments,
# or suggestions, please don't hesitate to contact me, Jacob Durrant, at
# jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for
# ligand-binding characterization." J Mol Graph Model 29(6): 888-893.

import __future__

import math
import os
import sys
import textwrap
import json
import io
import pandas as pd

VERSION = "1.3"


"""
Class Point defines a point on the grid
"""
class Point:
    x = 99999.0
    y = 99999.0
    z = 99999.0

    # Initialize nitialize a point
    # Param x (float): x coordinate
    # Param y (float): y coordinate
    # Param z (float): z coordinate
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # Returns a copy of a point
    # Param self (Point): point to be copied
    def copy_of(self):
        return Point(self.x, self.y, self.z)

    # Print the coordinates of a point
    # Param self (Point)
    def print_coors(self):
        print((str(self.x) + "\t" + str(self.y) + "\t" + str(self.z)))

    def snap(self, reso):  # snap the point to a grid
        self.x = round(self.x / reso) * reso
        self.y = round(self.y / reso) * reso
        self.z = round(self.z / reso) * reso

    # Returns the distance between two points
    # Param self (Point): this point
    # Param a_point (Point): the other point
    def dist_to(self, apoint):
        return math.sqrt(
            math.pow(self.x - apoint.x, 2)
            + math.pow(self.y - apoint.y, 2)
            + math.pow(self.z - apoint.z, 2)
        )

    # Returns a the coordinates of a point
    # Param self (Point)
    def description(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z)

    # Returns the magnitude of a point (distance from origin)
    # Param self (Point)
    def magnitude(self):
        return self.dist_to(Point(0, 0, 0))

    # Returns a PDB line for the point
    # Param self (Point)
    # Param index (integer): index of the point
    def create_PDB_line(self, index):
        output = "ATOM "
        output = output + str(index).rjust(6) + "X".rjust(5) + "XXX".rjust(4)
        output = output + ("%.3f" % self.x).rjust(18)
        output = output + ("%.3f" % self.y).rjust(8)
        output = output + ("%.3f" % self.z).rjust(8)
        output = output + "X".rjust(24)
        return output


"""
Class Atom defines an atom
"""


class Atom:
    # Initilaize an atom
    def __init__(self):
        self.atom_name = ""
        self.residue = ""
        self.coordinates = Point(99999, 99999, 99999)
        self.element = ""
        self.PDB_index = ""
        self.line = ""
        self.atom_type = ""
        self.indecies_of_atoms_connecting = []
        self.charge = 0
        self.resid = 0
        self.chain = ""
        self.structure = ""
        self.comment = ""

    # Returns a copy of an atom
    # Param self (Atom): atom to be copied
    def copy_of(self):
        theatom = Atom()
        theatom.atom_name = self.atom_name
        theatom.residue = self.residue
        theatom.coordinates = self.coordinates.copy_of()
        theatom.element = self.element
        theatom.PDB_index = self.PDB_index
        theatom.line = self.line
        theatom.atom_type = self.atom_type
        theatom.indecies_of_atoms_connecting = self.indecies_of_atoms_connecting[:]
        theatom.charge = self.charge
        theatom.resid = self.resid
        theatom.chain = self.chain
        theatom.structure = self.structure
        theatom.comment = self.comment

        return theatom

    # Returns a string identifying/describing the atom
    # Param self (Atom)
    def string_id(self):
        to_return = ""
        if self.chain.strip() != "":
            to_return = to_return + self.chain.strip() + ":"
        to_return = (
            to_return
            + self.residue.strip()
            + "("
            + str(self.resid)
            + "):"
            + self.atom_name.strip()
            + "("
            + str(self.PDB_index)
            + ")"
        )
        return to_return

    # Returns a PDB line for the atom
    # Param self (Atom)
    # Param index (integer): index of the point
    def create_PDB_line(self, index):
        output = "ATOM "
        output = output + str(index).rjust(6) + self.atom_name.rjust(5) + self.residue.rjust(4)
        output = output + ("%.3f" % self.coordinates.x).rjust(18)
        output = output + ("%.3f" % self.coordinates.y).rjust(8)
        output = output + ("%.3f" % self.coordinates.z).rjust(8)
        output = output + self.element.rjust(24)
        return output

    # Returns the number of an atom's nearest neighbors
    # Param self (Atom)
    def number_of_neighbors(self):
        return len(self.indecies_of_atoms_connecting)

    # Adds the index of an atom's neibor to the list of connecting atoms
    # Param self (Atom)
    # Param index (float): index of atom to be added
    def add_neighbor_atom_index(self, index):
        if not (index in self.indecies_of_atoms_connecting):
            self.indecies_of_atoms_connecting.append(index)

    # Returns whether the atom is part of the backbone or a sidechain
    # Only really applies to proteins, assuming standard atom names
    # Param self (Atom)
    def SideChainOrBackBone(self):
        if (
            self.atom_name.strip() == "CA"
            or self.atom_name.strip() == "C"
            or self.atom_name.strip() == "O"
            or self.atom_name.strip() == "N"
        ):
            return "BACKBONE"
        else:
            return "SIDECHAIN"

    # Reads name of atom in from a PDB line
    # Param self (Atom)
    # Param line (string): PDB line
    def read_PDB_line(self, line):
        self.line = line
        self.atom_name = line[11:16].strip()
        # Read atom name
        if len(self.atom_name) == 1:
            self.atom_name = self.atom_name + "  "
        elif len(self.atom_name) == 2:
            self.atom_name = self.atom_name + " "
        elif len(self.atom_name) == 3:
            # This line is necessary for babel to work, though many PDBs in
            # the PDB would have this line commented out
            self.atom_name = self.atom_name + " "

        self.coordinates = Point(
            float(line[30:38]), float(line[38:46]), float(line[46:54])
        )

        # Now read in atom type (for pdbqt)
        self.atom_type = line[77:79].strip().upper()

        # Read the atom's charge
        if line[69:76].strip() != "":
            self.charge = float(line[69:76])
        else:
            self.charge = 0.0

        # Try to guess at element from name
        if self.element == "":
            two_letters = self.atom_name[0:2].strip().upper()
            if two_letters == "BR":
                self.element = "BR"
            elif two_letters == "CL":
                self.element = "CL"
            elif two_letters == "BI":
                self.element = "BI"
            elif two_letters == "AS":
                self.element = "AS"
            elif two_letters == "AG":
                self.element = "AG"
            elif two_letters == "LI":
                self.element = "LI"
            # elif two_letters=='HG':
            #    self.element='HG'
            elif two_letters == "MG":
                self.element = "MG"
            elif two_letters == "MN":
                self.element = "MN"
            elif two_letters == "RH":
                self.element = "RH"
            elif two_letters == "ZN":
                self.element = "ZN"
            elif two_letters == "FE":
                self.element = "FE"
            else:
                # So, just assume it's the first letter. Any number needs to
                # be removed from the element name
                self.element = self.atom_name
                self.element = self.element.replace("0", "")
                self.element = self.element.replace("1", "")
                self.element = self.element.replace("2", "")
                self.element = self.element.replace("3", "")
                self.element = self.element.replace("4", "")
                self.element = self.element.replace("5", "")
                self.element = self.element.replace("6", "")
                self.element = self.element.replace("7", "")
                self.element = self.element.replace("8", "")
                self.element = self.element.replace("9", "")
                self.element = self.element.replace("@", "")

                self.element = self.element[0:1].strip().upper()

        self.PDB_index = line[6:12].strip()
        self.residue = line[16:20]
        # This only uses the rightmost three characters, essentially removing
        # unique rotamer identification
        self.residue = " " + self.residue[-3:]

        # It's possible the pdbqt might not have any resid entries.
        try:
            self.resid = int(line[23:26])
        except:
            pass

        self.chain = line[21:22]
        if self.residue.strip() == "":
            self.residue = " MOL"


"""
Class PDB handles PDB filing
"""


class PDB:

    # Initialize PDB
    def __init__(self):
        self.all_atoms = {}
        self.non_protein_atoms = {}
        self.max_x = -9999.99
        self.min_x = 9999.99
        self.max_y = -9999.99
        self.min_y = 9999.99
        self.max_z = -9999.99
        self.min_z = 9999.99
        self.rotateable_bonds_count = 0
        self.functions = MathFunctions()
        self.protein_resnames = [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "ASH",
            "ASX",
            "CYS",
            "CYM",
            "CYX",
            "GLN",
            "GLU",
            "GLH",
            "GLX",
            "GLY",
            "HIS",
            "HID",
            "HIE",
            "HIP",
            "ILE",
            "LEU",
            "LYS",
            "LYN",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ]
        self.aromatic_rings = []
        self.charges = []  # a list of points

    # Load PDB file
    # Param file_name (string)
    # Param min_x (float): minimum x coordinate
    # Param max_x (float): maximum x coordinate
    # Param min_y (float): minimum y coordinate
    # Param max_y (float): maximum y coordinate
    # Param min_z (float): minimum z coordinate
    # Param max_z (float): maximum z coordinate
    def load_PDB(
        self,
        file_name,
        min_x=-9999.99,
        max_x=9999.99,
        min_y=-9999.99,
        max_y=9999.99,
        min_z=-9999.99,
        max_z=9999.99,
    ):

        autoindex = 1

        self.__init__()

        # Now load the file into a list
        if 'REMARK' in file_name:
            lines = file_name.split('\n')
        else:
            file = open(file_name, "r")
            lines = file.readlines()
            file.close()

        # Keep track of atomname_resid_chain pairs, to make sure redundants
        # aren't loaded This basically gets rid of rotamers, I think.
        atom_already_loaded = []

        for t in range(0, len(lines)):
            line = lines[t]

            if line[:3] == "END" and line[:7] != "ENDROOT" and line[:9] != "ENDBRANCH":
                t = textwrap.wrap(
                    "WARNING: END or ENDMDL term found in "
                    + file_name
                    + ". Everything after the first instance of this term will be ignored. \
                    If any of your PDBQT files have multiple frames/poses, please partition them \
                    into separate files using vina_split and feed each of the the single-frame files into Binana separately.",
                    80,
                )
                print(("\n".join(t) + "\n"))
                print(line)
                break

            if "between atoms" in line and " A " in line:
                self.rotateable_bonds_count = self.rotateable_bonds_count + 1

            if len(line) >= 7:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    # Load atom data (coordinates, etc.)
                    TempAtom = Atom()
                    TempAtom.read_PDB_line(line)

                    if (
                        TempAtom.coordinates.x > min_x
                        and TempAtom.coordinates.x < max_x
                        and TempAtom.coordinates.y > min_y
                        and TempAtom.coordinates.y < max_y
                        and TempAtom.coordinates.z > min_z
                        and TempAtom.coordinates.z < max_z
                    ):

                        if self.max_x < TempAtom.coordinates.x:
                            self.max_x = TempAtom.coordinates.x
                        if self.max_y < TempAtom.coordinates.y:
                            self.max_y = TempAtom.coordinates.y
                        if self.max_z < TempAtom.coordinates.z:
                            self.max_z = TempAtom.coordinates.z

                        if self.min_x > TempAtom.coordinates.x:
                            self.min_x = TempAtom.coordinates.x
                        if self.min_y > TempAtom.coordinates.y:
                            self.min_y = TempAtom.coordinates.y
                        if self.min_z > TempAtom.coordinates.z:
                            self.min_z = TempAtom.coordinates.z

                        # this string uniquely identifies each atom
                        key = (
                            TempAtom.atom_name.strip()
                            + "_"
                            + str(TempAtom.resid)
                            + "_"
                            + TempAtom.residue.strip()
                            + "_"
                            + TempAtom.chain.strip()
                        )

                        if (
                            key in atom_already_loaded
                            and TempAtom.residue.strip() in self.protein_resnames
                        ):
                            # so this is a protein atom that has already been
                            # loaded once
                            '''
                            self.printout(
                                'Warning: Duplicate protein atom detected: "'
                                + TempAtom.line.strip()
                                + '". Not loading this duplicate.'
                            )
                            print("")
                            '''
                        if (
                            key not in atom_already_loaded
                            or not TempAtom.residue.strip() in self.protein_resnames
                        ):
                            # So either the atom hasn't been loaded, or else
                            # it's a non-protein atom So note that non-protein
                            # atoms can have redundant names, but protein
                            # atoms cannot. This is because protein residues
                            # often contain rotamers
                            atom_already_loaded.append(
                                key
                            )  # So each atom can only be loaded once. No rotamers.
                            self.all_atoms[
                                autoindex
                            ] = (
                                TempAtom
                            )  # So you're actually reindexing everything here.
                            if TempAtom.residue[-3:] not in self.protein_resnames:
                                self.non_protein_atoms[autoindex] = TempAtom

                            autoindex = autoindex + 1

        self.check_protein_format()

        # Only for the ligand, because bonds can be inferred based on
        # atomnames from PDB
        self.create_bonds_by_distance()

        self.assign_aromatic_rings()
        self.assign_charges()

    # Print the PDB line
    # Param self (PDB)
    # Param the_string (string)
    def printout(self, the_string):
        lines = textwrap.wrap(the_string, 80)
        for line in lines:
            print(line)

    # Write and save PDB line to a file
    # Param self (PDB)
    # Param file_name (string)
    def save_PDB(self, file_name):
        f = open(file_name, "w")
        towrite = self.save_PDB_String()
        if towrite.strip() == "":
            # So no PDB is empty, VMD will load them all
            towrite = "ATOM      1  X   XXX             0.000   0.000   0.000                       X"
        f.write(towrite)
        f.close()

    # Returns a new PDB line
    # Param self (PDB)
    def save_PDB_String(self):
        to_output = ""

        # Write coordinates of all atoms
        for atom_index in self.all_atoms:
            to_output = (
                to_output
                + self.all_atoms[atom_index].create_PDB_line(atom_index)
                + "\n"
            )

        return to_output

    # Adds a new atom to this PDB
    # Param self (PDB)
    # Param atom (Atom): new atom being added
    def add_new_atom(self, atom):
        # first get available index
        t = 1
        while t in list(self.all_atoms.keys()):
            t = t + 1

        # now add atom
        self.all_atoms[t] = atom

    # Assign residue name to atom
    # Param self (PDB)
    # Param  resname (string): residue name
    def set_resname(self, resname):
        for atom_index in self.all_atoms:
            self.all_atoms[atom_index].residue = resname

    # Returns a list of the indeces of atoms connected to a given element
    # Param self (PDB)
    # Param index (integer): index of atom
    # Param connected_atom_element (string): element in question
    def connected_atoms_of_given_element(self, index, connected_atom_element):
        atom = self.all_atoms[index]
        connected_atoms = []
        for index2 in atom.indecies_of_atoms_connecting:
            atom2 = self.all_atoms[index2]
            if atom2.element == connected_atom_element:
                connected_atoms.append(index2)
        return connected_atoms

    # Returns a list of the indeces of heavy atoms connected at specified index
    # Param self (PDB)
    # Param index (integer): index of atom
    def connected_heavy_atoms(self, index):
        atom = self.all_atoms[index]
        connected_atoms = []
        for index2 in atom.indecies_of_atoms_connecting:
            atom2 = self.all_atoms[index2]
            if atom2.element != "H":
                connected_atoms.append(index2)
        return connected_atoms

    # Correct format of the protein
    # Param self (PDB)
    def check_protein_format(self):
        curr_res = ""
        first = True
        residue = []

        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]

            key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

            if first is True:
                curr_res = key
                first = False

            if key != curr_res:

                self.check_protein_format_process_residue(residue, last_key)

                residue = []
                curr_res = key

            residue.append(atom.atom_name.strip())
            last_key = key

        self.check_protein_format_process_residue(residue, last_key)

    # Correct format of the protein and residues
    # Param self (PDB)
    # Param residue ()
    # Param last_key ()
    def check_protein_format_process_residue(self, residue, last_key):
        temp = last_key.strip().split("_")
        resname = temp[0]
        real_resname = resname[-3:]
        resid = temp[1]
        chain = temp[2]
        '''
        if real_resname in self.protein_resnames:  # so it's a protein residue
            if "N" not in residue:
                self.printout(
                    'Warning: There is no atom named "N" in the protein residue '
                    + last_key
                    + ". Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore."
                )
                print("")
            if "C" not in residue:
                self.printout(
                    'Warning: There is no atom named "C" in the protein residue '
                    + last_key
                    + ". Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore."
                )
                print("")
            if "CA" not in residue:
                self.printout(
                    'Warning: There is no atom named "CA" in the protein residue '
                    + last_key
                    + ". Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore."
                )
                print("")

            if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
                if "OE1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "OE1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "OE2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "OE2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if real_resname == "ASP" or real_resname == "ASH" or real_resname == "ASX":
                if "OD1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "OD1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "OD2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "OD2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if real_resname == "LYS" or real_resname == "LYN":
                if "NZ" not in residue:
                    self.printout(
                        'Warning: There is no atom named "NZ" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if real_resname == "ARG":
                if "NH1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "NH1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "NH2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "NH2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if (
                real_resname == "HIS"
                or real_resname == "HID"
                or real_resname == "HIE"
                or real_resname == "HIP"
            ):
                if "NE2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "NE2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "ND1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "ND1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if real_resname == "PHE":
                if "CG" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CG" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CZ" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CZ" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if real_resname == "TYR":
                if "CG" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CG" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CZ" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CZ" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if real_resname == "TRP":
                if "CG" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CG" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "NE1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "NE1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE3" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE3" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CZ2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CZ2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CZ3" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CZ3" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CH2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CH2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")

            if (
                real_resname == "HIS"
                or real_resname == "HID"
                or real_resname == "HIE"
                or real_resname == "HIP"
            ):
                if "CG" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CG" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "ND1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "ND1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CD2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CD2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "CE1" not in residue:
                    self.printout(
                        'Warning: There is no atom named "CE1" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
                if "NE2" not in residue:
                    self.printout(
                        'Warning: There is no atom named "NE2" in the protein residue '
                        + last_key
                        + ". Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore."
                    )
                    print("")
        '''
    # Functions to determine the bond connectivity based on distance
    # ==============================================================

    # Define bonds between atoms using distance on the grid
    # Param self (PDB)
    def create_bonds_by_distance(self):
        for AtomIndex1 in self.non_protein_atoms:
            atom1 = self.non_protein_atoms[AtomIndex1]
            if atom1.residue[-3:] not in self.protein_resnames:
                # so it's not a protein residue
                for AtomIndex2 in self.non_protein_atoms:
                    if AtomIndex1 != AtomIndex2:
                        atom2 = self.non_protein_atoms[AtomIndex2]
                        if not atom2.residue[-3:] in self.protein_resnames:
                            # so it's not a protein residue
                            dist = self.functions.distance(
                                atom1.coordinates, atom2.coordinates
                            )

                            if (
                                dist
                                < self.bond_length(atom1.element, atom2.element) * 1.2
                            ):
                                atom1.add_neighbor_atom_index(AtomIndex2)
                                atom2.add_neighbor_atom_index(AtomIndex1)

    # Retuns bond length between two elements
    # Param self (PDB)
    # Param element1 (string): symbol of first element
    # Param element2 (string): symbol of second element
    def bond_length(self, element1, element2):

        """Bond lengths taken from Handbook of Chemistry and Physics. The information provided there was very specific,
        so I tried to pick representative examples and used the bond lengths from those. Sitautions could arise where these
        lengths would be incorrect, probably slight errors (<0.06) in the hundreds."""

        distance = 0.0
        if element1 == "C" and element2 == "C":
            distance = 1.53
        if element1 == "N" and element2 == "N":
            distance = 1.425
        if element1 == "O" and element2 == "O":
            distance = 1.469
        if element1 == "S" and element2 == "S":
            distance = 2.048
        if (element1 == "C" and element2 == "H") or (
            element1 == "H" and element2 == "C"
        ):
            distance = 1.059
        if (element1 == "C" and element2 == "N") or (
            element1 == "N" and element2 == "C"
        ):
            distance = 1.469
        if (element1 == "C" and element2 == "O") or (
            element1 == "O" and element2 == "C"
        ):
            distance = 1.413
        if (element1 == "C" and element2 == "S") or (
            element1 == "S" and element2 == "C"
        ):
            distance = 1.819
        if (element1 == "N" and element2 == "H") or (
            element1 == "H" and element2 == "N"
        ):
            distance = 1.009
        if (element1 == "N" and element2 == "O") or (
            element1 == "O" and element2 == "N"
        ):
            distance = 1.463
        if (element1 == "O" and element2 == "S") or (
            element1 == "S" and element2 == "O"
        ):
            distance = 1.577
        if (element1 == "O" and element2 == "H") or (
            element1 == "H" and element2 == "O"
        ):
            distance = 0.967
        if (element1 == "S" and element2 == "H") or (
            element1 == "H" and element2 == "S"
        ):
            # This one not from source sited above. Not sure where it's from,
            # but it wouldn't ever be used in the current context ("AutoGrow")
            distance = 2.025 / 1.5
        if (element1 == "S" and element2 == "N") or (
            element1 == "H" and element2 == "N"
        ):
            distance = 1.633

        if (element1 == "C" and element2 == "F") or (
            element1 == "F" and element2 == "C"
        ):
            distance = 1.399
        if (element1 == "C" and element2 == "CL") or (
            element1 == "CL" and element2 == "C"
        ):
            distance = 1.790
        if (element1 == "C" and element2 == "BR") or (
            element1 == "BR" and element2 == "C"
        ):
            distance = 1.910
        if (element1 == "C" and element2 == "I") or (
            element1 == "I" and element2 == "C"
        ):
            distance = 2.162

        if (element1 == "S" and element2 == "BR") or (
            element1 == "BR" and element2 == "S"
        ):
            distance = 2.321
        if (element1 == "S" and element2 == "CL") or (
            element1 == "CL" and element2 == "S"
        ):
            distance = 2.283
        if (element1 == "S" and element2 == "F") or (
            element1 == "F" and element2 == "S"
        ):
            distance = 1.640
        if (element1 == "S" and element2 == "I") or (
            element1 == "I" and element2 == "S"
        ):
            distance = 2.687

        if (element1 == "P" and element2 == "BR") or (
            element1 == "BR" and element2 == "P"
        ):
            distance = 2.366
        if (element1 == "P" and element2 == "CL") or (
            element1 == "CL" and element2 == "P"
        ):
            distance = 2.008
        if (element1 == "P" and element2 == "F") or (
            element1 == "F" and element2 == "P"
        ):
            distance = 1.495
        if (element1 == "P" and element2 == "I") or (
            element1 == "I" and element2 == "P"
        ):
            distance = 2.490
        if (element1 == "P" and element2 == "O") or (
            element1 == "O" and element2 == "P"
        ):
            # estimate based on eye balling Handbook of Chemistry and Physics
            distance = 1.6

        if (element1 == "N" and element2 == "BR") or (
            element1 == "BR" and element2 == "N"
        ):
            distance = 1.843
        if (element1 == "N" and element2 == "CL") or (
            element1 == "CL" and element2 == "N"
        ):
            distance = 1.743
        if (element1 == "N" and element2 == "F") or (
            element1 == "F" and element2 == "N"
        ):
            distance = 1.406
        if (element1 == "N" and element2 == "I") or (
            element1 == "I" and element2 == "N"
        ):
            distance = 2.2

        if (element1 == "SI" and element2 == "BR") or (
            element1 == "BR" and element2 == "SI"
        ):
            distance = 2.284
        if (element1 == "SI" and element2 == "CL") or (
            element1 == "CL" and element2 == "SI"
        ):
            distance = 2.072
        if (element1 == "SI" and element2 == "F") or (
            element1 == "F" and element2 == "SI"
        ):
            distance = 1.636
        if (element1 == "SI" and element2 == "P") or (
            element1 == "P" and element2 == "SI"
        ):
            distance = 2.264
        if (element1 == "SI" and element2 == "S") or (
            element1 == "S" and element2 == "SI"
        ):
            distance = 2.145
        if (element1 == "SI" and element2 == "SI") or (
            element1 == "SI" and element2 == "SI"
        ):
            distance = 2.359
        if (element1 == "SI" and element2 == "C") or (
            element1 == "C" and element2 == "SI"
        ):
            distance = 1.888
        if (element1 == "SI" and element2 == "N") or (
            element1 == "N" and element2 == "SI"
        ):
            distance = 1.743
        if (element1 == "SI" and element2 == "O") or (
            element1 == "O" and element2 == "SI"
        ):
            distance = 1.631

        return distance

    # Functions to identify positive charges
    # ======================================

    # Assign Charges to atoms in protein
    # Param self (PDB)
    def assign_charges(self):
        # Get all the quartinary amines on non-protein residues (these are the
        # only non-protein groups that will be identified as positively
        # charged)
        all_charged = []
        for atom_index in self.non_protein_atoms:
            atom = self.non_protein_atoms[atom_index]
            if (
                atom.element == "MG"
                or atom.element == "MN"
                or atom.element == "RH"
                or atom.element == "ZN"
                or atom.element == "FE"
                or atom.element == "BI"
                or atom.element == "AS"
                or atom.element == "AG"
            ):
                chrg = self.Charged(atom.coordinates, [atom_index], True)
                self.charges.append(chrg)

            if atom.element == "N":
                if atom.number_of_neighbors() == 4:
                    # a quartinary amine, so it's easy
                    indexes = [atom_index]
                    indexes.extend(atom.indecies_of_atoms_connecting)

                    # so the indicies stored is just the index of the nitrogen
                    # and any attached atoms
                    chrg = self.Charged(atom.coordinates, indexes, True)
                    self.charges.append(chrg)
                elif atom.number_of_neighbors() == 3:
                    # maybe you only have two hydrogen's added, by they're sp3
                    # hybridized. Just count this as a quartinary amine, since
                    # I think the positive charge would be stabalized.
                    nitrogen = atom
                    atom1 = self.all_atoms[atom.indecies_of_atoms_connecting[0]]
                    atom2 = self.all_atoms[atom.indecies_of_atoms_connecting[1]]
                    atom3 = self.all_atoms[atom.indecies_of_atoms_connecting[2]]
                    angle1 = (
                        self.functions.angle_between_three_points(
                            atom1.coordinates, nitrogen.coordinates, atom2.coordinates
                        )
                        * 180.0
                        / math.pi
                    )
                    angle2 = (
                        self.functions.angle_between_three_points(
                            atom1.coordinates, nitrogen.coordinates, atom3.coordinates
                        )
                        * 180.0
                        / math.pi
                    )
                    angle3 = (
                        self.functions.angle_between_three_points(
                            atom2.coordinates, nitrogen.coordinates, atom3.coordinates
                        )
                        * 180.0
                        / math.pi
                    )
                    average_angle = (angle1 + angle2 + angle3) / 3
                    if math.fabs(average_angle - 109.0) < 5.0:
                        indexes = [atom_index]
                        indexes.extend(atom.indecies_of_atoms_connecting)
                        # so indexes added are the nitrogen and any attached
                        # atoms.
                        chrg = self.Charged(nitrogen.coordinates, indexes, True)
                        self.charges.append(chrg)

            if atom.element == "C":
                # let's check for guanidino-like groups (actually H2N-C-NH2,
                # where not CN3.)
                if atom.number_of_neighbors() == 3:
                    # the carbon has only three atoms connected to it
                    nitrogens = self.connected_atoms_of_given_element(atom_index, "N")
                    if len(nitrogens) >= 2:
                        # so carbon is connected to at least two nitrogens.
                        # now we need to count the number of nitrogens that
                        # are only connected to one heavy atom (the carbon).
                        nitrogens_to_use = []
                        all_connected = atom.indecies_of_atoms_connecting[:]
                        not_isolated = -1

                        for atmindex in nitrogens:
                            if len(self.connected_heavy_atoms(atmindex)) == 1:
                                nitrogens_to_use.append(atmindex)
                                all_connected.remove(atmindex)

                        if len(all_connected) > 0:
                            # get the atom that connects this charged group to
                            # the rest of the molecule, ultimately to make
                            # sure it's sp3 hybridized
                            not_isolated = all_connected[0]

                        if len(nitrogens_to_use) == 2 and not_isolated != -1:
                            # so there are at two nitrogens that are only
                            # connected to the carbon (and probably some
                            # hydrogens)

                            # now you need to make sure not_isolated atom is
                            # sp3 hybridized
                            not_isolated_atom = self.all_atoms[not_isolated]
                            if (
                                (
                                    not_isolated_atom.element == "C"
                                    and not_isolated_atom.number_of_neighbors() == 4
                                )
                                or (
                                    not_isolated_atom.element == "O"
                                    and not_isolated_atom.number_of_neighbors() == 2
                                )
                                or not_isolated_atom.element == "N"
                                or not_isolated_atom.element == "S"
                                or not_isolated_atom.element == "P"
                            ):

                                pt = self.all_atoms[
                                    nitrogens_to_use[0]
                                ].coordinates.copy_of()
                                pt.x = (
                                    pt.x
                                    + self.all_atoms[nitrogens_to_use[1]].coordinates.x
                                )
                                pt.y = (
                                    pt.y
                                    + self.all_atoms[nitrogens_to_use[1]].coordinates.y
                                )
                                pt.z = (
                                    pt.z
                                    + self.all_atoms[nitrogens_to_use[1]].coordinates.z
                                )
                                pt.x = pt.x / 2.0
                                pt.y = pt.y / 2.0
                                pt.z = pt.z / 2.0

                                indexes = [atom_index]
                                indexes.extend(nitrogens_to_use)
                                indexes.extend(
                                    self.connected_atoms_of_given_element(
                                        nitrogens_to_use[0], "H"
                                    )
                                )
                                indexes.extend(
                                    self.connected_atoms_of_given_element(
                                        nitrogens_to_use[1], "H"
                                    )
                                )

                                chrg = self.Charged(
                                    pt, indexes, True  # True because it's positive
                                )
                                self.charges.append(chrg)

            if atom.element == "C":
                # let's check for a carboxylate
                if atom.number_of_neighbors() == 3:
                    # a carboxylate carbon will have three items connected to
                    # it.
                    oxygens = self.connected_atoms_of_given_element(atom_index, "O")
                    if len(oxygens) == 2:
                        # a carboxylate will have two oxygens connected to it.
                        # now, each of the oxygens should be connected to only
                        # one heavy atom (so if it's connected to a hydrogen,
                        # that's okay)
                        if (
                            len(self.connected_heavy_atoms(oxygens[0])) == 1
                            and len(self.connected_heavy_atoms(oxygens[1])) == 1
                        ):
                            # so it's a carboxylate! Add a negative charge.
                            pt = self.all_atoms[oxygens[0]].coordinates.copy_of()
                            pt.x = pt.x + self.all_atoms[oxygens[1]].coordinates.x
                            pt.y = pt.y + self.all_atoms[oxygens[1]].coordinates.y
                            pt.z = pt.z + self.all_atoms[oxygens[1]].coordinates.z
                            pt.x = pt.x / 2.0
                            pt.y = pt.y / 2.0
                            pt.z = pt.z / 2.0
                            chrg = self.Charged(
                                pt, [oxygens[0], atom_index, oxygens[1]], False
                            )
                            self.charges.append(chrg)

            if atom.element == "P":
                # let's check for a phosphate or anything where a phosphorus
                # is bound to two oxygens where both oxygens are bound to only
                # one heavy atom (the phosphorus). I think this will get
                # several phosphorus substances.
                oxygens = self.connected_atoms_of_given_element(atom_index, "O")
                if len(oxygens) >= 2:
                    # the phosphorus is bound to at least two oxygens now
                    # count the number of oxygens that are only bound to the
                    # phosphorus
                    count = 0
                    for oxygen_index in oxygens:
                        if len(self.connected_heavy_atoms(oxygen_index)) == 1:
                            count = count + 1
                    if count >= 2:
                        # so there are at least two oxygens that are only
                        # bound to the phosphorus
                        indexes = [atom_index]
                        indexes.extend(oxygens)
                        chrg = self.Charged(atom.coordinates, indexes, False)
                        self.charges.append(chrg)

            if atom.element == "S":
                # let's check for a sulfonate or anything where a sulfur is
                # bound to at least three oxygens and at least three are bound
                # to only the sulfur (or the sulfur and a hydrogen).
                oxygens = self.connected_atoms_of_given_element(atom_index, "O")
                if len(oxygens) >= 3:
                    # the sulfur is bound to at least three oxygens. now count
                    # the number of oxygens that are only bound to the sulfur
                    count = 0
                    for oxygen_index in oxygens:
                        if len(self.connected_heavy_atoms(oxygen_index)) == 1:
                            count = count + 1
                    if count >= 3:
                        # so there are at least three oxygens that are only
                        # bound to the sulfur
                        indexes = [atom_index]
                        indexes.extend(oxygens)
                        chrg = self.Charged(atom.coordinates, indexes, False)
                        self.charges.append(chrg)

        # Now that you've found all the positive charges in non-protein
        # residues, it's time to look for aromatic rings in protein residues
        curr_res = ""
        first = True
        residue = []

        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]

            key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

            if first == True:
                curr_res = key
                first = False

            if key != curr_res:

                self.assign_charged_from_protein_process_residue(residue, last_key)

                residue = []
                curr_res = key

            residue.append(atom_index)
            last_key = key

        self.assign_charged_from_protein_process_residue(residue, last_key)

    # Assign charges but with protein residue
    # Param self (PDB)
    # Param residue ()
    # Param last_key ()
    def assign_charged_from_protein_process_residue(self, residue, last_key):
        temp = last_key.strip().split("_")
        resname = temp[0]
        real_resname = resname[-3:]
        resid = temp[1]
        chain = temp[2]

        if real_resname == "LYS" or real_resname == "LYN":
            # regardless of protonation state, assume it's charged.
            for index in residue:
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "NZ":

                    # quickly go through the residue and get the hydrogens
                    # attached to this nitrogen to include in the index list
                    indexes = [index]
                    for index2 in residue:
                        atom2 = self.all_atoms[index2]
                        if atom2.atom_name.strip() == "HZ1":
                            indexes.append(index2)
                        if atom2.atom_name.strip() == "HZ2":
                            indexes.append(index2)
                        if atom2.atom_name.strip() == "HZ3":
                            indexes.append(index2)

                    chrg = self.Charged(atom.coordinates, indexes, True)
                    self.charges.append(chrg)
                    break

        if real_resname == "ARG":
            charge_pt = Point(0.0, 0.0, 0.0)
            count = 0.0
            indices = []
            for index in residue:
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "NH1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "NH2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "2HH2":
                    indices.append(index)
                if atom.atom_name.strip() == "1HH2":
                    indices.append(index)
                if atom.atom_name.strip() == "CZ":
                    indices.append(index)
                if atom.atom_name.strip() == "2HH1":
                    indices.append(index)
                if atom.atom_name.strip() == "1HH1":
                    indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count

                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.Charged(charge_pt, indices, True)
                    self.charges.append(chrg)

        if real_resname in ["HIS", "HID", "HIE", "HIP"]:
            # regardless of protonation state, assume it's charged. This based
            # on "The Cation-Pi Interaction," which suggests protonated state
            # would be stabalized. But let's not consider HIS when doing salt
            # bridges.
            charge_pt = Point(0.0, 0.0, 0.0)
            count = 0.0
            indices = []
            for index in residue:
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "NE2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "ND1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "HE2":
                    indices.append(index)
                if atom.atom_name.strip() == "HD1":
                    indices.append(index)
                if atom.atom_name.strip() == "CE1":
                    indices.append(index)
                if atom.atom_name.strip() == "CD2":
                    indices.append(index)
                if atom.atom_name.strip() == "CG":
                    indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count
                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.Charged(charge_pt, indices, True)
                    self.charges.append(chrg)

        if real_resname in ["GLU", "GLH", "GLX"]:
            # regardless of protonation state, assume it's charged. This based
            # on "The Cation-Pi Interaction," which suggests protonated state
            # would be stabalized.
            charge_pt = Point(0.0, 0.0, 0.0)
            count = 0.0
            indices = []
            for index in residue:
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "OE1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "OE2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "CD":
                    indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count
                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.Charged(
                        charge_pt,
                        indices,
                        False,  # False because it's a negative charge
                    )
                    self.charges.append(chrg)

        if real_resname in ["ASP", "ASH", "ASX"]:
            # regardless of protonation state, assume it's charged. This based
            # on "The Cation-Pi Interaction," which suggests protonated state
            # would be stabalized.
            charge_pt = Point(0.0, 0.0, 0.0)
            count = 0.0
            indices = []
            for index in residue:
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "OD1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "OD2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atom_name.strip() == "CG":
                    indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count
                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.Charged(
                        charge_pt,
                        indices,
                        False,  # False because it's a negative charge
                    )
                    self.charges.append(chrg)

    """
    Class Charged defines the charge of atom
    """

    class Charged:
        # Initialize charge
        # Param self (Charged)
        # Param coordinates (Point): point on grid
        # Param indecies (list): indecies of atom
        # Param positive (boolean): True if atom is positively charged
        def __init__(self, coordinates, indices, positive):
            self.coordinates = coordinates
            self.indices = indices
            # true or false to specifiy if positive or negative charge
            self.positive = positive

    # Functions to identify aromatic rings
    # ====================================

    # Marks atoms present in an aromatic ring
    # Param indicies_of_ring (list): indecies of atoms in ring
    def add_aromatic_marker(self, indicies_of_ring):
        # first identify the center point
        points_list = []
        total = len(indicies_of_ring)
        x_sum = 0.0
        y_sum = 0.0
        z_sum = 0.0

        for index in indicies_of_ring:
            atom = self.all_atoms[index]
            points_list.append(atom.coordinates)
            x_sum = x_sum + atom.coordinates.x
            y_sum = y_sum + atom.coordinates.y
            z_sum = z_sum + atom.coordinates.z

        if total == 0:
            return  # to prevent errors in some cases

        center = Point(x_sum / total, y_sum / total, z_sum / total)

        # now get the radius of the aromatic ring
        radius = 0.0
        for index in indicies_of_ring:
            atom = self.all_atoms[index]
            dist = center.dist_to(atom.coordinates)
            if dist > radius:
                radius = dist

        # now get the plane that defines this ring
        if len(indicies_of_ring) < 3:
            # to prevent an error in some cases. If there aren't three point,
            # you can't define a plane
            return
        elif len(indicies_of_ring) == 3:
            A = self.all_atoms[indicies_of_ring[0]].coordinates
            B = self.all_atoms[indicies_of_ring[1]].coordinates
            C = self.all_atoms[indicies_of_ring[2]].coordinates
        elif len(indicies_of_ring) == 4:
            A = self.all_atoms[indicies_of_ring[0]].coordinates
            B = self.all_atoms[indicies_of_ring[1]].coordinates
            C = self.all_atoms[indicies_of_ring[3]].coordinates
        else:  # best, for 5 and 6 member rings
            A = self.all_atoms[indicies_of_ring[0]].coordinates
            B = self.all_atoms[indicies_of_ring[2]].coordinates
            C = self.all_atoms[indicies_of_ring[4]].coordinates

        AB = self.functions.vector_subtraction(B, A)
        AC = self.functions.vector_subtraction(C, A)
        ABXAC = self.functions.cross_product(AB, AC)

        # formula for plane will be ax + by + cz = d
        x1 = self.all_atoms[indicies_of_ring[0]].coordinates.x
        y1 = self.all_atoms[indicies_of_ring[0]].coordinates.y
        z1 = self.all_atoms[indicies_of_ring[0]].coordinates.z

        a = ABXAC.x
        b = ABXAC.y
        c = ABXAC.z
        d = a * x1 + b * y1 + c * z1

        ar_ring = self.AromaticRing(center, indicies_of_ring, [a, b, c, d], radius)
        self.aromatic_rings.append(ar_ring)

    """
    Class AromaticRing defines an a ring
    """

    class AromaticRing:
        # Initialize a ring
        # Param center (float): center of ring on grid
        # Param indices (list): indecies of atoms in ring
        # Param plane_coeff ():
        # Param radius (float): ring's radius
        def __init__(self, center, indices, plane_coeff, radius):
            self.center = center
            self.indices = indices
            self.plane_coeff = plane_coeff  # a*x + b*y + c*z = dI think that
            self.radius = radius

    # Denote aromatic rings
    def assign_aromatic_rings(self):
        # Get all the rings containing each of the atoms in the ligand
        all_rings = []
        for atom_index in self.non_protein_atoms:
            all_rings.extend(self.all_rings_containing_atom(atom_index))

        for ring_index_1 in range(len(all_rings)):
            ring1 = all_rings[ring_index_1]
            if len(ring1) != 0:
                for ring_index_2 in range(len(all_rings)):
                    if ring_index_1 != ring_index_2:
                        ring2 = all_rings[ring_index_2]
                        if len(ring2) != 0:
                            if self.set1_is_subset_of_set2(ring1, ring2) == True:
                                all_rings[ring_index_2] = []

        while [] in all_rings:
            all_rings.remove([])

        # Now we need to figure out which of these ligands are aromatic (planar)

        for ring_index in range(len(all_rings)):
            ring = all_rings[ring_index]
            is_flat = True
            for t in range(-3, len(ring) - 3):
                pt1 = self.non_protein_atoms[ring[t]].coordinates
                pt2 = self.non_protein_atoms[ring[t + 1]].coordinates
                pt3 = self.non_protein_atoms[ring[t + 2]].coordinates
                pt4 = self.non_protein_atoms[ring[t + 3]].coordinates

                # first, let's see if the last atom in this ring is a carbon connected to four atoms. That would be a quick way of telling this is not an aromatic ring
                cur_atom = self.non_protein_atoms[ring[t + 3]]
                if cur_atom.element == "C" and cur_atom.number_of_neighbors() == 4:
                    is_flat = False
                    break

                # now check the dihedral between the ring atoms to see if it's
                # flat
                angle = self.functions.dihedral(pt1, pt2, pt3, pt4) * 180 / math.pi
                if (angle > -165 and angle < -15) or (angle > 15 and angle < 165):
                    # 15 degress is the cutoff #, ring[t], ring[t+1],
                    # ring[t+2], ring[t+3] # range of this function is -pi to
                    # pi
                    is_flat = False
                    break

                # now check the dihedral between the ring atoms and an atom connected to the current atom to see if that's flat too.
                for substituent_atom_index in cur_atom.indecies_of_atoms_connecting:
                    pt_sub = self.non_protein_atoms[substituent_atom_index].coordinates
                    angle = (
                        self.functions.dihedral(pt2, pt3, pt4, pt_sub) * 180 / math.pi
                    )
                    if (angle > -165 and angle < -15) or (
                        angle > 15 and angle < 165
                    ):  # 15 degress is the cutoff #, ring[t], ring[t+1], ring[t+2], ring[t+3] # range of this function is -pi to pi
                        is_flat = False
                        break

            if is_flat == False:
                all_rings[ring_index] = []
            if len(ring) < 5:
                all_rings[
                    ring_index
                ] = []  # While I'm at it, three and four member rings are not aromatic
            if len(ring) > 6:
                all_rings[
                    ring_index
                ] = (
                    []
                )  # While I'm at it, if the ring has more than 6, also throw it out. So only 5 and 6 member rings are allowed.

        while [] in all_rings:
            all_rings.remove([])

        for ring in all_rings:
            self.add_aromatic_marker(ring)

        # Now that you've found all the rings in non-protein residues, it's
        # time to look for aromatic rings in protein residues
        curr_res = ""
        first = True
        residue = []

        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]

            key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

            if first == True:
                curr_res = key
                first = False

            if key != curr_res:

                self.assign_aromatic_rings_from_protein_process_residue(
                    residue, last_key
                )

                residue = []
                curr_res = key

            residue.append(atom_index)
            last_key = key

        self.assign_aromatic_rings_from_protein_process_residue(residue, last_key)

    # Denote aromatic ring with residue
    # Param residue ():
    # Param last_key ():
    def assign_aromatic_rings_from_protein_process_residue(self, residue, last_key):
        temp = last_key.strip().split("_")
        resname = temp[0]
        real_resname = resname[-3:]
        resid = temp[1]
        chain = temp[2]

        if real_resname == "PHE":
            indicies_of_ring = []

            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CG":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CZ":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD2":
                    indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

        if real_resname == "TYR":
            indicies_of_ring = []

            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CG":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CZ":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD2":
                    indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

        if (
            real_resname == "HIS"
            or real_resname == "HID"
            or real_resname == "HIE"
            or real_resname == "HIP"
        ):
            indicies_of_ring = []

            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CG":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "ND1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "NE2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD2":
                    indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

        if real_resname == "TRP":
            indicies_of_ring = []

            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CG":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "NE1":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD2":
                    indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

            indicies_of_ring = []

            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CD2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CE3":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CZ3":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CH2":
                    indicies_of_ring.append(index)
            for index in residue:  # written this way because order is important
                atom = self.all_atoms[index]
                if atom.atom_name.strip() == "CZ2":
                    indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

    # Return true if set one is a subset of set 2
    # Param set1 (set)
    # Param set2 (set)
    def set1_is_subset_of_set2(self, set1, set2):
        is_subset = True
        for item in set1:
            if item not in set2:
                is_subset = False
                break
        return is_subset

    # Return indecies of atoms in ring
    def all_rings_containing_atom(self, index):

        all_rings = []

        atom = self.all_atoms[index]
        for conneceted_atom in atom.indecies_of_atoms_connecting:
            self.ring_recursive(conneceted_atom, [index], index, all_rings)

        return all_rings

    # TODO not sure what this function does
    def ring_recursive(self, index, already_crossed, orig_atom, all_rings):

        if len(already_crossed) > 6:
            return  # since you're only considering aromatic rings containing 5 or 6 members anyway, save yourself some time.

        atom = self.all_atoms[index]

        temp = already_crossed[:]
        temp.append(index)

        for conneceted_atom in atom.indecies_of_atoms_connecting:
            if not conneceted_atom in already_crossed:
                self.ring_recursive(conneceted_atom, temp, orig_atom, all_rings)
            if conneceted_atom == orig_atom and orig_atom != already_crossed[-1]:
                all_rings.append(temp)

    # Functions to assign secondary structure to protein residues
    # ===========================================================

    def assign_secondary_structure(self):
        # first, we need to know what resid's are available
        resids = []
        last_key = "-99999_Z"
        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]
            key = str(atom.resid) + "_" + atom.chain
            if key != last_key:
                last_key = key
                resids.append(last_key)

        structure = {}
        for resid in resids:
            structure[resid] = "OTHER"

        atoms = []

        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]
            if atom.SideChainOrBackBone() == "BACKBONE":
                if len(atoms) < 8:
                    atoms.append(atom)
                else:
                    atoms.pop(0)
                    atoms.append(atom)

                    # now make sure the first four all have the same resid and
                    # the last four all have the same resid
                    if (
                        atoms[0].resid == atoms[1].resid
                        and atoms[0].resid == atoms[2].resid
                        and atoms[0].resid == atoms[3].resid
                        and atoms[0] != atoms[4].resid
                        and atoms[4].resid == atoms[5].resid
                        and atoms[4].resid == atoms[6].resid
                        and atoms[4].resid == atoms[7].resid
                        and atoms[0].resid + 1 == atoms[7].resid
                        and atoms[0].chain == atoms[7].chain
                    ):
                        resid1 = atoms[0].resid
                        resid2 = atoms[7].resid

                        # Now give easier to use names to the atoms
                        for atom in atoms:
                            if atom.resid == resid1 and atom.atom_name.strip() == "N":
                                first_N = atom
                            if atom.resid == resid1 and atom.atom_name.strip() == "C":
                                first_C = atom
                            if atom.resid == resid1 and atom.atom_name.strip() == "CA":
                                first_CA = atom

                            if atom.resid == resid2 and atom.atom_name.strip() == "N":
                                second_N = atom
                            if atom.resid == resid2 and atom.atom_name.strip() == "C":
                                second_C = atom
                            if atom.resid == resid2 and atom.atom_name.strip() == "CA":
                                second_CA = atom

                        # Now compute the phi and psi dihedral angles
                        phi = (
                            self.functions.dihedral(
                                first_C.coordinates,
                                second_N.coordinates,
                                second_CA.coordinates,
                                second_C.coordinates,
                            )
                            * 180.0
                            / math.pi
                        )
                        psi = (
                            self.functions.dihedral(
                                first_N.coordinates,
                                first_CA.coordinates,
                                first_C.coordinates,
                                second_N.coordinates,
                            )
                            * 180.0
                            / math.pi
                        )

                        # Now use those angles to determine if it's alpha or
                        # beta
                        if phi > -145 and phi < -35 and psi > -70 and psi < 50:
                            key1 = str(first_C.resid) + "_" + first_C.chain
                            key2 = str(second_C.resid) + "_" + second_C.chain
                            structure[key1] = "ALPHA"
                            structure[key2] = "ALPHA"
                        if (phi >= -180 and phi < -40 and psi <= 180 and psi > 90) or (
                            phi >= -180 and phi < -70 and psi <= -165
                        ):
                            # beta. This gets some loops (by my eye), but it's
                            # the best I could do.
                            key1 = str(first_C.resid) + "_" + first_C.chain
                            key2 = str(second_C.resid) + "_" + second_C.chain
                            structure[key1] = "BETA"
                            structure[key2] = "BETA"

        # Now update each of the atoms with this structural information
        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]
            key = str(atom.resid) + "_" + atom.chain
            atom.structure = structure[key]

        # Some more post processing.
        CA_list = []  # first build a list of the indices of all the alpha carbons
        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]
            if (
                atom.residue.strip() in self.protein_resnames
                and atom.atom_name.strip() == "CA"
            ):
                CA_list.append(atom_index)

        # some more post processing.
        change = True
        while change == True:
            change = False

            # A residue of index i is only going to be in an alpha helix its
            # CA is within 6 A of the CA of the residue i + 3
            for CA_atom_index in CA_list:
                CA_atom = self.all_atoms[CA_atom_index]
                if CA_atom.structure == "ALPHA":  # so it's in an alpha helix
                    another_alpha_is_close = False
                    for (
                        other_CA_atom_index
                    ) in CA_list:  # so now compare that CA to all the other CA's
                        other_CA_atom = self.all_atoms[other_CA_atom_index]
                        if other_CA_atom.structure == "ALPHA":
                            # so it's also in an alpha helix
                            if (
                                other_CA_atom.resid - 3 == CA_atom.resid
                                or other_CA_atom.resid + 3 == CA_atom.resid
                            ):  # so this CA atom is one of the ones the first atom might hydrogen bond with
                                if (
                                    other_CA_atom.coordinates.dist_to(
                                        CA_atom.coordinates
                                    )
                                    < 6.0
                                ):  # so these two CA atoms are close enough together that their residues are probably hydrogen bonded
                                    another_alpha_is_close = True
                                    break
                    if another_alpha_is_close == False:
                        self.set_structure_of_residue(
                            CA_atom.chain, CA_atom.resid, "OTHER"
                        )
                        change = True

            # Alpha helices are only alpha helices if they span at least 4
            # residues (to wrap around and hydrogen bond). I'm going to
            # require them to span at least 5 residues, based on examination
            # of many structures.
            for index_in_list in range(len(CA_list) - 5):

                index_in_pdb1 = CA_list[index_in_list]
                index_in_pdb2 = CA_list[index_in_list + 1]
                index_in_pdb3 = CA_list[index_in_list + 2]
                index_in_pdb4 = CA_list[index_in_list + 3]
                index_in_pdb5 = CA_list[index_in_list + 4]
                index_in_pdb6 = CA_list[index_in_list + 5]

                atom1 = self.all_atoms[index_in_pdb1]
                atom2 = self.all_atoms[index_in_pdb2]
                atom3 = self.all_atoms[index_in_pdb3]
                atom4 = self.all_atoms[index_in_pdb4]
                atom5 = self.all_atoms[index_in_pdb5]
                atom6 = self.all_atoms[index_in_pdb6]

                if (
                    atom1.resid + 1 == atom2.resid
                    and atom2.resid + 1 == atom3.resid
                    and atom3.resid + 1 == atom4.resid
                    and atom4.resid + 1 == atom5.resid
                    and atom5.resid + 1 == atom6.resid
                ):  # so they are sequential

                    if (
                        atom1.structure != "ALPHA"
                        and atom2.structure == "ALPHA"
                        and atom3.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        change = True

                    if (
                        atom2.structure != "ALPHA"
                        and atom3.structure == "ALPHA"
                        and atom4.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True

                    if (
                        atom3.structure != "ALPHA"
                        and atom4.structure == "ALPHA"
                        and atom5.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        change = True

                    if (
                        atom4.structure != "ALPHA"
                        and atom5.structure == "ALPHA"
                        and atom6.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

                    if (
                        atom1.structure != "ALPHA"
                        and atom2.structure == "ALPHA"
                        and atom3.structure == "ALPHA"
                        and atom4.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True

                    if (
                        atom2.structure != "ALPHA"
                        and atom3.structure == "ALPHA"
                        and atom4.structure == "ALPHA"
                        and atom5.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        change = True

                    if (
                        atom3.structure != "ALPHA"
                        and atom4.structure == "ALPHA"
                        and atom5.structure == "ALPHA"
                        and atom6.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

                    if (
                        atom1.structure != "ALPHA"
                        and atom2.structure == "ALPHA"
                        and atom3.structure == "ALPHA"
                        and atom4.structure == "ALPHA"
                        and atom5.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        change = True

                    if (
                        atom2.structure != "ALPHA"
                        and atom3.structure == "ALPHA"
                        and atom4.structure == "ALPHA"
                        and atom5.structure == "ALPHA"
                        and atom6.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

                    if (
                        atom1.structure != "ALPHA"
                        and atom2.structure == "ALPHA"
                        and atom3.structure == "ALPHA"
                        and atom4.structure == "ALPHA"
                        and atom5.structure == "ALPHA"
                        and atom6.structure != "ALPHA"
                    ):
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

            # now go through each of the BETA CA atoms. A residue is only
            # going to be called a beta sheet if CA atom is within 6.0 A of
            # another CA beta, same chain, but index difference > 2.
            for CA_atom_index in CA_list:
                CA_atom = self.all_atoms[CA_atom_index]
                if CA_atom.structure == "BETA":  # so it's in a beta sheet
                    another_beta_is_close = False
                    for other_CA_atom_index in CA_list:
                        if other_CA_atom_index != CA_atom_index:
                            # so not comparing an atom to itself
                            other_CA_atom = self.all_atoms[other_CA_atom_index]
                            if other_CA_atom.structure == "BETA":
                                # so you're comparing it only to other
                                # BETA-sheet atoms
                                if other_CA_atom.chain == CA_atom.chain:
                                    # so require them to be on the same chain.
                                    # needed to indecies can be fairly
                                    # compared
                                    if (
                                        math.fabs(other_CA_atom.resid - CA_atom.resid)
                                        > 2
                                    ):
                                        # so the two residues are not simply
                                        # adjacent to each other on the chain
                                        if (
                                            CA_atom.coordinates.dist_to(
                                                other_CA_atom.coordinates
                                            )
                                            < 6.0
                                        ):
                                            # so these to atoms are close to
                                            # each other
                                            another_beta_is_close = True
                                            break

                    if another_beta_is_close == False:
                        self.set_structure_of_residue(
                            CA_atom.chain, CA_atom.resid, "OTHER"
                        )
                        change = True

            # Now some more post-processing needs to be done. Do this again to
            # clear up mess that may have just been created (single residue
            # beta strand, for example) Beta sheets are usually at least 3
            # residues long

            for index_in_list in range(len(CA_list) - 3):

                index_in_pdb1 = CA_list[index_in_list]
                index_in_pdb2 = CA_list[index_in_list + 1]
                index_in_pdb3 = CA_list[index_in_list + 2]
                index_in_pdb4 = CA_list[index_in_list + 3]

                atom1 = self.all_atoms[index_in_pdb1]
                atom2 = self.all_atoms[index_in_pdb2]
                atom3 = self.all_atoms[index_in_pdb3]
                atom4 = self.all_atoms[index_in_pdb4]

                if (
                    atom1.resid + 1 == atom2.resid
                    and atom2.resid + 1 == atom3.resid
                    and atom3.resid + 1 == atom4.resid
                ):  # so they are sequential

                    if (
                        atom1.structure != "BETA"
                        and atom2.structure == "BETA"
                        and atom3.structure != "BETA"
                    ):
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        change = True

                    if (
                        atom2.structure != "BETA"
                        and atom3.structure == "BETA"
                        and atom4.structure != "BETA"
                    ):
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True

                    if (
                        atom1.structure != "BETA"
                        and atom2.structure == "BETA"
                        and atom3.structure == "BETA"
                        and atom4.structure != "BETA"
                    ):
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True

    def set_structure_of_residue(self, chain, resid, structure):
        for atom_index in self.all_atoms:
            atom = self.all_atoms[atom_index]
            if atom.chain == chain and atom.resid == resid:
                atom.structure = structure


"""
Class MathFunctions
"""


class MathFunctions:
    def planrity(self, point1, point2, point3, point4):

        x1 = point1.x
        y1 = point1.y
        z1 = point1.z
        x2 = point2.x
        y2 = point2.y
        z2 = point2.z
        x3 = point3.x
        y3 = point3.y
        z3 = point3.z
        x4 = point4.x
        y4 = point4.y
        z4 = point4.z

        A = (y1 * (z2 - z3)) + (y2 * (z3 - z1)) + (y3 * (z1 - z2))
        B = (z1 * (x2 - x3)) + (z2 * (x3 - x1)) + (z3 * (x1 - x2))
        C = (x1 * (y2 - y3)) + (x2 * (y3 - y1)) + (x3 * (y1 - y2))
        D = (
            ((-x1) * ((y2 * z3) - (y3 * z2)))
            + ((-x2) * ((y3 * z1) - (y1 * z3)))
            + ((-x3) * ((y1 * z2) - (y2 * z1)))
        )
        distance = (math.fabs((A * x4) + (B * y4) + (C * z4) + D)) / (
            math.sqrt(math.pow(A, 2) + math.pow(B, 2) + math.pow(C, 2))
        )

        A1 = (y1 * (z2 - z4)) + (y2 * (z4 - z1)) + (y4 * (z1 - z2))
        B1 = (z1 * (x2 - x4)) + (z2 * (x4 - x1)) + (z4 * (x1 - x2))
        C1 = (x1 * (y2 - y4)) + (x2 * (y4 - y1)) + (x4 * (y1 - y2))
        D1 = (
            ((-x1) * ((y2 * z4) - (y4 * z2)))
            + ((-x2) * ((y4 * z1) - (y1 * z4)))
            + ((-x4) * ((y1 * z2) - (y2 * z1)))
        )
        distance1 = (math.fabs((A1 * x3) + (B1 * y3) + (C1 * z3) + D1)) / (
            math.sqrt(math.pow(A1, 2) + math.pow(B1, 2) + math.pow(C1, 2))
        )

        A2 = (y1 * (z4 - z3)) + (y4 * (z3 - z1)) + (y3 * (z1 - z4))
        B2 = (z1 * (x4 - x3)) + (z4 * (x3 - x1)) + (z3 * (x1 - x4))
        C2 = (x1 * (y4 - y3)) + (x4 * (y3 - y1)) + (x3 * (y1 - y4))
        D2 = (
            ((-x1) * ((y4 * z3) - (y3 * z4)))
            + ((-x4) * ((y3 * z1) - (y1 * z3)))
            + ((-x3) * ((y1 * z4) - (y4 * z1)))
        )
        distance2 = (math.fabs((A2 * x2) + (B2 * y2) + (C2 * z2) + D2)) / (
            math.sqrt(math.pow(A2, 2) + math.pow(B2, 2) + math.pow(C2, 2))
        )

        A3 = (y4 * (z2 - z3)) + (y2 * (z3 - z4)) + (y3 * (z4 - z2))
        B3 = (z4 * (x2 - x3)) + (z2 * (x3 - x4)) + (z3 * (x4 - x2))
        C3 = (x4 * (y2 - y3)) + (x2 * (y3 - y4)) + (x3 * (y4 - y2))
        D3 = (
            ((-x4) * ((y2 * z3) - (y3 * z2)))
            + ((-x2) * ((y3 * z4) - (y4 * z3)))
            + ((-x3) * ((y4 * z2) - (y2 * z4)))
        )
        distance3 = (math.fabs((A3 * x1) + (B3 * y1) + (C3 * z1) + D3)) / (
            math.sqrt(math.pow(A3, 2) + math.pow(B3, 2) + math.pow(C3, 2))
        )

        final_dist = -1

        if distance < distance1 and distance < distance2 and distance < distance3:
            final_dist = distance
        elif distance1 < distance and distance1 < distance2 and distance1 < distance3:
            final_dist = distance1
        elif distance2 < distance and distance2 < distance1 and distance2 < distance3:
            final_dist = distance2
        elif distance3 < distance and distance3 < distance1 and distance3 < distance2:
            final_dist = distance3

        # Now normalize by the length of the longest bond

        return final_dist

    def vector_subtraction(self, vector1, vector2):  # vector1 - vector2
        return Point(
            vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z
        )

    def cross_product(self, Pt1, Pt2):  # never tested
        response = Point(0, 0, 0)

        response.x = Pt1.y * Pt2.z - Pt1.z * Pt2.y
        response.y = Pt1.z * Pt2.x - Pt1.x * Pt2.z
        response.z = Pt1.x * Pt2.y - Pt1.y * Pt2.x

        return response

    def vector_scalar_multiply(self, vector, scalar):
        return Point(vector.x * scalar, vector.y * scalar, vector.z * scalar)

    def dot_product(self, point1, point2):
        return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z

    def dihedral(self, point1, point2, point3, point4):  # never tested

        b1 = self.vector_subtraction(point2, point1)
        b2 = self.vector_subtraction(point3, point2)
        b3 = self.vector_subtraction(point4, point3)

        b2Xb3 = self.cross_product(b2, b3)
        b1Xb2 = self.cross_product(b1, b2)

        b1XMagb2 = self.vector_scalar_multiply(b1, b2.magnitude())
        radians = math.atan2(
            self.dot_product(b1XMagb2, b2Xb3), self.dot_product(b1Xb2, b2Xb3)
        )
        return radians

    def angle_between_three_points(
        self, point1, point2, point3
    ):  # As in three connected atoms
        vector1 = self.vector_subtraction(point1, point2)
        vector2 = self.vector_subtraction(point3, point2)
        return self.angle_between_points(vector1, vector2)

    def angle_between_points(self, point1, point2):
        new_point1 = self.return_normalized_vector(point1)
        new_point2 = self.return_normalized_vector(point2)
        dot_prod = self.dot_product(new_point1, new_point2)
        if dot_prod > 1.0:
            dot_prod = 1.0  # to prevent errors that can rarely occur
        if dot_prod < -1.0:
            dot_prod = -1.0
        return math.acos(dot_prod)

    def return_normalized_vector(self, vector):
        dist = self.distance(Point(0, 0, 0), vector)
        return Point(vector.x / dist, vector.y / dist, vector.z / dist)

    def distance(self, point1, point2):
        deltax = point1.x - point2.x
        deltay = point1.y - point2.y
        deltaz = point1.z - point2.z

        return math.sqrt(
            math.pow(deltax, 2) + math.pow(deltay, 2) + math.pow(deltaz, 2)
        )

    def project_point_onto_plane(self, a_point, plane_coefficients):
        # essentially finds the point on the plane that is closest to the
        # specified point the plane_coefficients are [a,b,c,d], where the
        # plane is ax + by + cz = d

        # First, define a plane using cooeficients a, b, c, d such that ax +
        # by + cz = d
        a = plane_coefficients[0]
        b = plane_coefficients[1]
        c = plane_coefficients[2]
        d = plane_coefficients[3]

        # Now, define a point in space (s,u,v)
        s = a_point.x
        u = a_point.y
        v = a_point.z

        # the formula of a line perpendicular to the plan passing through (s,u,v) is:
        # x = s + at
        # y = u + bt
        # z = v + ct

        t = (d - a * s - b * u - c * v) / (a * a + b * b + c * c)

        # here's the point closest on the plane
        x = s + a * t
        y = u + b * t
        z = v + c * t

        return Point(x, y, z)


"""
Class Binana
"""


class Binana:

    functions = MathFunctions()

    # supporting functions
    def list_alphebetize_and_combine(self, list_obj):
        list_obj.sort()
        return "_".join(list_obj)

    def hashtable_entry_add_one(self, hashtable, key, toadd=1):
        # note that dictionaries (hashtables) are passed by reference in
        # python
        if key in hashtable:
            hashtable[key] = hashtable[key] + toadd
        else:
            hashtable[key] = toadd

    def center(self, string, length):
        while len(string) < length:
            string = " " + string
            if len(string) < length:
                string = string + " "
        return string

    # The meat of the class
    def __init__(self, ligand_pdbqt_filename, receptor_pdbqt_filename):


        parameters = {'close_contacts_dist1_cutoff': 2.5,
                      'close_contacts_dist2_cutoff': 4.0, 
                      'electrostatic_dist_cutoff': 4.0, 
                      'active_site_flexibility_dist_cutoff': 4.0, 
                      'hydrophobic_dist_cutoff': 4.0, 
                      'hydrogen_bond_dist_cutoff': 4.0, 
                      'hydrogen_bond_angle_cutoff': 40.0, 
                      'pi_padding_dist': 0.75, 
                      'pi_pi_interacting_dist_cutoff': 7.5, 
                      'pi_stacking_angle_tolerance': 30.0, 
                      'T_stacking_angle_tolerance': 30.0, 
                      'T_stacking_closest_dist_cutoff': 5.0, 
                      'cation_pi_dist_cutoff': 6.0, 
                      'salt_bridge_dist_cutoff': 5.5,  
                      'output_dir': '', 
                      'output_file': ''}

        ligand = PDB()
        ligand.load_PDB(ligand_pdbqt_filename)

        receptor = PDB()
        receptor.load_PDB(receptor_pdbqt_filename)
        receptor.assign_secondary_structure()

        # Get distance measurements between protein and ligand atom types, as
        # well as some other measurements

        ligand_receptor_atom_type_pairs_less_than_two_half = {}
        ligand_receptor_atom_type_pairs_less_than_four = {}
        ligand_receptor_atom_type_pairs_electrostatic = {}
        active_site_flexibility = {}
        hbonds = {}
        hydrophobics = {}

        functions = MathFunctions()

        pdb_close_contacts = PDB()
        pdb_contacts = PDB()
        pdb_contacts_alpha_helix = PDB()
        pdb_contacts_beta_sheet = PDB()
        pdb_contacts_other_2nd_structure = PDB()
        pdb_side_chain = PDB()
        pdb_back_bone = PDB()
        pdb_hydrophobic = PDB()
        pdb_hbonds = PDB()

        close_contacts_labels = []
        contacts_labels = []
        hydrophobic_labels = []
        hbonds_labels = []

        for ligand_atom_index in ligand.all_atoms:
            for receptor_atom_index in receptor.all_atoms:
                ligand_atom = ligand.all_atoms[ligand_atom_index]
                receptor_atom = receptor.all_atoms[receptor_atom_index]

                dist = ligand_atom.coordinates.dist_to(receptor_atom.coordinates)
                if dist < parameters["close_contacts_dist1_cutoff"]:
                    # less than 2.5 A
                    list_ligand_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
                    self.hashtable_entry_add_one(
                        ligand_receptor_atom_type_pairs_less_than_two_half,
                        self.list_alphebetize_and_combine(list_ligand_atom),
                    )
                    pdb_close_contacts.add_new_atom(ligand_atom.copy_of())
                    pdb_close_contacts.add_new_atom(receptor_atom.copy_of())

                    close_contacts_labels.append(
                        (ligand_atom.string_id(), receptor_atom.string_id())
                    )

                elif dist < parameters["close_contacts_dist2_cutoff"]:
                    # less than 4 A
                    list_ligand_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
                    self.hashtable_entry_add_one(
                        ligand_receptor_atom_type_pairs_less_than_four,
                        self.list_alphebetize_and_combine(list_ligand_atom),
                    )
                    pdb_contacts.add_new_atom(ligand_atom.copy_of())
                    pdb_contacts.add_new_atom(receptor_atom.copy_of())

                    contacts_labels.append(
                        (ligand_atom.string_id(), receptor_atom.string_id())
                    )

                if dist < parameters["electrostatic_dist_cutoff"]:
                    # calculate electrostatic energies for all less than 4 A
                    ligand_charge = ligand_atom.charge
                    receptor_charge = receptor_atom.charge
                    coulomb_energy = (
                        ligand_charge * receptor_charge / dist
                    ) * 138.94238460104697e4  # to convert into J/mol # might be nice to double check this
                    list_ligand_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
                    self.hashtable_entry_add_one(
                        ligand_receptor_atom_type_pairs_electrostatic,
                        self.list_alphebetize_and_combine(list_ligand_atom),
                        coulomb_energy,
                    )

                if dist < parameters["active_site_flexibility_dist_cutoff"]:
                    # Now get statistics to judge active-site flexibility
                    flexibility_key = (
                        receptor_atom.SideChainOrBackBone()
                        + "_"
                        + receptor_atom.structure
                    )  # first can be sidechain or backbone, second back be alpha, beta, or other, so six catagories
                    if receptor_atom.structure == "ALPHA":
                        pdb_contacts_alpha_helix.add_new_atom(receptor_atom.copy_of())
                    elif receptor_atom.structure == "BETA":
                        pdb_contacts_beta_sheet.add_new_atom(receptor_atom.copy_of())
                    elif receptor_atom.structure == "OTHER":
                        pdb_contacts_other_2nd_structure.add_new_atom(
                            receptor_atom.copy_of()
                        )

                    if receptor_atom.SideChainOrBackBone() == "BACKBONE":
                        pdb_back_bone.add_new_atom(receptor_atom.copy_of())
                    elif receptor_atom.SideChainOrBackBone() == "SIDECHAIN":
                        pdb_side_chain.add_new_atom(receptor_atom.copy_of())

                    self.hashtable_entry_add_one(
                        active_site_flexibility, flexibility_key
                    )

                if dist < parameters["hydrophobic_dist_cutoff"]:
                    # Now see if there's hydrophobic contacts (C-C contacts)
                    if ligand_atom.element == "C" and receptor_atom.element == "C":
                        hydrophobic_key = (
                            receptor_atom.SideChainOrBackBone()
                            + "_"
                            + receptor_atom.structure
                        )
                        pdb_hydrophobic.add_new_atom(ligand_atom.copy_of())
                        pdb_hydrophobic.add_new_atom(receptor_atom.copy_of())

                        self.hashtable_entry_add_one(hydrophobics, hydrophobic_key)

                        hydrophobic_labels.append(
                            (ligand_atom.string_id(), receptor_atom.string_id())
                        )

                if dist < parameters["hydrogen_bond_dist_cutoff"]:
                    # Now see if there's some sort of hydrogen bond between
                    # these two atoms. distance cutoff = 4, angle cutoff = 40.
                    # Note that this is liberal.
                    if (ligand_atom.element == "O" or ligand_atom.element == "N") and (
                        receptor_atom.element == "O" or receptor_atom.element == "N"
                    ):

                        # now build a list of all the hydrogens close to these
                        # atoms
                        hydrogens = []

                        for atm_index in ligand.all_atoms:
                            if ligand.all_atoms[atm_index].element == "H":
                                # so it's a hydrogen
                                if (
                                    ligand.all_atoms[atm_index].coordinates.dist_to(
                                        ligand_atom.coordinates
                                    )
                                    < 1.3
                                ):
                                    # O-H distance is 0.96 A, N-H is 1.01 A.
                                    # See
                                    # http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
                                    ligand.all_atoms[atm_index].comment = "LIGAND"
                                    hydrogens.append(ligand.all_atoms[atm_index])

                        for atm_index in receptor.all_atoms:
                            if receptor.all_atoms[atm_index].element == "H":
                                # so it's a hydrogen
                                if (
                                    receptor.all_atoms[atm_index].coordinates.dist_to(
                                        receptor_atom.coordinates
                                    )
                                    < 1.3
                                ):
                                    # O-H distance is 0.96 A, N-H is 1.01 A.
                                    # See
                                    # http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
                                    receptor.all_atoms[atm_index].comment = "RECEPTOR"
                                    hydrogens.append(receptor.all_atoms[atm_index])

                        # now we need to check the angles
                        for hydrogen in hydrogens:
                            if (
                                math.fabs(
                                    180
                                    - functions.angle_between_three_points(
                                        ligand_atom.coordinates,
                                        hydrogen.coordinates,
                                        receptor_atom.coordinates,
                                    )
                                    * 180.0
                                    / math.pi
                                )
                                <= parameters["hydrogen_bond_angle_cutoff"]
                            ):
                                hbonds_key = (
                                    "HDONOR_"
                                    + hydrogen.comment
                                    + "_"
                                    + receptor_atom.SideChainOrBackBone()
                                    + "_"
                                    + receptor_atom.structure
                                )
                                pdb_hbonds.add_new_atom(ligand_atom.copy_of())
                                pdb_hbonds.add_new_atom(hydrogen.copy_of())
                                pdb_hbonds.add_new_atom(receptor_atom.copy_of())
                                self.hashtable_entry_add_one(hbonds, hbonds_key)

                                hbonds_labels.append(
                                    (
                                        ligand_atom.string_id(),
                                        hydrogen.string_id(),
                                        receptor_atom.string_id(),
                                    )
                                )

        # Get the total number of each atom type in the ligand
        ligand_atom_types = {}
        for ligand_atom_index in ligand.all_atoms:
            ligand_atom = ligand.all_atoms[ligand_atom_index]
            self.hashtable_entry_add_one(ligand_atom_types, ligand_atom.atom_type)

        # This is perhaps controversial. I noticed that often a pi-cation
        # interaction or other pi interaction was only slightly off, but
        # looking at the structure, it was clearly supposed to be a pi-cation
        # interaction. I've decided then to artificially expand the radius of
        # each pi ring. Think of this as adding in a VDW radius, or accounting
        # for poor crystal-structure resolution, or whatever you want to
        # justify it.
        pi_padding = parameters[
            "pi_padding_dist"
        ]

        # Count pi-pi stacking and pi-T stacking interactions
        PI_interactions = {}
        pdb_pistack = PDB()
        pdb_pi_T = PDB()
        pi_stacking_labels = []
        T_stacking_labels = []

        # "PI-Stacking Interactions ALIVE AND WELL IN PROTEINS" says distance
        # of 7.5 A is good cutoff. This seems really big to me, except that
        # pi-pi interactions (parallel) are actuall usually off centered.
        # Interesting paper. Note that adenine and tryptophan count as two
        # aromatic rings. So, for example, an interaction between these two,
        # if positioned correctly, could count for 4 pi-pi interactions.
        for aromatic1 in ligand.aromatic_rings:
            for aromatic2 in receptor.aromatic_rings:
                dist = aromatic1.center.dist_to(aromatic2.center)
                if dist < parameters["pi_pi_interacting_dist_cutoff"]:
                    # so there could be some pi-pi interactions. first, let's
                    # check for stacking interactions. Are the two pi's
                    # roughly parallel?
                    aromatic1_norm_vector = Point(
                        aromatic1.plane_coeff[0],
                        aromatic1.plane_coeff[1],
                        aromatic1.plane_coeff[2],
                    )

                    aromatic2_norm_vector = Point(
                        aromatic2.plane_coeff[0],
                        aromatic2.plane_coeff[1],
                        aromatic2.plane_coeff[2],
                    )

                    angle_between_planes = (
                        self.functions.angle_between_points(
                            aromatic1_norm_vector, aromatic2_norm_vector
                        )
                        * 180.0
                        / math.pi
                    )

                    if (
                        math.fabs(angle_between_planes - 0)
                        < parameters["pi_stacking_angle_tolerance"]
                        or math.fabs(angle_between_planes - 180)
                        < parameters["pi_stacking_angle_tolerance"]
                    ):
                        # so they're more or less parallel, it's probably
                        # pi-pi stackingoutput_dir now, pi-pi are not usually
                        # right on top of each other. They're often staggared.
                        # So I don't want to just look at the centers of the
                        # rings and compare. Let's look at each of the atoms.
                        # do atom of the atoms of one ring, when projected
                        # onto the plane of the other, fall within that other
                        # ring?

                        # start by assuming it's not a pi-pi stacking interaction
                        pi_pi = False

                        for ligand_ring_index in aromatic1.indices:
                            # project the ligand atom onto the plane of the
                            # receptor ring
                            pt_on_receptor_plane = self.functions.project_point_onto_plane(
                                ligand.all_atoms[ligand_ring_index].coordinates,
                                aromatic2.plane_coeff,
                            )
                            if (
                                pt_on_receptor_plane.dist_to(aromatic2.center)
                                <= aromatic2.radius + pi_padding
                            ):
                                pi_pi = True
                                break

                        if pi_pi == False:
                            # if you've already determined it's a pi-pi
                            # stacking interaction, no need to keep trying
                            for receptor_ring_index in aromatic2.indices:
                                # project the ligand atom onto the plane of the receptor ring
                                pt_on_ligand_plane = self.functions.project_point_onto_plane(
                                    receptor.all_atoms[receptor_ring_index].coordinates,
                                    aromatic1.plane_coeff,
                                )
                                if (
                                    pt_on_ligand_plane.dist_to(aromatic1.center)
                                    <= aromatic1.radius + pi_padding
                                ):
                                    pi_pi = True
                                    break

                        if pi_pi == True:
                            structure = receptor.all_atoms[
                                aromatic2.indices[0]
                            ].structure
                            if structure == "":
                                # since it could be interacting with a
                                # cofactor or something
                                structure = "OTHER"
                            key = "STACKING_" + structure

                            for index in aromatic1.indices:
                                pdb_pistack.add_new_atom(
                                    ligand.all_atoms[index].copy_of()
                                )
                            for index in aromatic2.indices:
                                pdb_pistack.add_new_atom(
                                    receptor.all_atoms[index].copy_of()
                                )

                            self.hashtable_entry_add_one(PI_interactions, key)

                            pi_stacking_labels.append(
                                (
                                    "["
                                    + " / ".join(
                                        [
                                            ligand.all_atoms[index].string_id()
                                            for index in aromatic1.indices
                                        ]
                                    )
                                    + "]",
                                    "["
                                    + " / ".join(
                                        [
                                            receptor.all_atoms[index].string_id()
                                            for index in aromatic2.indices
                                        ]
                                    )
                                    + "]",
                                )
                            )

                    elif (
                        math.fabs(angle_between_planes - 90)
                        < parameters["T_stacking_angle_tolerance"]
                        or math.fabs(angle_between_planes - 270)
                        < parameters["T_stacking_angle_tolerance"]
                    ):
                        # so they're more or less perpendicular, it's probably
                        # a pi-edge interaction

                        # having looked at many structures, I noticed the
                        # algorithm was identifying T-pi reactions when the
                        # two rings were in fact quite distant, often with
                        # other atoms in between. Eye-balling it, requiring
                        # that at their closest they be at least 5 A apart
                        # seems to separate the good T's from the bad
                        min_dist = 100.0
                        for ligand_ind in aromatic1.indices:
                            ligand_at = ligand.all_atoms[ligand_ind]
                            for receptor_ind in aromatic2.indices:
                                receptor_at = receptor.all_atoms[receptor_ind]
                                dist = ligand_at.coordinates.dist_to(
                                    receptor_at.coordinates
                                )
                                if dist < min_dist:
                                    min_dist = dist

                        if (
                            min_dist
                            <= parameters["T_stacking_closest_dist_cutoff"]
                        ):
                            # so at their closest points, the two rings come
                            # within 5 A of each other.

                            # okay, is the ligand pi pointing into the
                            # receptor pi, or the other way around? first,
                            # project the center of the ligand pi onto the
                            # plane of the receptor pi, and vs. versa

                            # This could be directional somehow, like a
                            # hydrogen bond.

                            pt_on_receptor_plane = self.functions.project_point_onto_plane(
                                aromatic1.center, aromatic2.plane_coeff
                            )
                            pt_on_lignad_plane = self.functions.project_point_onto_plane(
                                aromatic2.center, aromatic1.plane_coeff
                            )

                            # now, if it's a true pi-T interaction, this
                            # projected point should fall within the ring
                            # whose plane it's been projected into.
                            if (
                                pt_on_receptor_plane.dist_to(aromatic2.center)
                                <= aromatic2.radius + pi_padding
                            ) or (
                                pt_on_lignad_plane.dist_to(aromatic1.center)
                                <= aromatic1.radius + pi_padding
                            ):
                                # so it is in the ring on the projected plane.
                                structure = receptor.all_atoms[
                                    aromatic2.indices[0]
                                ].structure
                                if structure == "":
                                    # since it could be interacting with a
                                    # cofactor or something
                                    structure = "OTHER"

                                key = "T-SHAPED_" + structure

                                for index in aromatic1.indices:
                                    pdb_pi_T.add_new_atom(
                                        ligand.all_atoms[index].copy_of()
                                    )
                                for index in aromatic2.indices:
                                    pdb_pi_T.add_new_atom(
                                        receptor.all_atoms[index].copy_of()
                                    )

                                self.hashtable_entry_add_one(PI_interactions, key)

                                T_stacking_labels.append(
                                    (
                                        "["
                                        + " / ".join(
                                            [
                                                ligand.all_atoms[index].string_id()
                                                for index in aromatic1.indices
                                            ]
                                        )
                                        + "]",
                                        "["
                                        + " / ".join(
                                            [
                                                receptor.all_atoms[index].string_id()
                                                for index in aromatic2.indices
                                            ]
                                        )
                                        + "]",
                                    )
                                )

        # Now identify pi-cation interactions
        pdb_pi_cat = PDB()
        pi_cat_labels = []

        for aromatic in receptor.aromatic_rings:
            for charged in ligand.charges:
                if charged.positive == True:
                    # so only consider positive charges
                    if (
                        charged.coordinates.dist_to(aromatic.center)
                        < parameters["cation_pi_dist_cutoff"]
                    ):
                        # distance cutoff based on "Cation-pi interactions in
                        # structural biology." project the charged onto the
                        # plane of the aromatic
                        charge_projected = self.functions.project_point_onto_plane(
                            charged.coordinates, aromatic.plane_coeff
                        )

                        if (
                            charge_projected.dist_to(aromatic.center)
                            < aromatic.radius + pi_padding
                        ):
                            structure = receptor.all_atoms[
                                aromatic.indices[0]
                            ].structure
                            if structure == "":
                                # since it could be interacting with a
                                # cofactor or something
                                structure = "OTHER"

                            key = "PI-CATION_LIGAND-CHARGED_" + structure

                            for index in aromatic.indices:
                                pdb_pi_cat.add_new_atom(
                                    receptor.all_atoms[index].copy_of()
                                )
                            for index in charged.indices:
                                pdb_pi_cat.add_new_atom(
                                    ligand.all_atoms[index].copy_of()
                                )

                            self.hashtable_entry_add_one(PI_interactions, key)

                            pi_cat_labels.append(
                                (
                                    "["
                                    + " / ".join(
                                        [
                                            ligand.all_atoms[index].string_id()
                                            for index in charged.indices
                                        ]
                                    )
                                    + "]",
                                    "["
                                    + " / ".join(
                                        [
                                            receptor.all_atoms[index].string_id()
                                            for index in aromatic.indices
                                        ]
                                    )
                                    + "]",
                                )
                            )

        for aromatic in ligand.aromatic_rings:
            # now it's the ligand that has the aromatic group
            for charged in receptor.charges:
                if charged.positive == True:
                    # so only consider positive charges
                    if (
                        charged.coordinates.dist_to(aromatic.center)
                        < parameters["cation_pi_dist_cutoff"]
                    ):
                        # distance cutoff based on "Cation-pi interactions in
                        # structural biology." project the charged onto the
                        # plane of the aromatic
                        charge_projected = self.functions.project_point_onto_plane(
                            charged.coordinates, aromatic.plane_coeff
                        )

                        if (
                            charge_projected.dist_to(aromatic.center)
                            < aromatic.radius + pi_padding
                        ):
                            structure = receptor.all_atoms[charged.indices[0]].structure
                            if structure == "":
                                # since it could be interacting with a
                                # cofactor or something
                                structure = "OTHER"

                            key = "PI-CATION_RECEPTOR-CHARGED_" + structure

                            for index in aromatic.indices:
                                pdb_pi_cat.add_new_atom(
                                    ligand.all_atoms[index].copy_of()
                                )
                            for index in charged.indices:
                                pdb_pi_cat.add_new_atom(
                                    receptor.all_atoms[index].copy_of()
                                )

                            self.hashtable_entry_add_one(PI_interactions, key)

                            pi_cat_labels.append(
                                (
                                    "["
                                    + " / ".join(
                                        [
                                            ligand.all_atoms[index].string_id()
                                            for index in aromatic.indices
                                        ]
                                    )
                                    + "]",
                                    "["
                                    + " / ".join(
                                        [
                                            receptor.all_atoms[index].string_id()
                                            for index in charged.indices
                                        ]
                                    )
                                    + "]",
                                )
                            )

        # now count the number of salt bridges
        pdb_salt_bridges = PDB()
        salt_bridges = {}
        salt_bridge_labels = []
        for receptor_charge in receptor.charges:
            for ligand_charge in ligand.charges:
                if ligand_charge.positive != receptor_charge.positive:
                    # so they have oppositve charges
                    if (
                        ligand_charge.coordinates.dist_to(receptor_charge.coordinates)
                        < parameters["salt_bridge_dist_cutoff"]
                    ):
                        # 4  is good cutoff for salt bridges according to
                        # "Close-Range Electrostatic Interactions in
                        # Proteins", but looking at complexes, I decided to go
                        # with 5.5 A
                        structure = receptor.all_atoms[
                            receptor_charge.indices[0]
                        ].structure
                        if structure == "":
                            # since it could be interacting with a cofactor or
                            # something
                            structure = "OTHER"

                        key = "SALT-BRIDGE_" + structure

                        for index in receptor_charge.indices:
                            pdb_salt_bridges.add_new_atom(
                                receptor.all_atoms[index].copy_of()
                            )
                        for index in ligand_charge.indices:
                            pdb_salt_bridges.add_new_atom(
                                ligand.all_atoms[index].copy_of()
                            )

                        self.hashtable_entry_add_one(salt_bridges, key)

                        salt_bridge_labels.append(
                            (
                                "["
                                + " / ".join(
                                    [
                                        ligand.all_atoms[index].string_id()
                                        for index in ligand_charge.indices
                                    ]
                                )
                                + "]",
                                "["
                                + " / ".join(
                                    [
                                        receptor.all_atoms[index].string_id()
                                        for index in receptor_charge.indices
                                    ]
                                )
                                + "]",
                            )
                        )

        # custom output for SCORCH

        # these are the aromatics
        stacking = []
        t_shaped = []
        pi_cation = []

        for key in PI_interactions:
            value = PI_interactions[key]
            together = key + "_" + str(value)
            if "STACKING" in together:
                stacking.append(together)
            if "CATION" in together:
                pi_cation.append(together)
            if "SHAPED" in together:
                t_shaped.append(together)

        # stacking
        stacking_dict = dict()
        for pc in stacking:
            stack_name = '_'.join(pc.split('_')[:-1])
            stacking_dict[stack_name] = pc.split('_')[-1]
        
        # t-shaped
        t_shaped_dict = dict()
        for pc in t_shaped:
            stack_name = '_'.join(pc.split('_')[:-1])
            t_shaped_dict[stack_name] = pc.split('_')[-1]
        
        # pi cations
        pi_cation_dict = dict()
        for pc in pi_cation:
            stack_name = '_'.join(pc.split('_')[:-1])
            pi_cation_dict[stack_name] = pc.split('_')[-1]
        
        # build outputs into nested dictionary
        output = {'closest':ligand_receptor_atom_type_pairs_less_than_two_half,
                  'close':ligand_receptor_atom_type_pairs_less_than_four,
                  'elsums':ligand_receptor_atom_type_pairs_electrostatic,
                  'bpfs':active_site_flexibility,
                  'hbonds':hbonds,
                  'hydrophobics':hydrophobics,
                  'ligand_atoms':ligand_atom_types,
                  'stacking':stacking_dict,
                  't_stacking':t_shaped_dict,
                  'pi_cation':pi_cation_dict,
                  'salt_bridges':salt_bridges,
                  'nrot':ligand.rotateable_bonds_count
                  }

        # assign the output
        self.out = output

        '''
        if parameters["output_dir"] != "":
            # so an output directory has been specified. Write the pdb files
            # out separately

            pdb_close_contacts.save_PDB(
                parameters["output_dir"] + "/close_contacts.pdb"
            )
            pdb_contacts.save_PDB(parameters["output_dir"] + "/contacts.pdb")
            pdb_contacts_alpha_helix.save_PDB(
                parameters["output_dir"] + "/contacts_alpha_helix.pdb"
            )
            pdb_contacts_beta_sheet.save_PDB(
                parameters["output_dir"] + "/contacts_beta_sheet.pdb"
            )
            pdb_contacts_other_2nd_structure.save_PDB(
                parameters["output_dir"]
                + "/contacts_other_secondary_structure.pdb"
            )
            pdb_back_bone.save_PDB(parameters["output_dir"] + "/back_bone.pdb")
            pdb_side_chain.save_PDB(parameters["output_dir"] + "/side_chain.pdb")
            pdb_hydrophobic.save_PDB(
                parameters["output_dir"] + "/hydrophobic.pdb"
            )
            pdb_hbonds.save_PDB(parameters["output_dir"] + "/hydrogen_bonds.pdb")
            pdb_pistack.save_PDB(
                parameters["output_dir"] + "/pi_pi_stacking.pdb"
            )
            pdb_pi_T.save_PDB(parameters["output_dir"] + "/T_stacking.pdb")
            pdb_pi_cat.save_PDB(parameters["output_dir"] + "/cat_pi.pdb")
            pdb_salt_bridges.save_PDB(
                parameters["output_dir"] + "/salt_bridges.pdb"
            )
            ligand.save_PDB(parameters["output_dir"] + "/ligand.pdb")
            receptor.save_PDB(parameters["output_dir"] + "/receptor.pdb")

            f = open(parameters["output_dir"] + "log.txt", "w")
            f.write(output.replace("REMARK ", ""))
            f.close()

            f = open(parameters["output_dir"] + "state.vmd", "w")
            f.write(self.vmd_state_file())
            f.close()
            # TODO: JY put json_file() call here?
        if (
            parameters["output_file"] == ""
            and parameters["output_dir"] == ""
        ):
            # so you're not outputing to either a file or a directory
            print((output.replace("REMARK ", "")))

        if parameters["output_file"] != "":
            # so it's writing to a single file.

            # first, make an explaination.

            explain = (
                'The residue named "CCN" illustrates close contacts where the protein and ligand atoms come within '
                + str(parameters["close_contacts_dist1_cutoff"])
                + ' of each other. "CON" illustrates close contacts where the protein and ligand atoms come within '
                + str(parameters["close_contacts_dist2_cutoff"])
                + ' of each other. "ALP", "BET", and "OTH" illustrates receptor contacts whose respective protein residues have the alpha-helix, beta-sheet, or "other" secondary structure. "BAC" and "SID" illustrate receptor contacts that are part of the protein backbone and sidechain, respectively. "HYD" illustrates hydrophobic contacts between the protein and ligand. "HBN" illustrates hydrogen bonds. "SAL" illustrates salt bridges. "PIS" illustrates pi-pi stacking interactions, "PIT" illustrates T-stacking interactions, and "PIC" illustrates cation-pi interactions. Protein residue names are unchanged, but the ligand residue is now named "LIG".'
            )

            output = output + "REMARK\n"

            lines = textwrap.wrap(explain, 71)
            for line in lines:
                output = output + "REMARK " + line + "\n"

            output = output + "REMARK\n"

            output = (
                output
                + receptor.save_PDB_String()
                + "TER\n"
                + ligand.save_PDB_String()
                + "TER\n"
                + pdb_close_contacts.save_PDB_String()
                + "TER\n"
            )
            output = (
                output
                + pdb_contacts.save_PDB_String()
                + "TER\n"
                + pdb_contacts_alpha_helix.save_PDB_String()
                + "TER\n"
                + pdb_contacts_beta_sheet.save_PDB_String()
                + "TER\n"
            )
            output = (
                output
                + pdb_contacts_other_2nd_structure.save_PDB_String()
                + "TER\n"
                + pdb_back_bone.save_PDB_String()
                + "TER\n"
                + pdb_side_chain.save_PDB_String()
                + "TER\n"
            )
            output = (
                output
                + pdb_hydrophobic.save_PDB_String()
                + "TER\n"
                + pdb_hbonds.save_PDB_String()
                + "TER\n"
                + pdb_pistack.save_PDB_String()
                + "TER\n"
                + pdb_pi_T.save_PDB_String()
                + "TER\n"
            )
            output = (
                output
                + pdb_pi_cat.save_PDB_String()
                + "TER\n"
                + pdb_salt_bridges.save_PDB_String()
                + "TER\n"
            )

            f = open(parameters["output_file"], "w")
            f.write(output)
            f.close()
        '''

    def vmd_state_file(self):
        vmd = []
        vmd.append("set viewplist {}")
        vmd.append("set fixedlist {}")
        vmd.append("# Display settings")
        vmd.append("display projection   Orthographic")
        vmd.append("display depthcue   on")
        vmd.append("display cuestart   0.500000")
        vmd.append("display cueend     10.000000")
        vmd.append("display cuedensity 0.200000")
        vmd.append("display cuemode    Exp2")

        vmd.append(
            "mol new back_bone.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 1.000000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top back_bone.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new side_chain.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 1.000000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top side_chain.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new close_contacts.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 1.000000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top close_contacts.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new contacts.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 0.500000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top contacts.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new contacts_alpha_helix.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 1.000000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top contacts_alpha_helix.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new contacts_beta_sheet.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 1.000000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top contacts_beta_sheet.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new contacts_other_secondary_structure.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 1.000000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top contacts_other_secondary_structure.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new hydrophobic.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation VDW 0.500000 8.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top hydrophobic.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new hydrogen_bonds.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation Licorice 0.300000 10.000000 10.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top hydrogen_bonds.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new salt_bridges.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation Licorice 0.300000 10.000000 10.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top salt_bridges.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new cat_pi.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation Licorice 0.300000 10.000000 10.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top cat_pi.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new pi_pi_stacking.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation Licorice 0.300000 10.000000 10.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top pi_pi_stacking.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new T_stacking.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation Licorice 0.300000 10.000000 10.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top T_stacking.pdb")
        vmd.append("molinfo top set drawn 0")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append(
            "mol new ligand.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation CPK 1.000000 0.300000 8.000000 6.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top ligand.pdb")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")
        vmd.append("set topmol [molinfo top]")

        vmd.append(
            "mol new receptor.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
        )
        vmd.append("mol delrep 0 top")
        vmd.append("mol representation Lines 3.000000")
        vmd.append("mol color Name")
        vmd.append("mol selection {all}")
        vmd.append("mol material Opaque")
        vmd.append("mol addrep top")
        vmd.append("mol selupdate 0 top 0")
        vmd.append("mol colupdate 0 top 0")
        vmd.append("mol scaleminmax top 0 0.000000 0.000000")
        vmd.append("mol smoothrep top 0 0")
        vmd.append("mol drawframes top 0 {now}")
        vmd.append("mol rename top receptor.pdb")
        vmd.append(
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
        )
        vmd.append("lappend viewplist [molinfo top]")

        vmd.append("foreach v $viewplist {")
        vmd.append(
            "  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)"
        )
        vmd.append("}")
        vmd.append("foreach v $fixedlist {")
        vmd.append("  molinfo $v set fixed 1")
        vmd.append("}")
        vmd.append("unset viewplist")
        vmd.append("unset fixedlist")
        vmd.append("mol top $topmol")
        vmd.append("unset topmol")
        vmd.append("color Display {Background} white")
        return "\n".join(vmd)

    '''
    # json output
    def json_file(self):
        json_output = {}
        # display information
        json_output["display"] = [
            "set viewplist {}",
            "set fixedlist {}",
            "# Display settings",
            "display projection   Orthographic",
            "display depthcue   on",
            "display cuestart   0.500000",
            "display cueend     10.000000",
            "display cuedensity 0.200000",
            "display cuemode    Exp2",
        ]

        # new back_bone
        json_output["new back_bone"] = [
            "mol new back_bone.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top back_bone.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new side_chain
        json_output["new side_chain"] = [
            "mol new side_chain.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top side_chain.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new close_contacts
        json_output["new close_contacts"] = [
            "mol new close_contacts.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top close_contacts.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new contacts
        json_output["new contacts"] = [
            "mol new contacts.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top contacts.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new contacts_alpha_helix
        json_output["new contacts_alpha_helix"] = [
            "mol new contacts_alpha_helix.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top contacts_alpha_helix.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new contacts_beta_sheet
        json_output["new contacts_beta_sheet"] = [
            "mol new contacts_beta_sheet.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top contacts_beta_sheet.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new contacts_other_secondary_structure
        json_output["new contacts_other_secondary_structure"] = [
            "mol new contacts_other_secondary_structure.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top contacts_other_secondary_structure.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new hydrophobic
        json_output["new hydrophobic"] = [
            "mol new hydrophobic.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top hydrophobic.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new hydrogen_bonds
        json_output["new hydrogen_bonds"] = [
            "mol new hydrogen_bonds.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top hydrogen_bonds.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new salt_bridges
        json_output["new salt_bridges"] = [
            "mol new salt_bridges.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top salt_bridges.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new cat_pi
        json_output["new cat_pi"] = [
            "mol new cat_pi.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top cat_pi.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new pi_pi_stacking
        json_output["new pi_pi_stacking"] = [
            "mol new pi_pi_stacking.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top pi_pi_stacking.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new T_stacking
        json_output["new T_stacking"] = [
            "mol new T_stacking.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top T_stacking.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new ligand
        json_output["new ligand"] = [
            "mol new ligand.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top ligand.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # new receptor
        json_output["new receptor"] = [
            "mol new receptor.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all",
            "mol delrep 0 top",
            "mol representation VDW 1.000000 8.000000",
            "mol color Name",
            "mol selsection {all}",
            "mol material Opaque",
            "mol addrep top",
            "mol selupdate 0 top 0",
            "mol colupdate 0 to 0",
            "mol scaleminmax top 0 0.000000 0.000000",
            "mol smoothrep top 0 0",
            "mol drawfrarmes top 0 {now}",
            "mol rename top receptor.pdb",
            "molinfo top set drawn 0",
            "set viewpoints([molinfo top]) {{{1 0 0 -75.1819} {0 1 0 -83.0219} {0 0 1 -119.981} {0 0 0 1}} {{-0.0620057 0.672762 -0.737291 0} {0.428709 0.685044 0.589035 0} {0.90135 -0.279568 -0.33089 0} {0 0 0 1}} {{0.11999 0 0 0} {0 0.11999 0 0} {0 0 0.11999 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}"
            "lappend viewplist [molinfo top]",
        ]

        # more display
        json_output = [
            "foreach v $viewplist {",
            "  mol $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)",
            "}",
            "foreach v $fixedlist {",
            "  molinfo $v set fixed 1",
            "}",
            "unset viewplist",
            "unset fixedlist",
            "mol top $topmol",
            "unset topmol",
            "color Display {Background} white",
        ]

        # return output as dictionary
        return json.dumps(json_output)
    '''


class CommandLineParameters:
    params = {}

    def is_num(self, num):
        try:
            t = float(num)
            return t
        except ValueError:
            return num

    def __init__(self, parameters):

        # first, set defaults
        self.params["close_contacts_dist1_cutoff"] = 2.5
        self.params["close_contacts_dist2_cutoff"] = 4.0
        self.params["electrostatic_dist_cutoff"] = 4.0
        self.params["active_site_flexibility_dist_cutoff"] = 4.0
        self.params["hydrophobic_dist_cutoff"] = 4.0
        self.params["hydrogen_bond_dist_cutoff"] = 4.0
        self.params["hydrogen_bond_angle_cutoff"] = 40.0
        self.params["pi_padding_dist"] = 0.75
        self.params["pi_pi_interacting_dist_cutoff"] = 7.5
        self.params["pi_stacking_angle_tolerance"] = 30.0
        self.params["T_stacking_angle_tolerance"] = 30.0
        self.params["T_stacking_closest_dist_cutoff"] = 5.0
        self.params["cation_pi_dist_cutoff"] = 6.0
        self.params["salt_bridge_dist_cutoff"] = 5.5
        self.params["receptor"] = ""
        self.params["ligand"] = ""
        self.params["output_dir"] = ""
        self.params["output_file"] = ""

        # now get user inputed values

        for index in range(len(parameters)):
            item = parameters[index]
            if item[:1] == "-":
                # so it's a parameter key value
                key = item.replace("-", "")
                value = self.is_num(parameters[index + 1])
                if key in list(self.params.keys()):
                    self.params[key] = value
                    parameters[index] = ""
                    parameters[index + 1] = ""

        # make a list of all the command-line parameters not used
        self.error = ""
        for index in range(1, len(parameters)):
            item = parameters[index]
            if item != "":
                self.error = self.error + item + " "

        # Make sure the output directory, if specified, ends in a /
        if self.params["output_dir"] != "":
            if self.params["output_dir"][-1:] != os.sep:
                self.params["output_dir"] = self.params["output_dir"] + os.sep

        # If an output directory is specified but a log file isn't, set a
        # default logfile
        if self.params["output_dir"] != "" and self.params["output_file"] == "":
            self.params["output_file"] = self.params["output_dir"] + "output.pdb"

    def okay_to_proceed(self):
        # at the very least, you need the ligand and the receptor
        if self.params["receptor"] != "" and self.params["ligand"] != "":
            return True
        else:
            return False

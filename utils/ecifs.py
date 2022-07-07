#######################################################################
# This script contains functions from the GitHub referenced in the    #
# published paper: https://doi.org/10.1093/bioinformatics/btaa982     #
# from Difacquim at UNAM - https://www.difacquim.com/. It is intended #
# to replicate their calculations of ECIFs to provide features for    #
# our scoring functions. The csv file "PDB_Atom_Keys.csv" is also the #
# work of Difacquim.                                                  #
#                                                                     #
# Format conversion functions at the end of the script added by       #
# @milesmcgibbon                                                      #
#                                                                     #
#######################################################################

import numpy as np
import pandas as pd
import os
from os import listdir
from rdkit import Chem
from scipy.spatial.distance import cdist
from itertools import product
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
import openbabel as ob

# Possible predefined protein atoms
ECIF_ProteinAtoms = ['C;4;1;3;0;0', 'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1',
                     'C;4;3;0;0;0', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                     'C;5;3;0;0;0', 'C;6;3;0;0;0', 'N;3;1;2;0;0', 'N;3;2;0;1;1',
                     'N;3;2;1;0;0', 'N;3;2;1;1;1', 'N;3;3;0;0;1', 'N;4;1;2;0;0',
                     'N;4;1;3;0;0', 'N;4;2;1;0;0', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                     'S;2;1;1;0;0', 'S;2;2;0;0;0']

# Possible ligand atoms according to the PDBbind 2016 "refined set"
ECIF_LigandAtoms = ['Br;1;1;0;0;0', 'C;3;3;0;1;1', 'C;4;1;1;0;0', 'C;4;1;2;0;0',
                     'C;4;1;3;0;0', 'C;4;2;0;0;0', 'C;4;2;1;0;0', 'C;4;2;1;0;1',
                     'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1', 'C;4;3;0;0;0',
                     'C;4;3;0;0;1', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                     'C;4;4;0;0;0', 'C;4;4;0;0;1', 'C;5;3;0;0;0', 'C;5;3;0;1;1',
                     'C;6;3;0;0;0', 'Cl;1;1;0;0;0', 'F;1;1;0;0;0', 'I;1;1;0;0;0',
                     'N;3;1;0;0;0', 'N;3;1;1;0;0', 'N;3;1;2;0;0', 'N;3;2;0;0;0',
                     'N;3;2;0;0;1', 'N;3;2;0;1;1', 'N;3;2;1;0;0', 'N;3;2;1;0;1',
                     'N;3;2;1;1;1', 'N;3;3;0;0;0', 'N;3;3;0;0;1', 'N;3;3;0;1;1',
                     'N;4;1;2;0;0', 'N;4;1;3;0;0', 'N;4;2;1;0;0', 'N;4;2;2;0;0',
                     'N;4;2;2;0;1', 'N;4;3;0;0;0', 'N;4;3;0;0;1', 'N;4;3;1;0;0',
                     'N;4;3;1;0;1', 'N;4;4;0;0;0', 'N;4;4;0;0;1', 'N;5;2;0;0;0',
                     'N;5;3;0;0;0', 'N;5;3;0;1;1', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                     'O;2;2;0;0;0', 'O;2;2;0;0;1', 'O;2;2;0;1;1', 'P;5;4;0;0;0',
                     'P;6;4;0;0;0', 'P;6;4;0;0;1', 'P;7;4;0;0;0', 'S;2;1;0;0;0',
                     'S;2;1;1;0;0', 'S;2;2;0;0;0', 'S;2;2;0;0;1', 'S;2;2;0;1;1',
                     'S;3;3;0;0;0', 'S;3;3;0;0;1', 'S;4;3;0;0;0', 'S;6;4;0;0;0',
                     'S;6;4;0;0;1', 'S;7;4;0;0;0']

# Exhaustive ECIF list used in SCORCH
exhaustive_ECIFs = ['C;4;1;3;0;0-Br;1;1;0;0;0',
                    'C;4;1;3;0;0-C;4;1;3;0;0',
                    'C;4;1;3;0;0-C;4;2;0;0;0',
                    'C;4;1;3;0;0-C;4;2;1;0;0',
                    'C;4;1;3;0;0-C;4;2;1;0;1',
                    'C;4;1;3;0;0-C;4;2;2;0;0',
                    'C;4;1;3;0;0-C;4;2;2;0;1',
                    'C;4;1;3;0;0-C;4;3;0;0;0',
                    'C;4;1;3;0;0-C;4;3;0;0;1',
                    'C;4;1;3;0;0-C;4;3;1;0;0',
                    'C;4;1;3;0;0-C;4;3;1;0;1',
                    'C;4;1;3;0;0-C;4;4;0;0;0',
                    'C;4;1;3;0;0-C;4;4;0;0;1',
                    'C;4;1;3;0;0-Cl;1;1;0;0;0',
                    'C;4;1;3;0;0-F;1;1;0;0;0',
                    'C;4;1;3;0;0-N;3;1;2;0;0',
                    'C;4;1;3;0;0-N;3;2;1;0;0',
                    'C;4;1;3;0;0-N;3;2;1;0;1',
                    'C;4;1;3;0;0-N;3;3;0;0;0',
                    'C;4;1;3;0;0-N;3;3;0;0;1',
                    'C;4;1;3;0;0-N;5;3;0;0;0',
                    'C;4;1;3;0;0-O;2;1;0;0;0',
                    'C;4;1;3;0;0-O;2;1;1;0;0',
                    'C;4;1;3;0;0-O;2;2;0;0;0',
                    'C;4;1;3;0;0-O;2;2;0;0;1',
                    'C;4;1;3;0;0-P;5;4;0;0;0',
                    'C;4;1;3;0;0-S;2;2;0;0;1',
                    'C;4;1;3;0;0-S;6;4;0;0;0',
                    'C;4;2;1;1;1-C;4;1;3;0;0',
                    'C;4;2;1;1;1-C;4;2;0;0;0',
                    'C;4;2;1;1;1-C;4;2;1;0;0',
                    'C;4;2;1;1;1-C;4;2;1;0;1',
                    'C;4;2;1;1;1-C;4;2;2;0;0',
                    'C;4;2;1;1;1-C;4;2;2;0;1',
                    'C;4;2;1;1;1-C;4;3;0;0;0',
                    'C;4;2;1;1;1-C;4;3;0;0;1',
                    'C;4;2;1;1;1-C;4;3;1;0;0',
                    'C;4;2;1;1;1-C;4;3;1;0;1',
                    'C;4;2;1;1;1-C;4;4;0;0;0',
                    'C;4;2;1;1;1-C;4;4;0;0;1',
                    'C;4;2;1;1;1-Cl;1;1;0;0;0',
                    'C;4;2;1;1;1-F;1;1;0;0;0',
                    'C;4;2;1;1;1-N;3;1;2;0;0',
                    'C;4;2;1;1;1-N;3;2;1;0;0',
                    'C;4;2;1;1;1-N;3;2;1;0;1',
                    'C;4;2;1;1;1-N;3;3;0;0;0',
                    'C;4;2;1;1;1-N;3;3;0;0;1',
                    'C;4;2;1;1;1-N;4;1;3;0;0',
                    'C;4;2;1;1;1-N;4;2;2;0;1',
                    'C;4;2;1;1;1-N;4;4;0;0;0',
                    'C;4;2;1;1;1-N;5;3;0;0;0',
                    'C;4;2;1;1;1-O;2;1;0;0;0',
                    'C;4;2;1;1;1-O;2;1;1;0;0',
                    'C;4;2;1;1;1-O;2;2;0;0;0',
                    'C;4;2;1;1;1-O;2;2;0;0;1',
                    'C;4;2;1;1;1-P;5;4;0;0;0',
                    'C;4;2;1;1;1-S;2;2;0;0;1',
                    'C;4;2;1;1;1-S;6;4;0;0;0',
                    'C;4;2;2;0;0-C;4;1;3;0;0',
                    'C;4;2;2;0;0-C;4;2;0;0;0',
                    'C;4;2;2;0;0-C;4;2;1;0;0',
                    'C;4;2;2;0;0-C;4;2;1;0;1',
                    'C;4;2;2;0;0-C;4;2;2;0;0',
                    'C;4;2;2;0;0-C;4;2;2;0;1',
                    'C;4;2;2;0;0-C;4;3;0;0;0',
                    'C;4;2;2;0;0-C;4;3;0;0;1',
                    'C;4;2;2;0;0-C;4;3;1;0;0',
                    'C;4;2;2;0;0-C;4;3;1;0;1',
                    'C;4;2;2;0;0-C;4;4;0;0;0',
                    'C;4;2;2;0;0-C;4;4;0;0;1',
                    'C;4;2;2;0;0-Cl;1;1;0;0;0',
                    'C;4;2;2;0;0-F;1;1;0;0;0',
                    'C;4;2;2;0;0-N;3;1;0;0;0',
                    'C;4;2;2;0;0-N;3;1;1;0;0',
                    'C;4;2;2;0;0-N;3;1;2;0;0',
                    'C;4;2;2;0;0-N;3;2;1;0;0',
                    'C;4;2;2;0;0-N;3;2;1;0;1',
                    'C;4;2;2;0;0-N;3;3;0;0;1',
                    'C;4;2;2;0;0-N;4;1;3;0;0',
                    'C;4;2;2;0;0-O;2;1;0;0;0',
                    'C;4;2;2;0;0-O;2;1;1;0;0',
                    'C;4;2;2;0;0-O;2;2;0;0;0',
                    'C;4;2;2;0;0-O;2;2;0;0;1',
                    'C;4;2;2;0;0-P;5;4;0;0;0',
                    'C;4;2;2;0;0-S;2;2;0;0;1',
                    'C;4;2;2;0;0-S;6;4;0;0;0',
                    'C;4;2;2;0;1-C;4;1;3;0;0',
                    'C;4;2;2;0;1-C;4;2;1;0;1',
                    'C;4;2;2;0;1-C;4;2;2;0;0',
                    'C;4;2;2;0;1-C;4;2;2;0;1',
                    'C;4;2;2;0;1-C;4;3;0;0;0',
                    'C;4;2;2;0;1-C;4;3;0;0;1',
                    'C;4;2;2;0;1-C;4;3;1;0;1',
                    'C;4;2;2;0;1-N;3;2;1;0;0',
                    'C;4;2;2;0;1-N;3;2;1;0;1',
                    'C;4;2;2;0;1-O;2;1;0;0;0',
                    'C;4;2;2;0;1-O;2;1;1;0;0',
                    'C;4;2;2;0;1-O;2;2;0;0;0',
                    'C;4;2;2;0;1-S;6;4;0;0;0',
                    'C;4;3;0;0;0-C;4;1;3;0;0',
                    'C;4;3;0;0;0-C;4;2;0;0;0',
                    'C;4;3;0;0;0-C;4;2;1;0;0',
                    'C;4;3;0;0;0-C;4;2;1;0;1',
                    'C;4;3;0;0;0-C;4;2;2;0;0',
                    'C;4;3;0;0;0-C;4;2;2;0;1',
                    'C;4;3;0;0;0-C;4;3;0;0;0',
                    'C;4;3;0;0;0-C;4;3;0;0;1',
                    'C;4;3;0;0;0-C;4;3;1;0;0',
                    'C;4;3;0;0;0-C;4;3;1;0;1',
                    'C;4;3;0;0;0-C;4;4;0;0;0',
                    'C;4;3;0;0;0-C;4;4;0;0;1',
                    'C;4;3;0;0;0-Cl;1;1;0;0;0',
                    'C;4;3;0;0;0-F;1;1;0;0;0',
                    'C;4;3;0;0;0-N;3;1;1;0;0',
                    'C;4;3;0;0;0-N;3;1;2;0;0',
                    'C;4;3;0;0;0-N;3;2;1;0;0',
                    'C;4;3;0;0;0-N;3;2;1;0;1',
                    'C;4;3;0;0;0-N;3;3;0;0;0',
                    'C;4;3;0;0;0-N;3;3;0;0;1',
                    'C;4;3;0;0;0-N;4;1;3;0;0',
                    'C;4;3;0;0;0-O;2;1;0;0;0',
                    'C;4;3;0;0;0-O;2;1;1;0;0',
                    'C;4;3;0;0;0-O;2;2;0;0;0',
                    'C;4;3;0;0;0-O;2;2;0;0;1',
                    'C;4;3;0;0;0-P;5;4;0;0;0',
                    'C;4;3;0;0;0-S;2;2;0;0;1',
                    'C;4;3;0;0;0-S;6;4;0;0;0',
                    'C;4;3;0;1;1-C;4;1;3;0;0',
                    'C;4;3;0;1;1-C;4;2;1;0;0',
                    'C;4;3;0;1;1-C;4;2;1;0;1',
                    'C;4;3;0;1;1-C;4;2;2;0;0',
                    'C;4;3;0;1;1-C;4;3;0;0;0',
                    'C;4;3;0;1;1-C;4;3;0;0;1',
                    'C;4;3;0;1;1-C;4;3;1;0;0',
                    'C;4;3;0;1;1-C;4;3;1;0;1',
                    'C;4;3;0;1;1-C;4;4;0;0;0',
                    'C;4;3;0;1;1-C;4;4;0;0;1',
                    'C;4;3;0;1;1-F;1;1;0;0;0',
                    'C;4;3;0;1;1-N;3;1;2;0;0',
                    'C;4;3;0;1;1-N;3;2;1;0;0',
                    'C;4;3;0;1;1-N;3;2;1;0;1',
                    'C;4;3;0;1;1-N;3;3;0;0;1',
                    'C;4;3;0;1;1-N;4;1;3;0;0',
                    'C;4;3;0;1;1-N;4;4;0;0;0',
                    'C;4;3;0;1;1-O;2;1;0;0;0',
                    'C;4;3;0;1;1-O;2;2;0;0;0',
                    'C;4;3;0;1;1-O;2;2;0;0;1',
                    'C;4;3;0;1;1-S;2;2;0;0;1',
                    'C;4;3;1;0;0-C;4;1;3;0;0',
                    'C;4;3;1;0;0-C;4;2;0;0;0',
                    'C;4;3;1;0;0-C;4;2;1;0;0',
                    'C;4;3;1;0;0-C;4;2;1;0;1',
                    'C;4;3;1;0;0-C;4;2;2;0;0',
                    'C;4;3;1;0;0-C;4;2;2;0;1',
                    'C;4;3;1;0;0-C;4;3;0;0;0',
                    'C;4;3;1;0;0-C;4;3;0;0;1',
                    'C;4;3;1;0;0-C;4;3;1;0;0',
                    'C;4;3;1;0;0-C;4;3;1;0;1',
                    'C;4;3;1;0;0-C;4;4;0;0;0',
                    'C;4;3;1;0;0-C;4;4;0;0;1',
                    'C;4;3;1;0;0-Cl;1;1;0;0;0',
                    'C;4;3;1;0;0-F;1;1;0;0;0',
                    'C;4;3;1;0;0-N;3;1;2;0;0',
                    'C;4;3;1;0;0-N;3;2;1;0;0',
                    'C;4;3;1;0;0-N;3;2;1;0;1',
                    'C;4;3;1;0;0-N;3;3;0;0;0',
                    'C;4;3;1;0;0-N;3;3;0;0;1',
                    'C;4;3;1;0;0-N;4;1;3;0;0',
                    'C;4;3;1;0;0-N;5;3;0;0;0',
                    'C;4;3;1;0;0-O;2;1;0;0;0',
                    'C;4;3;1;0;0-O;2;1;1;0;0',
                    'C;4;3;1;0;0-O;2;2;0;0;0',
                    'C;4;3;1;0;0-O;2;2;0;0;1',
                    'C;4;3;1;0;0-P;5;4;0;0;0',
                    'C;4;3;1;0;0-S;2;2;0;0;1',
                    'C;4;3;1;0;0-S;6;4;0;0;0',
                    'C;4;3;1;0;1-C;4;1;3;0;0',
                    'C;4;3;1;0;1-C;4;2;1;0;1',
                    'C;4;3;1;0;1-C;4;3;0;0;1',
                    'C;4;3;1;0;1-O;2;1;0;0;0',
                    'C;5;3;0;0;0-C;4;1;3;0;0',
                    'C;5;3;0;0;0-C;4;2;1;0;1',
                    'C;5;3;0;0;0-C;4;2;2;0;0',
                    'C;5;3;0;0;0-C;4;2;2;0;1',
                    'C;5;3;0;0;0-C;4;3;0;0;0',
                    'C;5;3;0;0;0-C;4;3;0;0;1',
                    'C;5;3;0;0;0-C;4;3;1;0;0',
                    'C;5;3;0;0;0-C;4;3;1;0;1',
                    'C;5;3;0;0;0-C;4;4;0;0;0',
                    'C;5;3;0;0;0-N;3;1;1;0;0',
                    'C;5;3;0;0;0-N;3;1;2;0;0',
                    'C;5;3;0;0;0-N;3;2;1;0;0',
                    'C;5;3;0;0;0-N;3;2;1;0;1',
                    'C;5;3;0;0;0-N;3;3;0;0;1',
                    'C;5;3;0;0;0-N;4;1;3;0;0',
                    'C;5;3;0;0;0-N;4;2;2;0;1',
                    'C;5;3;0;0;0-O;2;1;0;0;0',
                    'C;5;3;0;0;0-O;2;1;1;0;0',
                    'C;5;3;0;0;0-O;2;2;0;0;0',
                    'C;5;3;0;0;0-O;2;2;0;0;1',
                    'C;5;3;0;0;0-P;5;4;0;0;0',
                    'C;6;3;0;0;0-C;4;1;3;0;0',
                    'C;6;3;0;0;0-C;4;2;1;0;1',
                    'C;6;3;0;0;0-C;4;2;2;0;0',
                    'C;6;3;0;0;0-C;4;2;2;0;1',
                    'C;6;3;0;0;0-C;4;3;0;0;0',
                    'C;6;3;0;0;0-C;4;3;0;0;1',
                    'C;6;3;0;0;0-C;4;3;1;0;1',
                    'C;6;3;0;0;0-N;3;2;1;0;1',
                    'C;6;3;0;0;0-O;2;1;0;0;0',
                    'C;6;3;0;0;0-O;2;1;1;0;0',
                    'C;6;3;0;0;0-P;5;4;0;0;0',
                    'N;3;1;2;0;0-C;4;1;3;0;0',
                    'N;3;1;2;0;0-C;4;2;1;0;1',
                    'N;3;1;2;0;0-C;4;2;2;0;0',
                    'N;3;1;2;0;0-C;4;2;2;0;1',
                    'N;3;1;2;0;0-C;4;3;0;0;0',
                    'N;3;1;2;0;0-C;4;3;0;0;1',
                    'N;3;1;2;0;0-C;4;3;1;0;0',
                    'N;3;1;2;0;0-C;4;3;1;0;1',
                    'N;3;1;2;0;0-N;3;2;1;0;0',
                    'N;3;1;2;0;0-N;3;2;1;0;1',
                    'N;3;1;2;0;0-N;3;3;0;0;1',
                    'N;3;1;2;0;0-O;2;1;0;0;0',
                    'N;3;1;2;0;0-O;2;1;1;0;0',
                    'N;3;1;2;0;0-O;2;2;0;0;0',
                    'N;3;2;0;1;1-C;4;1;3;0;0',
                    'N;3;2;0;1;1-C;4;2;1;0;1',
                    'N;3;2;0;1;1-C;4;2;2;0;0',
                    'N;3;2;0;1;1-C;4;2;2;0;1',
                    'N;3;2;0;1;1-C;4;3;0;0;0',
                    'N;3;2;0;1;1-C;4;3;0;0;1',
                    'N;3;2;0;1;1-C;4;3;1;0;1',
                    'N;3;2;0;1;1-N;3;2;1;0;1',
                    'N;3;2;0;1;1-O;2;1;0;0;0',
                    'N;3;2;0;1;1-O;2;1;1;0;0',
                    'N;3;2;1;0;0-C;4;2;1;0;1',
                    'N;3;2;1;0;0-C;4;2;2;0;0',
                    'N;3;2;1;0;0-C;4;2;2;0;1',
                    'N;3;2;1;0;0-C;4;3;0;0;1',
                    'N;3;2;1;0;0-C;4;3;1;0;0',
                    'N;3;2;1;0;0-C;4;3;1;0;1',
                    'N;3;2;1;0;0-C;4;4;0;0;0',
                    'N;3;2;1;0;0-N;3;1;1;0;0',
                    'N;3;2;1;0;0-N;3;1;2;0;0',
                    'N;3;2;1;0;0-N;3;2;1;0;0',
                    'N;3;2;1;0;0-N;3;2;1;0;1',
                    'N;3;2;1;0;0-N;3;3;0;0;1',
                    'N;3;2;1;0;0-N;4;1;3;0;0',
                    'N;3;2;1;0;0-O;2;1;0;0;0',
                    'N;3;2;1;0;0-O;2;1;1;0;0',
                    'N;3;2;1;0;0-O;2;2;0;0;0',
                    'N;3;2;1;0;0-O;2;2;0;0;1',
                    'N;3;2;1;0;0-P;5;4;0;0;0',
                    'N;3;2;1;0;0-S;2;2;0;0;1',
                    'N;3;2;1;0;0-S;6;4;0;0;0',
                    'N;3;2;1;1;1-C;4;1;3;0;0',
                    'N;3;2;1;1;1-C;4;2;1;0;1',
                    'N;3;2;1;1;1-C;4;2;2;0;0',
                    'N;3;2;1;1;1-C;4;2;2;0;1',
                    'N;3;2;1;1;1-C;4;3;0;0;0',
                    'N;3;2;1;1;1-C;4;3;0;0;1',
                    'N;3;2;1;1;1-C;4;3;1;0;0',
                    'N;3;2;1;1;1-C;4;3;1;0;1',
                    'N;3;2;1;1;1-N;3;1;2;0;0',
                    'N;3;2;1;1;1-N;3;2;1;0;0',
                    'N;3;2;1;1;1-N;3;2;1;0;1',
                    'N;3;2;1;1;1-N;3;3;0;0;1',
                    'N;3;2;1;1;1-O;2;1;0;0;0',
                    'N;3;2;1;1;1-O;2;1;1;0;0',
                    'N;3;2;1;1;1-O;2;2;0;0;0',
                    'N;3;2;1;1;1-S;2;1;1;0;0',
                    'N;3;3;0;0;1-C;4;1;3;0;0',
                    'N;3;3;0;0;1-C;4;2;1;0;1',
                    'N;3;3;0;0;1-C;4;3;0;0;1',
                    'N;3;3;0;0;1-N;3;2;1;0;1',
                    'N;3;3;0;0;1-O;2;1;0;0;0',
                    'N;4;1;2;0;0-C;4;1;3;0;0',
                    'N;4;1;2;0;0-C;4;2;1;0;1',
                    'N;4;1;2;0;0-C;4;2;2;0;1',
                    'N;4;1;2;0;0-C;4;3;1;0;0',
                    'N;4;1;2;0;0-O;2;1;1;0;0',
                    'N;4;1;2;0;0-P;5;4;0;0;0',
                    'N;4;1;3;0;0-C;4;2;1;0;1',
                    'N;4;1;3;0;0-C;4;2;2;0;0',
                    'N;4;1;3;0;0-C;4;3;0;0;0',
                    'N;4;1;3;0;0-C;4;3;0;0;1',
                    'N;4;1;3;0;0-C;4;3;1;0;0',
                    'N;4;1;3;0;0-C;4;3;1;0;1',
                    'N;4;1;3;0;0-N;3;2;1;0;1',
                    'N;4;1;3;0;0-O;2;1;0;0;0',
                    'N;4;1;3;0;0-O;2;1;1;0;0',
                    'N;4;1;3;0;0-P;5;4;0;0;0',
                    'N;4;2;1;0;0-C;4;1;3;0;0',
                    'N;4;2;1;0;0-C;4;2;1;0;1',
                    'N;4;2;1;0;0-C;4;2;2;0;0',
                    'N;4;2;1;0;0-C;4;3;0;0;0',
                    'N;4;2;1;0;0-C;4;3;0;0;1',
                    'N;4;2;1;0;0-C;4;3;1;0;1',
                    'N;4;2;1;0;0-O;2;1;0;0;0',
                    'N;4;2;1;0;0-O;2;1;1;0;0',
                    'N;4;2;1;0;0-P;5;4;0;0;0',
                    'O;2;1;0;0;0-C;4;1;3;0;0',
                    'O;2;1;0;0;0-C;4;2;1;0;0',
                    'O;2;1;0;0;0-C;4;2;1;0;1',
                    'O;2;1;0;0;0-C;4;2;2;0;0',
                    'O;2;1;0;0;0-C;4;2;2;0;1',
                    'O;2;1;0;0;0-C;4;3;0;0;0',
                    'O;2;1;0;0;0-C;4;3;0;0;1',
                    'O;2;1;0;0;0-C;4;3;1;0;0',
                    'O;2;1;0;0;0-C;4;3;1;0;1',
                    'O;2;1;0;0;0-C;4;4;0;0;0',
                    'O;2;1;0;0;0-C;4;4;0;0;1',
                    'O;2;1;0;0;0-N;3;1;1;0;0',
                    'O;2;1;0;0;0-N;3;1;2;0;0',
                    'O;2;1;0;0;0-N;3;2;1;0;0',
                    'O;2;1;0;0;0-N;3;2;1;0;1',
                    'O;2;1;0;0;0-N;3;3;0;0;0',
                    'O;2;1;0;0;0-N;3;3;0;0;1',
                    'O;2;1;0;0;0-N;4;1;3;0;0',
                    'O;2;1;0;0;0-O;2;1;0;0;0',
                    'O;2;1;0;0;0-O;2;1;1;0;0',
                    'O;2;1;0;0;0-O;2;2;0;0;0',
                    'O;2;1;0;0;0-O;2;2;0;0;1',
                    'O;2;1;0;0;0-P;5;4;0;0;0',
                    'O;2;1;0;0;0-S;2;2;0;0;1',
                    'O;2;1;0;0;0-S;6;4;0;0;0',
                    'O;2;1;1;0;0-C;4;1;3;0;0',
                    'O;2;1;1;0;0-C;4;2;1;0;1',
                    'O;2;1;1;0;0-C;4;2;2;0;0',
                    'O;2;1;1;0;0-C;4;2;2;0;1',
                    'O;2;1;1;0;0-C;4;3;0;0;0',
                    'O;2;1;1;0;0-C;4;3;0;0;1',
                    'O;2;1;1;0;0-C;4;3;1;0;0',
                    'O;2;1;1;0;0-C;4;3;1;0;1',
                    'O;2;1;1;0;0-C;4;4;0;0;1',
                    'O;2;1;1;0;0-Cl;1;1;0;0;0',
                    'O;2;1;1;0;0-F;1;1;0;0;0',
                    'O;2;1;1;0;0-N;3;1;1;0;0',
                    'O;2;1;1;0;0-N;3;1;2;0;0',
                    'O;2;1;1;0;0-N;3;2;1;0;0',
                    'O;2;1;1;0;0-N;3;2;1;0;1',
                    'O;2;1;1;0;0-N;3;3;0;0;1',
                    'O;2;1;1;0;0-N;4;1;3;0;0',
                    'O;2;1;1;0;0-O;2;1;0;0;0',
                    'O;2;1;1;0;0-O;2;1;1;0;0',
                    'O;2;1;1;0;0-O;2;2;0;0;0',
                    'O;2;1;1;0;0-O;2;2;0;0;1',
                    'O;2;1;1;0;0-P;5;4;0;0;0',
                    'O;2;1;1;0;0-S;6;4;0;0;0',
                    'S;2;1;1;0;0-C;4;1;3;0;0',
                    'S;2;1;1;0;0-C;4;2;1;0;1',
                    'S;2;1;1;0;0-C;4;2;2;0;0',
                    'S;2;1;1;0;0-C;4;2;2;0;1',
                    'S;2;1;1;0;0-C;4;3;0;0;0',
                    'S;2;1;1;0;0-C;4;3;0;0;1',
                    'S;2;1;1;0;0-C;4;3;1;0;1',
                    'S;2;1;1;0;0-N;3;2;1;0;1',
                    'S;2;1;1;0;0-O;2;1;0;0;0',
                    'S;2;1;1;0;0-O;2;1;1;0;0',
                    'S;2;2;0;0;0-C;4;1;3;0;0',
                    'S;2;2;0;0;0-C;4;2;1;0;1',
                    'S;2;2;0;0;0-C;4;2;2;0;0',
                    'S;2;2;0;0;0-C;4;2;2;0;1',
                    'S;2;2;0;0;0-C;4;3;0;0;1',
                    'S;2;2;0;0;0-C;4;3;1;0;0',
                    'S;2;2;0;0;0-C;4;3;1;0;1',
                    'S;2;2;0;0;0-Cl;1;1;0;0;0',
                    'S;2;2;0;0;0-N;3;2;1;0;1',
                    'S;2;2;0;0;0-N;3;3;0;0;1',
                    'S;2;2;0;0;0-O;2;1;0;0;0',
                    'S;2;2;0;0;0-O;2;1;1;0;0']

exhaustive_Ligand_Atoms = [p.split('-')[1] for p in exhaustive_ECIFs]

exhaustive_Protein_Atoms = [p.split('-')[0] for p in exhaustive_ECIFs]

PossibleECIF = [i[0]+"-"+i[1] for i in product(ECIF_ProteinAtoms, ECIF_LigandAtoms)]

ELEMENTS_ProteinAtoms = ["C","N","O", "S"]
ELEMENTS_LigandAtoms = ["Br", "C", "Cl", "F", "I", "N", "O", "P", "S"]
PossibleELEMENTS = [i[0]+"-"+i[1] for i in product(ELEMENTS_ProteinAtoms, ELEMENTS_LigandAtoms)]

LigandDescriptors = ['MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex', 'MinAbsEStateIndex',
                      'qed', 'MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons',
                      'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BalabanJ',
                      'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n',
                      'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Kappa1',
                      'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA14', 'SMR_VSA1', 'SMR_VSA10',
                      'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7',
                      'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12',
                      'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6',
                      'SlogP_VSA7', 'SlogP_VSA8', 'TPSA', 'EState_VSA1', 'EState_VSA10',
                      'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5',
                      'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1',
                      'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5',
                      'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3',
                      'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
                      'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles',
                      'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors',
                      'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
                      'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount', 'MolLogP',
                      'MolMR', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_N',
                      'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO',
                      'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O',
                      'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde',
                      'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide',
                      'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azo', 'fr_barbitur',
                      'fr_benzene', 'fr_bicyclic', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester',
                      'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone',
                      'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone',
                      'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine',
                      'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitroso', 'fr_oxazole',
                      'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond',
                      'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_pyridine', 'fr_quatN',
                      'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole',
                      'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_urea']

DescCalc = MolecularDescriptorCalculator(LigandDescriptors)

def GetAtomType(atom):
# This function takes an atom in a molecule and returns its type as defined for ECIF

    AtomType = [atom.GetSymbol(),
                str(atom.GetExplicitValence()),
                str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() != "H"])),
                str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "H"])),
                str(int(atom.GetIsAromatic())),
                str(int(atom.IsInRing())),
               ]

    return(";".join(AtomType))

def LoadSDFasDF(SDF):
# This function takes an SDF for a ligand as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF

    suppl = Chem.SDMolSupplier()
    suppl.SetData(SDF, sanitize=False)
    m = next(suppl)
    m.UpdatePropertyCache(strict=False)

    ECIF_atoms = []

    for atom in m.GetAtoms():
        if atom.GetSymbol() != "H": # Include only non-hydrogen atoms
            entry = [int(atom.GetIdx())]
            entry.append(GetAtomType(atom))
            pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
            entry.append(float("{0:.4f}".format(pos.x)))
            entry.append(float("{0:.4f}".format(pos.y)))
            entry.append(float("{0:.4f}".format(pos.z)))
            ECIF_atoms.append(entry)

    df = pd.DataFrame(ECIF_atoms)
    df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE","X","Y","Z"]
    if len(set(df["ECIF_ATOM_TYPE"]) - set(ECIF_LigandAtoms)) > 0:
        pass
        # print("WARNING: Ligand contains unsupported atom types. Only supported atom-type pairs are counted.")
    return(df)

Atom_Keys=pd.read_csv(os.path.join(os.path.dirname(__file__),"PDB_Atom_Keys.csv"), sep=",")

def LoadPDBasDF(PDB):
# This function takes a PDB for a protein as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF

    ECIF_atoms = []

    f = PDB.split('\n')
    for i in f:
        if i[:4] == "ATOM":
            # Include only non-hydrogen atoms
            if (len(i[12:16].replace(" ","")) < 4 and i[12:16].replace(" ","")[0] != "H") or (len(i[12:16].replace(" ","")) == 4 and i[12:16].replace(" ","")[1] != "H" and i[12:16].replace(" ","")[0] != "H"):
                ECIF_atoms.append([int(i[6:11]),
                         i[17:20]+"-"+i[12:16].replace(" ",""),
                         float(i[30:38]),
                         float(i[38:46]),
                         float(i[46:54])
                        ])

    df = pd.DataFrame(ECIF_atoms, columns=["ATOM_INDEX","PDB_ATOM","X","Y","Z"])
    df = df.merge(Atom_Keys, left_on='PDB_ATOM', right_on='PDB_ATOM')[["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]].sort_values(by="ATOM_INDEX").reset_index(drop=True)
    if list(df["ECIF_ATOM_TYPE"].isna()).count(True) > 0:
        print("WARNING: Protein contains unsupported atom types. Only supported atom-type pairs are counted.")
    return(df)

def GetPLPairs(PDB_protein, SDF_ligand, distance_cutoff):
# This function returns the protein-ligand atom-type pairs for a given distance cutoff

    # Load both structures as pandas DataFrames
    Target = LoadPDBasDF(PDB_protein)
    Ligand = LoadSDFasDF(SDF_ligand)

    # Take all atoms from the target within a cubic box around the ligand considering the "distance_cutoff criterion"
    for i in ["X","Y","Z"]:
        Target = Target[Target[i] < float(Ligand[i].max())+distance_cutoff]
        Target = Target[Target[i] > float(Ligand[i].min())-distance_cutoff]

    # Remove any atoms we don't use to speed efficiency
    Target = Target.loc[Target.ECIF_ATOM_TYPE.isin(exhaustive_Protein_Atoms)]
    Ligand = Ligand.loc[Ligand.ECIF_ATOM_TYPE.isin(exhaustive_Ligand_Atoms)]

    # Get all possible pairs
    Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
    Pairs = [x[0]+"-"+x[1] for x in Pairs]
    Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])

    # calculate distances
    Distances = cdist(Target[["X","Y","Z"]], Ligand[["X","Y","Z"]], metric="euclidean")
    Distances = Distances.reshape(Distances.shape[0]*Distances.shape[1],1)
    Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

    Pairs = pd.concat([Pairs,Distances], axis=1)
    Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)

    # remove element pair counts as we don't use them
    # Pairs from ELEMENTS could be easily obtained froms pairs from ECIF
    # Pairs["ELEMENTS_PAIR"] = [x.split("-")[0].split(";")[0]+"-"+x.split("-")[1].split(";")[0] for x in Pairs["ECIF_PAIR"]]
    return Pairs

def GetECIF(protein, ligand, distance_cutoff):
# Main function for the calculation of ECIF
    SDF_ligand = SDF(ligand)
    PDB_protein = PDB(protein)
    Pairs = GetPLPairs(PDB_protein, SDF_ligand, distance_cutoff=distance_cutoff)
    ECIF = [list(Pairs["ECIF_PAIR"]).count(x) for x in PossibleECIF]

    return ECIF


def SDF(pdbqt): # uses openbabel to convert pdbqt files to sdf format

    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", "sdf")
    mol = ob.OBMol()
    obConversion.ReadString(mol, pdbqt)

    # add all hydrogens to match atom types as this was done in the publication
    mol.AddHydrogens()

    # make pdbqt block from pose
    outMDL = obConversion.WriteString(mol)

    return outMDL


def PDB(pdbqt): # uses openbabel to convert pdbqt files to pdb format

    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", "pdb")
    mol = ob.OBMol()
    obConversion.ReadFile(mol, pdbqt)

    # make pdbqt block from pose
    outMDL = obConversion.WriteString(mol)

    return outMDL

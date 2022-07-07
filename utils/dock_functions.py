#######################################################################
# This script contains functions used by scoring.py to convert smiles #
# to a dockable pdbqt 3D ligand file and dock these into the given    #
# receptor. These functions are only called when the -dock argument   #
# is supplied to scoring.py                                           #
#                                                                     #
# Script Authors:                                                     #
# @milesmcgibbon                                                      #
#                                                                     #
#######################################################################

import os
from biopandas.pdb import PandasPdb
import logging
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
import sys
import subprocess
import pybel

RDLogger.DisableLog('rdApp.*')

stem_path = os.getcwd()

def get_filepaths(folder_path):
    files = os.listdir(folder_path)
    files = [f'{folder_path}{file}' for file in files]
    return files

def get_smiles(smi_datafile):

    smi_dict = dict()

    lines = open(smi_datafile,'r').read().split('\n')

    lines = [l for l in lines if l != '']

    for count, line in enumerate(lines):
        count = count + 1
        line = line.split()
        if len(line) == 1:
            line = (line[0], count)
        elif line[1] == '':
            line = (line[0], count)
        smi_dict[str(line[1])] = line[0]

    return smi_dict

def clean_smile_with_pybel(smile):

    pybel_mol = pybel.readstring("smi", smile)
    fixed_smile = pybel_mol.write("smi")
    return fixed_smile.strip()


def convert_to_rdmol(file, file_extension, sanitize_bool):

    if file_extension == 'mol':
        ref_lig = Chem.MolFromMolFile(file, sanitize=sanitize_bool, strictParsing=sanitize_bool)
    elif file_extension == 'mol2':
        ref_lig = Chem.MolFromMol2File(file, sanitize=sanitize_bool, cleanupSubstructures=sanitize_bool)
    elif 'sdf' in file_extension:
        ref_lig = Chem.SDMolSupplier(file, sanitize=sanitize_bool)[0]

    return ref_lig


def get_coordinates(file, padding): # find the center x, y, z coordinates of the active crystal ligand and construct a binding site cuboid around this center point

    file_extension = file.split('.')[-1].strip().lstrip()

    if 'pdb' in file_extension:
        pmol = PandasPdb().read_pdb(file)

    else:
        ref_lig = convert_to_rdmol(file, True)

        if ref_lig is None:
            ref_lig = convert_to_rdmol(file, False)

        try:
            pdb_block = Chem.MolToPDBBlock(ref_lig)
        except:
            logging.critical(f"""\nERROR: Could not parse supplied reference ligand '{file}'. Try:\n- supplying coordinates\n- converting this molecule to a pdb file\n- suppling a different reference ligand\n- making sure any mol2 files are Corina produced, not Tripos""")
            sys.exit()

        # read the pdb file as a dataframe
        pmol = PandasPdb()

        pmol.pdb_text = pdb_block

        pmol._df = pmol._construct_df(pdb_lines=pmol.pdb_text.splitlines(True))

    # combine atoms and hetatm coordinates into one master dataframe
    atom_df = pmol.df['ATOM']
    hetatm_df = pmol.df['HETATM']
    full_atom_df = pd.concat([atom_df, hetatm_df])

    # find the minimum and maximum x, y, z coordinates for atoms in crystal ligand file
    x_min, x_max = min(full_atom_df.x_coord), max(full_atom_df.x_coord)
    y_min, y_max = min(full_atom_df.y_coord), max(full_atom_df.y_coord)
    z_min, z_max = min(full_atom_df.z_coord), max(full_atom_df.z_coord)

    # find center point for each axis
    x_center = x_min + abs(x_min-x_max)/2
    y_center = y_min + abs(y_min-y_max)/2
    z_center = z_min + abs(z_min-z_max)/2

    # calculate lengths for each axis based on difference between min and max values multiplied by user defined padding variable
    x_range = (abs(x_min-x_max)+padding)
    y_range = (abs(y_min-y_max)+padding)
    z_range = (abs(z_min-z_max)+padding)

    return x_center, y_center, z_center, x_range, y_range, z_range

def dock_file(docker_command, protein_filepath, ligand_filepath, center_x, center_y, center_z, size_x, size_y, size_z, ex, nw, nm, er, outfile): # dock the decoy pdbqt to the receptor pdbqt using GWOVina CLI
    os.system(f'{docker_command} --receptor {protein_filepath} --ligand {ligand_filepath}  --center_x  {center_x} --center_y {center_y} --center_z {center_z} --size_x  {size_x} --size_y {size_y}  --size_z {size_z}' \
              f' --exhaustiveness={ex} --num_wolves={nw} --num_modes={nm} --energy_range={er} --out {outfile} >&2')

def autodock_convert(ligand_filename, mgl_tools_path): # converts files from .pdb format to .pdbqt format using AutoDockTools

    ligand_filepath = os.path.join(os.path.dirname(__file__),'temp','pdb_files',ligand_filename)

    destination_path = os.path.join(os.path.dirname(__file__),'temp','pdbqt_files','')

    # define mgl_tools script paths
    pythonsh_path = f'{mgl_tools_path}bin/pythonsh'
    prepare_ligand_path = f'{mgl_tools_path}MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'

    # use AutoDockTools CLI to convert pdb file adding gasteiger charges and polar hydrogens
    try:
        output = subprocess.check_output(f'{pythonsh_path} {prepare_ligand_path} -l {ligand_filepath} -A hydrogens -o {destination_path}{ligand_filename}qt -U nphs', shell=True, stderr=subprocess.STDOUT)
    except:
        logging.warning(f'WARNING: Skipping {ligand_filepath} - MGLTools conversion error')

def merge_args(filepath1, filepaths):

    filepath1_list = [filepath1 for path in filepaths]

    merged_dict = dict(zip(filepaths, filepath1_list))

    return merged_dict


def make_pdbs_from_smiles(smi_tuple): # make and save pdb copies of active smiles files

    id = smi_tuple[0]

    smi = smi_tuple[1]

    result = 0

    try:

        # give the smile atoms coordinates and temporary hydrogens for addition of 3D structure
        ligand = Chem.MolFromSmiles(smi)

        if ligand is None:
            clean_smile = clean_smile_with_pybel(smi)
            ligand = Chem.MolFromSmiles(clean_smile)
            result = 1

        h_ligand = Chem.AddHs(ligand)
        AllChem.EmbedMolecule(h_ligand,randomSeed=0xf00d)

        # remove the temporary hydrogens from now 3D molecule
        ligand = Chem.RemoveHs(h_ligand)

        Chem.MolToPDBFile(ligand, f'utils/temp/pdb_files/{id}.pdb')

        return result

    except Exception as e:

       logging.warning(f'{smi}')
       logging.warning(f'WARNING: Could not convert {id} (above) to pdb: {str(e)}')

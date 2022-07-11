"""
SCORCH script version 1.0 - Run python scoring.py -h for help    
                                                                    
Script Authors:                                                     
@sammoneykyrle                                                      
@milesmcgibbon                                                      
                                                                    
School of Biological Sciences                                       
The University of Edinburgh                                         
"""

#######################################################################
# Imports

# import os and set tensorflow verbosity
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['NUMEXPR_MAX_THREADS'] = '1'

# import other libraries
import sys
import math
import json
import psutil
import shutil
import pickle
import logging
import argparse
import textwrap
import contextlib
import numpy as np
import pandas as pd
from tqdm import tqdm
import xgboost as xgb
import tensorflow as tf
from utils.ecifs import *
from functools import reduce
from utils import binana, kier
tf.get_logger().setLevel('ERROR')
from utils.dock_functions import *
from functools import partialmethod
from warnings import filterwarnings
from tensorflow.keras.models import load_model
from joblib import Parallel, parallel, delayed, load

#######################################################################
# Global Variables

# filter pandas warnings
filterwarnings('ignore')

# get working directory where scoring function is being deployed
stem_path = os.getcwd()

#######################################################################
# Functions

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):

    """
    Context manager to allow joblib
    progress monitoring
    """

    def tqdm_print_progress(self):
        if self.n_completed_tasks > tqdm_object.n:
            n_completed = self.n_completed_tasks - tqdm_object.n
            tqdm_object.update(n=n_completed)

    original_print_progress = parallel.Parallel.print_progress
    parallel.Parallel.print_progress = tqdm_print_progress

    try:
        yield tqdm_object
    finally:
        parallel.Parallel.print_progress = original_print_progress
        tqdm_object.close()


def run_binana(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get BINANA descriptors for    
    protein-ligand complex                  
                                        
    Inputs:  ligand as a pdbqt string block,         
    receptor pdbqt filepath                 
                                            
    Output: BINANA protein-ligand complex   
    descriptor features as a dictionary     
    """

    # empty dictionary to populate with features
    binana_features = dict()

    # get dictionary of binana features
    main_binana_out = binana.Binana(ligand_pdbqt_block, receptor_filepath).out


    # define the features we want
    keep_closest_contacts = ["2.5 (HD, OA)", 
                             "2.5 (HD, HD)", 
                             "2.5 (HD, N)", 
                             "2.5 (C, HD)", 
                             "2.5 (OA, ZN)", 
                             "2.5 (HD, ZN)", 
                             "2.5 (A, HD)"]
    
    keep_close_contacts = ["4.0 (C, C)", 
                           "4.0 (HD, OA)", 
                           "4.0 (C, HD)", 
                           "4.0 (C, N)", 
                           "4.0 (A, C)",
                           "4.0 (A, OA)", 
                           "4.0 (N, OA)", 
                           "4.0 (A, N)", 
                           "4.0 (HD, N)", 
                           "4.0 (HD, HD)", 
                           "4.0 (A, HD)", 
                           "4.0 (OA, OA)", 
                           "4.0 (C, OA)", 
                           "4.0 (N, N)",
                           "4.0 (C, SA)", 
                           "4.0 (HD, SA)", 
                           "4.0 (OA, SA)", 
                           "4.0 (N, SA)", 
                           "4.0 (A, A)", 
                           "4.0 (HD, S)", 
                           "4.0 (S, ZN)", 
                           "4.0 (N, ZN)", 
                           "4.0 (HD, ZN)", 
                           "4.0 (A, SA)", 
                           "4.0 (OA, ZN)", 
                           "4.0 (C, ZN)", 
                           "4.0 (C, NA)", 
                           "4.0 (NA, OA)", 
                           "4.0 (HD, NA)", 
                           "4.0 (N, NA)", 
                           "4.0 (A, NA)", 
                           "4.0 (BR, C)", 
                           "4.0 (HD, P)", 
                           "4.0 (F, N)", 
                           "4.0 (F, HD)", 
                           "4.0 (C, CL)", 
                           "4.0 (CL, HD)"]

    keep_ligand_atoms = ["LA N",
                         "LA HD"]
    
    keep_elsums = [ "ElSum (C, C)",
                    "ElSum (HD, OA)",
                    "ElSum (C, HD)",
                    "ElSum (C, N)",
                    "ElSum (A, C)",
                    "ElSum (A, OA)",
                    "ElSum (N, OA)",
                    "ElSum (A, N)",
                    "ElSum (HD, HD)",
                    "ElSum (A, HD)",
                    "ElSum (OA, OA)",
                    "ElSum (C, OA)",
                    "ElSum (N, N)",
                    "ElSum (C, SA)",
                    "ElSum (HD, SA)",
                    "ElSum (OA, SA)",
                    "ElSum (N, SA)",
                    "ElSum (A, A)",
                    "ElSum (N, S)",
                    "ElSum (HD, S)",
                    "ElSum (OA, S)",
                    "ElSum (A, SA)",
                    "ElSum (C, NA)",
                    "ElSum (NA, OA)",
                    "ElSum (HD, NA)",
                    "ElSum (N, NA)",
                    "ElSum (A, NA)",
                    "ElSum (BR, C)",
                    "ElSum (HD, P)",
                    "ElSum (OA, P)",
                    "ElSum (N, P)",
                    "ElSum (C, F)",
                    "ElSum (F, N)",
                    "ElSum (A, F)",
                    "ElSum (CL, OA)",
                    "ElSum (C, CL)",
                    "ElSum (CL, N)",
                    "ElSum (A, CL)"]

    # add closest contacts to binana_features dict
    for contact in keep_closest_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['closest'].get(binana_name)
    
    # add close contacts to binana_features dict
    for contact in keep_close_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['close'].get(binana_name)
    
    # add ligand atoms to binana_features dict as binary tallies
    for atom in keep_ligand_atoms:
        binana_name = atom.split()[-1]
        if main_binana_out['ligand_atoms'].get(binana_name) is None:
            binana_features[atom] = 0
        else:
            binana_features[atom] = 1
    
    # add electrostatics to binana_features dict
    for elsum in keep_elsums:
        binana_name = elsum.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[elsum] = main_binana_out['elsums'].get(binana_name)

    # add active site flexibility features to binana_features
    binana_features["BPF ALPHA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_ALPHA")
    binana_features["BPF ALPHA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_ALPHA")
    binana_features["BPF BETA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_BETA")
    binana_features["BPF BETA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_BETA")
    binana_features["BPF OTHER SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_OTHER")
    binana_features["BPF OTHER BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_OTHER")

    # add hydrophobic features to binana_features
    binana_features["HC ALPHA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_ALPHA")
    binana_features["HC ALPHA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_ALPHA")
    binana_features["HC BETA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_BETA")
    binana_features["HC BETA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_BETA")
    binana_features["HC OTHER SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_OTHER")
    binana_features["HC OTHER BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_OTHER")

    # add hydrogen bond features to binana_features
    binana_features["HB ALPHA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_ALPHA")
    binana_features["HB BETA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_OTHER")
    binana_features["HB ALPHA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_ALPHA")
    binana_features["HB ALPHA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_ALPHA")
    binana_features["HB BETA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_OTHER")

    # add salt bridge features to binana_features
    binana_features["SB ALPHA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_ALPHA")
    binana_features["SB BETA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_BETA")
    binana_features["SB OTHER"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_OTHER")

    # add aromatic stacking features to binana_features
    binana_features["piStack ALPHA"] = main_binana_out['stacking'].get("STACKING ALPHA")
    binana_features["piStack BETA"] = main_binana_out['stacking'].get("STACKING BETA")
    binana_features["piStack OTHER"] = main_binana_out['stacking'].get("STACKING OTHER")
    binana_features["tStack ALPHA"] = main_binana_out['t_stacking'].get("T-SHAPED_ALPHA")
    binana_features["tStack BETA"] = main_binana_out['t_stacking'].get("T-SHAPED_BETA")
    binana_features["tStack OTHER"] = main_binana_out['t_stacking'].get("T-SHAPED_OTHER")

    # add cation pi features to binana_features
    binana_features["catPi BETA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binana_features["catPi OTHER LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")

    # add rotatable bond count to binana features
    binana_features["nRot"] = main_binana_out['nrot']

    # return dictionary
    return binana_features

def kier_flexibility(ligand_pdbqt_block):

    """
    Function: Calculate Kier flexibility    
    for ligand                              
                                            
    Inputs: ligand as a pdbqt string block  
                                            
    Output: Kier flexibility                
    """

    # parse pdbqt block
    mol = kier.SmilePrep(ligand_pdbqt_block)

    # calculate flexibility
    return kier.CalculateFlexibility(mol)

def calculate_ecifs(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get ECIFs for protein-ligand  
    complex                                 
                                            
    Inputs: ligand as a pdbqt string block, 
    receptor pdbqt filepath                 
                                            
    Output: ECIF protein-ligand complex     
    descriptor features as a pandas DataFrame      
    """

    # get ECIFs with default cutoff using imported functions
    # from utils/ecifs.py
    ECIF_data = GetECIF(receptor_filepath, ligand_pdbqt_block, distance_cutoff=6.0)

    # replace the semicolons to make valid dataframe headers
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]

    # zip into a dictionary and convert to dataframe
    ECIF_data = dict(zip(ECIFHeaders,ECIF_data))
    ECIF_df = pd.DataFrame(ECIF_data,index=[0])

    # return the dataframe
    return ECIF_df

def extract(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get all descriptor features   
    for protein-ligand complex              
                                            
    Inputs: ligand as a pdbqt string block, 
    receptor pdbqt filepath     
                                            
    Output: All protein-ligand complex      
    descriptor features as a DataFrame      
    """
    
    # get the kier flexibility
    k = kier_flexibility(ligand_pdbqt_block)

    # get the binana descriptors and build into dataframe
    binana_dict = run_binana(ligand_pdbqt_block,receptor_filepath)
    binana_df = pd.DataFrame([binana_dict])

    # get the ECIFs
    ECIF = calculate_ecifs(ligand_pdbqt_block, receptor_filepath)

    # concatenate all feature columns to one row dataframe
    df = pd.concat([ECIF,binana_df],axis=1)

    # add the kier flexibility as a column
    df['Kier Flexibility'] = k

    # return the features
    return df

def prune_df_headers(df):

    """
    Function: Condense features for model input                             
                                            
    Inputs: Full Dataframe of               
    protein-ligand complex descriptors     
                                            
    Output: DataFrame of features for model input                                   
    """

    # load features we want in the correct order
    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    headers_58 = reference_headers.get('492_models_58')

    # subset the dataframe of features
    df = df[headers_58]

    # return the dataframe
    return df

def multiple_pose_check(ligand_filepath):

    """
    Function: Transform ligand.pdbqt file of      
    poses/models into pdbqt string blocks   
                                            
    Inputs: ligand.pdbqt filepath           
                                            
    Output: List of model/pose pdbqt string 
    blocks as a tuple with the pose number
    e.g. [('_pose_1','REMARK....)]                                 
    """

    # make empty list for populating
    pdbqt_pose_blocks = list()

    # open the input ligand file
    lig_text = open(ligand_filepath, 'r').read()

    # split by poses
    lig_poses = lig_text.split('MODEL')

    # for each pose clean up any whitespace or empty lines
    for pose in lig_poses:
        lines = pose.split('\n')
        clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]

        # if there are less than three lines then
        # its an artefact of the splitting and we ignore 
        if len(clean_lines) < 3:
            pass

        # otherwise join up the correct poses
        else:
            pose = '\n'.join(clean_lines)

            # add the pose to the list of poses
            pdbqt_pose_blocks.append(pose)

    # map the poses into a list of tuples
    pdbqt_pose_blocks = list(map(lambda x: (f'_pose_{pdbqt_pose_blocks.index(x) + 1}', x), pdbqt_pose_blocks))

    # return the list of tuples
    return pdbqt_pose_blocks

def run_networks(df, models_to_load, model_name):

    """
    Function: Get single mean prediction from consensus of 
    neural network models for a protein-ligand complex              
                                            
    Inputs: Dataframe of features for model,
            list of model filenames to load,
            model name string for prediction columns

    Output: Dataframe of model predictions with 
            final average mean column              
    """

    # make empty dataframe to populate wit predictions
    predictions = pd.DataFrame()

    # for each model
    for i in tqdm(range(len(models_to_load))):

        # load the model and get its prediction on the data
        model = load_model(models_to_load[i])
        y_pred = model.predict(df)

        # add the prediction as a column to the predictions df
        predictions[f'{model_name}_{i + 1}'] = y_pred.flatten()

    # take an average of all the predictions
    predictions[f'{model_name}_models_average'] = predictions.mean(axis=1)

    # return the df of predictions
    return predictions


def binary_concat(dfs, headers):

    """
    Function: Concatenate list of           
    dataframes into a single dataframe by   
    sequentially writing to a single binary 
    file (removes pd.concat bottleneck)     
                                            
    Inputs: List of dataframes, dataframe   
    headers as a list                       
                                            
    Output: Single combined dataframe       
    """

    # set total rows count
    total_rows = 0

    # make temporary directory for binary file if it
    # doesn't exist yet
    if not os.path.isdir(os.path.join('utils','temp')):
        os.makedirs(os.path.join('utils','temp'))

    # create a temporary binary file
    with open(os.path.join('utils','temp','features.bin'),'wb+') as binary_store:

        # then for each dataframe write it to binary
        for df in dfs:

            # make sure nRot is numeric in type
            df['nRot'] = pd.to_numeric(df['nRot'])

            # get the rows and columns of the dataframe
            rows, fixed_total_columns = df.shape

            # add the rows to total rows counter
            total_rows += rows

            # write the df to the binary store file
            binary_store.write(df.values.tobytes())

            # store the dtypes
            typ = df.values.dtype

    # open the temporary binary file
    with open(os.path.join('utils','temp','features.bin'),'rb') as binary_store:

        # read the binary file
        buffer = binary_store.read()

        # shape the data into a dataframe using the info from
        # when we wrote all the files to the binary store
        data = np.frombuffer(buffer, dtype=typ).reshape(total_rows, fixed_total_columns)
        master_df = pd.DataFrame(data = data, columns = headers)

    # make sure the binary store is deleted in case we call the function again
    os.remove(os.path.join('utils','temp','features.bin'))

    # return the concatenated dataframe
    return master_df

def parse_module_args(args_dict):

    """
    Function: Parse user arguments when     
    script is imported as a module          
                                            
    Inputs: User arguments as a dictionary  
                                            
    Output: Populated params dictionary     
    """

    # empty list to use as a spoof sys.argv result
    command_input = list()

    # check if any boolean flag arguments have been passed
    boolean_args = ['verbose','return_pose_scores']
    for key, value in args_dict.items():
        if key in boolean_args:
            if value:
                command_input.append(f'-{key}')
        
        # otherwise add them as normal args with their values
        else:
            command_input.append(f'-{key}')
            command_input.append(str(value))

    # then parse the args as if they were command line
    parsed_args = parse_args(command_input)

    # return the arguments
    return parsed_args

def parse_args(args):

    """
    Function: Parse user defined command    
    line arguments                          
                                            
    Inputs: Command line arguments          
                                            
    Output: Populated params dictionary     
    """

    parser = argparse.ArgumentParser(description="SCORCH 1.0")

    requiredNamed = parser.add_argument_group('required named arguments')
    
    requiredNamed.add_argument('-l','--ligand', help="""Ligands to score against the supplied receptor. 
                                                 Can be a .smi or .txt filepath, a .pdbqt filepath, or path to a folder of pdbqt files. 
                                                 If .smi file is supplied, --range and --center args or --ref_lig args are 
                                                 also required.""", required=True)
    requiredNamed.add_argument('-r','--receptor', help="Receptor to score ligands against. Must be a filepath to a .pdbqt file", required=True)
    parser.add_argument('-rl','--ref_lig', help="Filepath to example ligand in receptor binding site (mol, mol2, sdf, pdb or pdbqt)")
    parser.add_argument('-t','--threads', default=1, help="Number of CPU threads to parallelise SCORCH over", type=int)
    parser.add_argument('-c','--center', help="'[x, y, z]' coordinates of the center of the binding site for docking")
    parser.add_argument('-ra','--range', help="'[x, y, z]' axis lengths to define a box around --center coordinates for docking")
    parser.add_argument('-o','--out', help="Filepath for output csv (If not supplied, scores are written to stdout)")
    parser.add_argument('-p','--return_pose_scores', action='store_true', help="If supplied, scoring values for individual poses in each ligand file are returned")
    parser.add_argument('-v','--verbose', action='store_true', help="If supplied, progress bars and indicators are displayed while scoring")
    parser.add_argument('-s','--pose_1', action='store_true', help="Consider only the first pose in each pdbqt file to score - NOT RECOMMENDED")
    parser.add_argument('-d','--dock', action='store_true', help="""If supplied, input ligands are assumed to be text SMILES and will be docked 
                                                                    using GWOVina before scoring. This will be autodetected if a .smi or .txt file is supplied""")
    params = parser.parse_args()

    if params.verbose:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    else:
        tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)
        logging.basicConfig(level=logging.WARNING, format='%(message)s')

    if os.path.isdir(params.ligand):
        params.ligand = [os.path.join(params.ligand, file) for file in os.listdir(params.ligand)]
        receptors = [params.receptor for i in range(len(params.ligand))]
        params.receptor = receptors

    elif '.smi' in params.ligand or '.txt' in params.ligand:
        params.dock = True
    else:
        params.ligand = [params.ligand]
        params.receptor = [params.receptor]

    return params

def prepare_and_dock_inputs(params):

    """
    Function: Prepare smiles inputs for scoring by converting
              to pdbqt files and docking with GWOVina

    Inputs:   Parsed command line parameters

    Outputs:  Updated parameters where the ['ligand'] is docked pdbqt versions
              of supplied smiles strings
              Dictionary of smiles strings with their identifiers
    """

    # load the docking settings
    dock_settings = json.load(open(os.path.join('utils','params','dock_settings.json')))

    # check there has been a binding site location passed to use for docking
    if params.ref_lig is None:
        if params.center is None and params.range is None:

            # if not then complain and exit
            logging.critical("ERROR: No reference ligand or binding site coordinates supplied. Try:\n- ensuring center and range values are entered correctly\n- supplying a reference ligand")
            sys.exit()

        # if center and range have been supplied then parse them into
        # coordinates variable
        else:
            try:
                coords = (float(params.center[0]),
                            float(params.center[1]),
                            float(params.center[2]),
                            float(params.range[0]),
                            float(params.range[1]),
                            float(params.range[2]))
            
            # if this doesn't work then complain and exit
            except:
                logging.critical("\nERROR: Binding site coordinates for docking are missing or incorrectly defined. Try:\n- ensuring center and range values are entered correctly\n- using a reference ligand instead")
                sys.exit()

    # if the user has passed a reference ligand 
    else:
        coords = get_coordinates(params.ref_lig, dock_settings['padding'])

    # make sure all temporary folders exist
    if not os.path.isdir(os.path.join('utils','temp','pdb_files')):
        os.makedirs(os.path.join('utils','temp','pdb_files'))
        os.makedirs(os.path.join('utils','temp','pdbqt_files'))
        os.makedirs(os.path.join('utils','temp','docked_pdbqt_files'))

    # make sure the folders are empty if they exist
    pdbs = get_filepaths(os.path.join('utils','temp','pdb_files',''))
    for pdb in pdbs:
        os.remove(pdb)

    pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))
    for pdbqt in pdbqts:
        os.remove(pdbqt)

    docked_pdbqts = get_filepaths(os.path.join('utils','temp','docked_pdbqt_files',''))
    for docked_pdbqt in docked_pdbqts:
        os.remove(docked_pdbqt)

    # load the smiles from the input file as a dictionary
    # with an identifier as the key and the smile string as the value
    smi_dict = get_smiles(params.ligand)

    logging.info('Generating 3D pdbs from SMILES...')

    # parallelise the conversion of smiles strings to 3D pdbqts with RDKit
    with tqdm_joblib(tqdm(desc="Generating...", total=len(smi_dict))) as progress_bar:
        Parallel(n_jobs=params.threads)(delayed(make_pdbs_from_smiles)(smi) for smi in smi_dict.items())

    # get a list of all successfully converted pdbs
    pdbs = os.listdir(os.path.join('utils','temp','pdb_files',''))

    logging.info('Converting pdbs to pdbqts...')

    # parallelise conversion of pdbs to pdbqt files
    with tqdm_joblib(tqdm(desc="Converting...", total=len(pdbs))) as progress_bar:
        Parallel(n_jobs=params.threads)(delayed(autodock_convert)(pdb_file, os.path.join('utils','MGLTools-1.5.6','')) for pdb_file in pdbs)

    # get list of successfully converted pdbqts
    pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))

    # get os name
    if sys.platform.lower() == 'darwin':
        os_name = 'mac'
    elif 'linux' in sys.platform.lower():
        os_name = 'linux'

    logging.info("Docking pdbqt ligands...")

    # then for each pdbqt
    for pdbqt in tqdm(pdbqts):

        # then dock each one
        dock_file(
                    os.path.join('utils','gwovina-1.0','build',os_name,'release','gwovina'),
                    params.receptor,
                    pdbqt,
                    *coords,
                    dock_settings['gwovina_settings']['exhaustiveness'],
                    dock_settings['gwovina_settings']['num_wolves'],
                    dock_settings['gwovina_settings']['num_modes'],
                    dock_settings['gwovina_settings']['energy_range'],
                    outfile=os.path.join(f'{stem_path}','utils','temp','docked_pdbqt_files',f'{os.path.split(pdbqt)[1]}')
                    )

    # get name of the input smiles file
    if '.' in params.ligand:
        docked_ligands_folder = os.path.basename(params.ligand).split('.')[0]
    else:
        docked_ligands_folder = os.path.basename(params.ligand)

    # define folder with the input file name to store the docked ligands
    docked_ligands_path = os.path.join('docked_ligands',docked_ligands_folder,'')

    # add each docked ligand in temp file as a ligand to score to the list of ligands in params.ligand
    params.ligand = [os.path.join('utils','temp','docked_pdbqt_files', file) for file in os.listdir(os.path.join('utils','temp','docked_pdbqt_files'))]
    
    # build receptor as a repeating list into params dict
    receptors = [params.receptor for i in range(len(params.ligand))]
    params.receptor = receptors

    # make sure the docked ligands folder exists
    if not os.path.isdir('docked_ligands'):
        os.mkdir('docked_ligands')
    if not os.path.isdir(docked_ligands_path):
        os.makedirs(docked_ligands_path)
    
    # copy the docked ligands into a main folder so they are accessible
    for file in params.ligand:
        shutil.copy(file, docked_ligands_path)   
    
    # return the updated parameters with docked ligands to score and the
    # dictionary of ids to smiles strings
    return params, smi_dict

def count_input_poses(list_of_ligands):

    """
    Function: Count the total number of poses across
              all ligands that need scoring

    Inputs:   List of ligand pdbqt files

    Outputs:  Total poses to score as integer
    """

    # set up count
    total_poses = 0

    # then for each ligand
    for ligand in tqdm(list_of_ligands):
        
        # open it and count the number of models
        ligand_file = open(ligand).read()
        poses = ligand_file.count('MODEL')

        # add the poses to total pose count
        if poses == 0:
            total_poses += 1
        else:
            total_poses += poses

    # return integer of total poses
    return total_poses

def ligand_pose_generator(params, lower_index, upper_index):

    """
    Function: Generate a list of receptor and ligand arguments
              between specific indices to pass to a function to be scored.
              This is necessary because if we do this at the start, if
              many ligands are being scored we cannot store them all in 
              memory and the script crashes

    Inputs:   Command line arguments, the lower and upper indexes
              of ligand poses to score

    Outputs:  List of tuples to score in the form
              [(receptor filename, ligand filename, (ligand pose number, ligand pdbqt block))]
    """

    # set up list to populate
    requested_receptor_ligand_args = list()

    # track the pose index
    pose_index = 0

    # then for each ligand
    for ligand_index, ligand_filepath in enumerate(params.ligand):

         
         # don't waste any time if we're already over the upper index
        if pose_index > upper_index:
            break

        # load the poses from the current ligand being considered
        pdbqt_pose_blocks = list()
        lig_text = open(ligand_filepath, 'r').read()
        lig_poses = lig_text.split('MODEL')
        for pose in lig_poses:
            lines = pose.split('\n')
            clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
            if len(clean_lines) < 3:
                pass
            else:
                pose = '\n'.join(clean_lines)
                pdbqt_pose_blocks.append(pose)

                # stop if we only want one pose
                if params.pose_1:
                    break

        # make a tuple with pdbqt block and pose name
        poses = [(f'_pose_{pdbqt_pose_blocks.index(pose) + 1}', pose) for pose in pdbqt_pose_blocks]

        # for each pose
        for pose in poses:
            
            # if the pose is one we want
            if lower_index <= pose_index < upper_index:
            
                # add it to the receptor ligand arguments
                receptor_ligand_args = (params.receptor[ligand_index], ligand_filepath, pose)

                requested_receptor_ligand_args.append(receptor_ligand_args)
    
            # update the pose index
            pose_index += 1

    # return the requested poses between the specified indexes
    return requested_receptor_ligand_args

def calculate_batches_needed(total_poses):

    """
    Function: Calculate how many batches the ligands
              need to be scored in with the available
              memory
    
    Input:    Total number of poses as an integer

    Output:   Number of batches to split the ligand poses into
    """

    # estimate the ram usage from empirical data
    estimated_ram_usage = (360540*total_poses) + 644792975
    available_ram = psutil.virtual_memory().total
    safe_ram_available = available_ram*0.7

    # use this info to split the batches into least possible
    # batches without getting a memory error
    if estimated_ram_usage > safe_ram_available:
        batches_needed = math.ceil(estimated_ram_usage/safe_ram_available)
    else:
        batches_needed = 1

    # return batch number as integer
    return batches_needed


def list_to_chunk_indexes(list_length, number_of_chunks):

    """
    Function: Create nested list of upper and lower indexes
              which relate to batches of input list

    Inputs:   List length (integer) and number of batches

    Outputs:  Nested list of chained indexes e.g.
              [(0, 573),(573, 1092)]
    """

    # set up list to populate
    indices = list()

    # get the size of each batch
    chunksize = math.ceil(list_length/number_of_chunks)

    # then for each chunk taking chunksize steps
    for i in range(0, list_length, chunksize):

        # if its the last chunk then it ends at the end of the list
        if i+chunksize < list_length:
            indices.append((i, i+chunksize))

        # otherwise it spans one chunksize chunk
        else:
            indices.append((i, list_length))

    # return list of indices
    return indices

def prepare_features(receptor_ligand_args):

    filterwarnings('ignore')

    """
    Function: Wrapper to prepare all requested protein-ligand            
              complexes/poses for scoring             
                                            
    Inputs:   A tuple of a receptor ligand pair to score                            
                                            
    Output:   A single row dataframe of protein-ligand complex
              features                      
    """

    # grab the ligand pose number and pdbqt string block
    ligand_pose_number = receptor_ligand_args[2][0]
    ligand_pdbqt_block = receptor_ligand_args[2][1]
     
    # grab the ligand and receptor filepaths and filenames
    receptor_filepath = receptor_ligand_args[0]
    ligand_filepath = receptor_ligand_args[1]
    ligand_basename = os.path.basename(ligand_filepath)
    ligand_basename = ligand_basename.replace('.pdbqt', ligand_pose_number)
    receptor_basename = os.path.basename(receptor_filepath)

    # extract the interaction features
    features = extract(ligand_pdbqt_block, receptor_filepath)

    # prune the headers down to those needed for model scoring
    multi_pose_features = prune_df_headers(features)

    # fill None values with 0 for binana features
    multi_pose_features.fillna(0, inplace=True)

    # add receptor and ligand info to features
    multi_pose_features['Receptor'] = receptor_basename
    multi_pose_features['Ligand'] = ligand_basename
    
    # return the dataframe
    return multi_pose_features

def scale_multipose_features(df):

    """
    Function: Scale features using prefitted feature scaler

    Input:    Dataframe of features to scale

    Output:   Scaled dataframe of features
    """

    # open the reference features file
    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    scaler_58 = reference_headers.get('for_scaler_58')
    headers_58 = reference_headers.get('492_models_58')

    # store the ligand and receptor information
    ligands, receptors = df['Ligand'], df['Receptor']

    # get the missing columns that the scaler expects
    missing_columns = list(set(scaler_58) - set(list(df)))

    # fill in dummy columns for the scaler
    for col in missing_columns:
        df[col] = 0
    df = df[scaler_58]

    # load the scaler and scale the data
    scaler = load(os.path.join('utils','params','58_maxabs_scaler_params.save'))
    scaled = scaler.transform(df)
    df[df.columns] = scaled

    # then only keep the columns we want for the models
    df = df[headers_58]

    # add the receptor and ligand info back in
    df['Ligand'], df['Receptor'] = ligands, receptors

    # return the dataframe
    return df

def score(models, features):

    """
    Function: Score supplied ligands with   
    an individual model                     
                                            
    Inputs: Tuple of (model_name,           
                      model_binary_file,    
                      feature dataframes)   
                                            
    Output: Dataframe of model predictions  
    """

    # get the variables out of the tuple
    model_name = models[0]
    model_file = models[1]

    logging.info(f'Scoring with {model_name}...')

    # create results and features dataframe
    results = features[['Ligand','Receptor']].copy().reset_index(drop=True)
    df = features.drop(['Ligand','Receptor'], axis=1)

    # get the model scores for each row in the dataframe
    if 'xgb' in model_name:
        dtest = xgb.DMatrix(df, feature_names=df.columns)
        results[model_name] = model_file.predict(dtest)

    else:
        network_predictions = run_networks(df, model_file, model_name)
        results[network_predictions.columns] = network_predictions[network_predictions.columns]

    # return prediction results
    return results

def print_intro(params):

    """
    Function: Prints chosen arguments to    
    stdout                                  
                                            
    Inputs: User command line parameters    
    dictionary                              
                                            
    Output: None                            
    """


    logging.info('\n')
    logging.info('**************************************************************************')

    logging.info('SCORCH v1.0')
    logging.info('Miles McGibbon, Samuel Money-Kyrle, Vincent Blay & Douglas R. Houston\n')

    logging.info('**************************************************************************\n')

    if not params.dock:
        logging.info(f'Found {len(params.ligand)} ligand(s) for scoring against a single receptor...\n')

    else:
        ligand_count = open(params.ligand).read().split("\n")
        ligand_count = len([l for l in ligand_count if l != ''])
        logging.info(f'Parsed {ligand_count} ligand smiles for docking and scoring against a single receptor...\n')

        logging.info('**************************************************************************\n')

def prepare_models(params):

    """
    Function: Loads machine-learning model  
    binaries                                
                                            
    Inputs: User command line parameters    
    dictionary                              
                                            
    Output: Dictionary of {model_name:      
                           model_binary}    
    """

    logging.info('\n**************************************************************************\n')
    logging.info('Model Request Summary:\n')

    # empty dictionary to store models
    models = {}

    # load xgboost model
    logging.info('XGBoost Model: Yes')
    xgb_path = os.path.join('utils','models','xgboost_models','495_models_58_booster.pkl')
    models['xgboost_model'] = pickle.load(open(xgb_path,'rb'))

    # load best 15 neural network models for WD models and feedforward models
    logging.info('Feedforward NN Models : Yes')
    models['ff_nn'] = os.path.join('utils','models','ff_nn_models')
    model_ranks = pickle.load(open(os.path.join(models['ff_nn'],'rankings.pkl'),'rb'))
    model_ranks = model_ranks[:15]
    models['ff_nn'] = [os.path.join(models['ff_nn'], 'models',f'{model[1]}.hdf5') for model in model_ranks]

    logging.info('W&D NN Models : Yes')
    models['wd_nn'] = os.path.join('utils','models','wd_nn_models')
    model_ranks = pickle.load(open(os.path.join(models['wd_nn'],'rankings.pkl'),'rb'))
    model_ranks = model_ranks[:15]
    models['wd_nn'] = [os.path.join(models['wd_nn'], 'models',f'{model[1]}.hdf5') for model in model_ranks]
    logging.info('\n')

    if params.pose_1:

        logging.info('Calculating scores for first pose only in pdbqt file(s)\n')

    # return dictionary of model binaries
    return models

def score_ligand_batch(params, ligand_batch, model_binaries):

    """
    Function:  Make a dataframe of scores and stats for a batch of
               ligands

    Inputs:    Parsed command line parameters, batch of receptor ligand tuples,
               and the binary files for the models to use
    
    Outputs:   Dataframe of scores and stats for ligands in the batch
    """

    # multiprocess the extracting of features from the protein ligand pairs
    with tqdm_joblib(tqdm(desc="Preparing features", total=len(ligand_batch))) as progress_bar:
        multi_pose_features = Parallel(n_jobs=params.threads)(delayed(prepare_features)(ligand) for ligand in ligand_batch)

    # concatenate all the produced features
    multi_pose_features = pd.concat(multi_pose_features)

    # scale the features
    multi_pose_features = scale_multipose_features(multi_pose_features)

    # load the models
    models = [(m[0], m[1]) for m in model_binaries]

    model_results = list()

    logging.info('**************************************************************************\n')

    # score the features and add dataframe of results to list
    for model in models:
        model_results.append(score(model, multi_pose_features))
        logging.info('Done!')

    logging.info('**************************************************************************\n')

    # merge the results into one dataframe
    merged_results = reduce(lambda x, y: pd.merge(x, y, on = ['Receptor','Ligand']), model_results)

    # create main scorch score by taking mean of model scores
    multi_models = ['xgboost_model',
                    'ff_nn_models_average',
                    'wd_nn_models_average']
    merged_results['SCORCH_pose_score'] = merged_results[multi_models].mean(axis=1)

    # calculate scorch certainty score
    max_std = 0.4714 # result from [0, 0, 1] or [1, 1, 0]
    minimum_val = 1 - max_std
    merged_results['SCORCH_stdev'] = merged_results[multi_models].std(axis=1, ddof=0)
    merged_results['SCORCH_certainty'] = ((1-merged_results['SCORCH_stdev'])-minimum_val)/max_std

    # subset only columns we want
    merged_results = merged_results[['Receptor',
                                        'Ligand',
                                        'SCORCH_pose_score',
                                        'SCORCH_certainty']].copy()

    # add ligand id and pose number columns
    merged_results['Ligand_ID'] = merged_results['Ligand'].apply(lambda x: x.split('_pose')[0])
    merged_results['Pose_Number'] = merged_results['Ligand'].apply(lambda x: x.split('_pose_')[-1])

    del merged_results['Ligand']

    # return final dataframe of results
    return merged_results

def create_final_results(params, ligand_scores):

    """
    Function: Create final results output after batches have been scored

    Inputs:   Parsed command line parameters, list of dataframes of ligand scores

    Outputs:  Single dataframe of final results from whole scoring run
    """

    # concatenate the ligand scores and calculate the best scoring pose per ligand
    final_ligand_scores = pd.concat(ligand_scores)
    final_ligand_scores['SCORCH_score'] = final_ligand_scores.groupby(['Ligand_ID'])['SCORCH_pose_score'].transform('max')
    final_ligand_scores['best_pose'] = np.where(final_ligand_scores.SCORCH_score == final_ligand_scores.SCORCH_pose_score, 1, 0)

    # if the user doesnt want pose scores then remove all but the best scoring pose
    if not params.return_pose_scores:
        final_ligand_scores = final_ligand_scores.loc[final_ligand_scores.best_pose == 1]

        # only return specified columns
        final_ligand_scores = final_ligand_scores[['Receptor',
                                                    'Ligand_ID',
                                                    'Pose_Number',
                                                    'SCORCH_score',
                                                    'SCORCH_certainty']].copy()

        # sort the dataframe in order of best ligands first
        final_ligand_scores = final_ligand_scores.sort_values(by='SCORCH_score', ascending=False)

    # if user wants pose scores returned then sort in order of best ligands, and each ligand pose also sorted in order
    else:
        final_ligand_scores = final_ligand_scores.sort_values(by=['SCORCH_score','SCORCH_pose_score'], ascending=False)
    
    # then round all the scores to 3 decimal places
    numerics = list(final_ligand_scores.select_dtypes(include=[np.number]))
    final_ligand_scores[numerics] = final_ligand_scores[numerics].round(5)

    # return the final results
    return final_ligand_scores

def scoring(params):

    """
    Function: Score protein-ligand complex(es)                             
                                            
    Inputs: User command line parameters dictionary                              
                                            
    Output: Dataframe of scoring function predictions                             
    """

    # print help/intro if requested
    print_intro(params)

    # prepare and dock smiles if smiles ligands supplied
    if params.dock:

       params, smi_dict = prepare_and_dock_inputs(params)

    logging.info('**************************************************************************\n')
    
    logging.info('Counting input poses...')

    
    # count input poses
    if params.pose_1:
        total_poses = len(params.ligand)
    else:
        total_poses = count_input_poses(params.ligand)

    # calculate how many scoring batches we need with the system memory available
    batches_needed = calculate_batches_needed(total_poses)

    # calculate pose indexes to score in each batch
    ligand_batch_indexes = list_to_chunk_indexes(total_poses, batches_needed)

    # prepare the models
    model_dict = prepare_models(params)
    model_binaries = list(model_dict.items())

    # empty list for scores from each batch
    ligand_scores = list()

    logging.info('**************************************************************************\n')
 
    # then score ligands in batches
    for batch_number, index_range in enumerate(ligand_batch_indexes):

        logging.info(f"Scoring ligand batch {batch_number + 1} of {batches_needed}\n")

        logging.info('**************************************************************************\n')

        logging.info('Loading ligand batch poses...')

        # load in ligands to score for this batch
        ligand_batch = ligand_pose_generator(params, index_range[0], index_range[1])
        
        # score the batch and get the results
        merged_results = score_ligand_batch(params, ligand_batch, model_binaries)

        # add the results to the ligand_scores list
        ligand_scores.append(merged_results)

    # make final results from ligand scores
    final_ligand_scores = create_final_results(params, ligand_scores)

    # if we have docked the ligands then add column for their smiles strings
    if params.dock:
        final_ligand_scores['Ligand_SMILE'] = final_ligand_scores['Ligand_ID'].map(smi_dict)


    # return the final dataframe of ligand scores
    return final_ligand_scores

#######################################################################
# Main script

if __name__ == "__main__":

    # parse user arguments
    params = parse_args(sys.argv)

    # score complexes
    scoring_function_results = scoring(params)

    # output results
    if not params.out:

        # send to stdout if no outfile given
        sys.stdout.write(scoring_function_results.to_csv(index=False))

    else:

        # otherwise save to user specified csv
        scoring_function_results.to_csv(params.out, index=False)


# scorch.py end
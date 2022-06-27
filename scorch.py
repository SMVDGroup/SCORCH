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

# import other libraries
import sys
import math
import json
import psutil
import shutil
import pickle
import logging
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
from itertools import product, chain
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

def get_ligand_id(ligand):
    """
    Function: Get ligand ID from pose scores
                                            
    Inputs: Ligand pose name                
                                            
    Output: Ligand ID base name             
    """

    Ligand_ID = ligand.split('_pose')[0]
    return Ligand_ID

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

    mol = kier.SmilePrep(ligand_pdbqt_block)
    return kier.CalculateFlexibility(mol)

def calculate_ecifs(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get ECIFs for protein-ligand  
    complex                                 
                                            
    Inputs: ligand as a pdbqt string block, 
    receptor pdbqt filepath                 
                                            
    Output: ECIF protein-ligand complex     
    descriptor features as a DataFrame      
    """

    ECIF_data = GetECIF(receptor_filepath, ligand_pdbqt_block, distance_cutoff=6.0)
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]
    ECIF_data = dict(zip(ECIFHeaders,ECIF_data))
    ECIF_df = pd.DataFrame(ECIF_data,index=[0])

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
    
    k = kier_flexibility(ligand_pdbqt_block)
    binana_dict = run_binana(ligand_pdbqt_block,receptor_filepath)
    binana_df = pd.DataFrame([binana_dict])
    ECIF = calculate_ecifs(ligand_pdbqt_block, receptor_filepath)
    df = pd.concat([ECIF,binana_df],axis=1)
    df['Kier Flexibility'] = k

    return df

def prune_df_headers(df):

    """
    Function: Condense and features for     
    model input                             
                                            
    Inputs: Full Dataframe of               
    protein-ligand complex descriptors,     
    boolean for single pose model type,     
    boolean for further condensation with   
    principle component analysis            
                                            
    Output: DataFrame of features for model 
    input                                   
    """

    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    headers_58 = reference_headers.get('492_models_58')
    df = df[headers_58]

    return df

def multiple_pose_check(ligand_filepath, pose_1):

    """
    Function: Transform ligand.pdbqt        
    poses/models into pdbqt string blocks   
                                            
    Inputs: ligand.pdbqt filepath           
                                            
    Output: List of model/pose pdbqt string 
    blocks                                  
    """

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

    pdbqt_pose_blocks = list(map(lambda x: (f'_pose_{pdbqt_pose_blocks.index(x) + 1}', x), pdbqt_pose_blocks))

    return pdbqt_pose_blocks

def run_networks(df, model_file, model_name):

    models_to_load = model_file

    """
    Function: Get prediction from MLP model 
    for protein-ligand complex              
                                            
    Inputs: Number of networks as integer,  
    condensed protein-ligand complex        
    features as DataFrame                   
                                            
    Output: Float prediction                
    """


    predictions = pd.DataFrame()
    model_columns = list()

    for i in tqdm(range(len(models_to_load))):
        model = load_model(models_to_load[i])
        y_pred = model.predict(df)
        model_columns.append(f'{model_name}_{i + 1}')
        predictions[f'{model_name}_{i + 1}'] = y_pred.flatten()

    predictions[f'{model_name}_models_average'] = predictions[model_columns].mean(axis=1)

    return predictions.reset_index(drop=True)

def run_xgboost_models(df):

    """
    Function: Get prediction from XGB model 
    for protein-ligand complex              
                                            
    Inputs: Condensed protein-ligand        
    complex features as DataFrame,          
    boolean for single pose model           
                                            
    Output: Float prediction                
    """

    global xgboost_models
    dtest = xgb.DMatrix(df, feature_names=df.columns)
    prediction = xgboost_models.predict(dtest)
    return prediction


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

    total_rows = 0
    if not os.path.isdir(os.path.join('utils','temp')):
        os.makedirs(os.path.join('utils','temp'))

    with open(os.path.join('utils','temp','features.bin'),'wb+') as binary_store:
        for df in dfs:
            df['nRot'] = pd.to_numeric(df['nRot'])
            rows, fixed_total_columns = df.shape
            total_rows += rows
            binary_store.write(df.values.tobytes())
            typ = df.values.dtype

    with open(os.path.join('utils','temp','features.bin'),'rb') as binary_store:
        buffer = binary_store.read()
        data = np.frombuffer(buffer, dtype=typ).reshape(total_rows, fixed_total_columns)
        master_df = pd.DataFrame(data = data, columns = headers)
    os.remove(os.path.join('utils','temp','features.bin'))
    return master_df

def parse_module_args(args_dict):

    """
    Function: Parse user arguments when     
    script is imported as a module          
                                            
    Inputs: User arguments as a dictionary  
                                            
    Output: Populated params dictionary     
    """


    command_input = list()

    boolean_args = ['verbose','return_pose_scores']

    for key, value in args_dict.items():
        if key in boolean_args:
            if value:
                command_input.append(f'-{key}')
        else:
            command_input.append(f'-{key}')
            command_input.append(str(value))

    parsed_args = parse_args(command_input)

    return parsed_args

def parse_args(args):

    """
    Function: Parse user defined command    
    line arguments                          
                                            
    Inputs: Command line arguments          
                                            
    Output: Populated params dictionary     
    """

    params = {}

    if '-h' in args:
        prefix = "\t\t"
        expanded_indent = textwrap.fill(prefix+'$', replace_whitespace=False)[:-1]
        subsequent_indent = ' ' * len(expanded_indent)
        wrapper = textwrap.TextWrapper(initial_indent=prefix,
                                       subsequent_indent=subsequent_indent)
        with open(os.path.join('utils','help_string.txt')) as help_string:
            help = help_string.read()
            for line in help.split('\n'):
                if line.isupper():
                    print(line)
                elif  '-' in line:
                    print(line)
                else:
                    print(wrapper.fill(line))
        sys.exit()

    try:
        params['ligand'] = args[args.index('-ligand') + 1]
        params['receptor'] = args[args.index('-receptor') + 1]
        try:
            params['threads'] = int(args[args.index('-threads') + 1])
        except:
            params['threads'] = 1
        try:
            params['ref_lig'] = args[args.index('-ref_lig') + 1]
        except:
            params['ref_lig'] = None
        try:
            params['center'] = json.loads(args[args.index('-center') + 1])
        except:
            params['center'] = None
        try:
            params['range'] = json.loads(args[args.index('-range') + 1])
        except:
            params['range'] = None
        try:
            params['out'] = args[args.index('-out') + 1]
        except:
            params['out'] = False

        params['screen'] = False
        params['single'] = False
        params['pose_1'] = False
        params['return_pose_scores'] = False
        params['dock'] = False
        params['concise'] = True
        params['num_networks'] = 15

        if '-verbose' in args:
            params['verbose'] = True
            logging.basicConfig(level=logging.INFO, format='%(message)s')
        else:
            params['verbose'] = False
            tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)
            logging.basicConfig(level=logging.WARNING, format='%(message)s')

        if '-return_pose_scores' in args:
            params['return_pose_scores'] = True

        if '-pose_1' in args:
            params['pose_1'] = True

        if '-detailed' in args:
            params['concise'] = False

        elif os.path.isdir(params['ligand']) == True:
            params['ligand'] = [os.path.join(params['ligand'], file) for file in os.listdir(params['ligand'])]
            receptors = [params['receptor'] for i in range(len(params['ligand']))]
            params['receptor'] = receptors
            params['screen'] = True

        elif '.smi' in params['ligand'] or '.txt' in params['ligand']:
            params['dock'] = True

        else:
            params['ligand'] = [params['ligand']]
            params['receptor'] = [params['receptor']]
            params['single'] = True

        if '-num_networks' in args:
            params['num_networks'] = int(args[args.index('-num_networks') + 1])


    except ValueError as e:
        if 'is not in list' in str(e):
            missing = str(e).replace("'",'').replace(' is not in list','')
            logging.critical(f'Error: essential argument {missing} not supplied')
            logging.critical('Run "python scoring.py -h" for usage instructions')
            sys.exit()

    return params

def prepare_and_dock_inputs(params):

    dock_settings = json.load(open(os.path.join('utils','params','dock_settings.json')))

    if params['ref_lig'] is None:
        if params['center'] is None and params['range'] is None:
            logging.critical("ERROR: No reference ligand or binding site coordinates supplied. Try:\n- ensuring center and range values are entered correctly\n- supplying a reference ligand")
            sys.exit()
        else:
            try:
                coords = (float(params['center'][0]),
                            float(params['center'][1]),
                            float(params['center'][2]),
                            float(params['range'][0]),
                            float(params['range'][1]),
                            float(params['range'][2]))
            except:
                logging.critical("\nERROR: Binding site coordinates for docking are missing or incorrectly defined. Try:\n- ensuring center and range values are entered correctly\n- using a reference ligand instead")
                sys.exit()
    else:
        coords = get_coordinates(params['ref_lig'], dock_settings['padding'])

    if not os.path.isdir(os.path.join('utils','temp','pdb_files')):
        os.makedirs(os.path.join('utils','temp','pdb_files'))
        os.makedirs(os.path.join('utils','temp','pdbqt_files'))
        os.makedirs(os.path.join('utils','temp','docked_pdbqt_files'))

    pdbs = get_filepaths(os.path.join('utils','temp','pdb_files',''))
    for pdb in pdbs:
        os.remove(pdb)

    pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))
    for pdbqt in pdbqts:
        os.remove(pdbqt)

    docked_pdbqts = get_filepaths(os.path.join('utils','temp','docked_pdbqt_files',''))
    for docked_pdbqt in docked_pdbqts:
        os.remove(docked_pdbqt)

    smi_dict = get_smiles(params['ligand'])

    logging.info('Generating 3D pdbs from SMILES...')

    with tqdm_joblib(tqdm(desc="Generating...", total=len(smi_dict))) as progress_bar:
        Parallel(n_jobs=params['threads'])(delayed(make_pdbs_from_smiles)(smi) for smi in smi_dict.items())

    pdbs = os.listdir(os.path.join('utils','temp','pdb_files',''))

    logging.info('Converting pdbs to pdbqts...')

    merged_pdb_args = merge_args(os.path.join('utils','MGLTools-1.5.6',''), pdbs)

    with tqdm_joblib(tqdm(desc="Converting...", total=len(merged_pdb_args))) as progress_bar:
        Parallel(n_jobs=params['threads'])(delayed(autodock_convert)(pdb_arg) for pdb_arg in merged_pdb_args.items())

    pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))

    if sys.platform.lower() == 'darwin':
        os_name = 'mac'
    elif 'linux' in sys.platform.lower():
        os_name = 'linux'

    logging.info("Docking pdbqt ligands...")
    for pdbqt in tqdm(pdbqts):
        dock_file(
                    os.path.join('utils','gwovina-1.0','build',os_name,'release','gwovina'),
                    params['receptor'],
                    pdbqt,
                    *coords,
                    dock_settings['gwovina_settings']['exhaustiveness'],
                    dock_settings['gwovina_settings']['num_wolves'],
                    dock_settings['gwovina_settings']['num_modes'],
                    dock_settings['gwovina_settings']['energy_range'],
                    outfile=os.path.join(f'{stem_path}','utils','temp','docked_pdbqt_files',f'{os.path.split(pdbqt)[1]}')
                    )


    if '.' in params['ligand']:
        docked_ligands_folder = os.path.basename(params['ligand']).split('.')[0]
    else:
        docked_ligands_folder = os.path.basename(params['ligand'])

    docked_ligands_path = os.path.join('docked_ligands',docked_ligands_folder,'')


    params['ligand'] = [os.path.join('utils','temp','docked_pdbqt_files', file) for file in os.listdir(os.path.join('utils','temp','docked_pdbqt_files'))]
    receptors = [params['receptor'] for i in range(len(params['ligand']))]
    params['receptor'] = receptors

    if not os.path.isdir('docked_ligands'):
        os.mkdir('docked_ligands')
    if not os.path.isdir(docked_ligands_path):
        os.makedirs(docked_ligands_path)

    for file in params['ligand']:
        shutil.copy(file, docked_ligands_path)   
    
    return params, smi_dict

def parse_ligand_poses(params):

    poses = list(map(lambda x: multiple_pose_check(x, params['pose_1']), params['ligand']))

    if params['pose_1']:
        poses = [pose[0] for pose in poses]
        receptor_ligand_args = list(zip(params['receptor'], params['ligand'], poses))
    else:
        receptor_ligand_args = list(map(lambda x,y,z: product([x],[y],z),params['receptor'],params['ligand'],poses))
        receptor_ligand_args = list(chain.from_iterable(receptor_ligand_args))
    
    return receptor_ligand_args

def calculate_batches_needed(receptor_ligand_args):

    total_poses = len(receptor_ligand_args)
    estimated_ram_usage = (360540*total_poses) + 644792975
    available_ram = psutil.virtual_memory().total
    safe_ram_available = available_ram*0.8

    if estimated_ram_usage > safe_ram_available:
        batches_needed = math.ceil(estimated_ram_usage/safe_ram_available)
    else:
        batches_needed = 1

    return batches_needed

def prepare_features(receptor_ligand_args):

    filterwarnings('ignore')

    """
    Function: Wrapper to prepare            
    all requested protein-ligand            
    complexes/poses for scoring             
                                            
    Inputs: User defined params dictionary  
    (as global),                            
    dictionary of paired ligand/receptor    
    filepaths                               
                                            
    Output: Writes results as row to output 
    csv file                                
    """

    ligand_pose_number = receptor_ligand_args[2][0]
    ligand_pdbqt_block = receptor_ligand_args[2][1]

    receptor_filepath = receptor_ligand_args[0]
    ligand_filepath = receptor_ligand_args[1]
    ligand_basename = os.path.basename(ligand_filepath)
    ligand_basename = ligand_basename.replace('.pdbqt', ligand_pose_number)
    receptor_basename = os.path.basename(receptor_filepath)

    features = extract(ligand_pdbqt_block, receptor_filepath)

    multi_pose_features = prune_df_headers(features)

    multi_pose_features.fillna(0, inplace=True)

    multi_pose_features['Receptor'] = receptor_basename
    multi_pose_features['Ligand'] = ligand_basename
    
    return multi_pose_features

def scale_multipose_features(df):

    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    scaler_58 = reference_headers.get('for_scaler_58')
    headers_58 = reference_headers.get('492_models_58')

    ligands, receptors = df['Ligand'], df['Receptor']

    missing_columns = list(set(scaler_58) - set(list(df)))

    for col in missing_columns:
        df[col] = 0

    df = df[scaler_58]
    scaler = load(os.path.join('utils','params','58_maxabs_scaler_params.save'))
    scaled = scaler.transform(df)
    df[df.columns] = scaled
    df = df[headers_58]

    df['Ligand'], df['Receptor'] = ligands, receptors

    return df

def score(models):

    """
    Function: Score supplied ligands with   
    an individual model                     
                                            
    Inputs: Tuple of (model_name,           
                      model_binary_file,    
                      feature dataframes)   
                                            
    Output: Dataframe of model predictions  
    """

    model_name = models[0]

    model_file = models[1]

    features = models[2]

    logging.info(f'Scoring with {model_name}...')

    results = features[['Ligand','Receptor']].copy().reset_index(drop=True)
    df = features.drop(['Ligand','Receptor'], axis=1)

    if 'xgb' in model_name:
        dtest = xgb.DMatrix(df, feature_names=df.columns)
        results[model_name] = model_file.predict(dtest)

    else:
        network_predictions = run_networks(df, model_file, model_name)
        results[network_predictions.columns] = network_predictions[network_predictions.columns]

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

    if not params['dock']:
        logging.info(f'Found {len(params["ligand"])} ligand(s) for scoring against a single receptor...\n')

    else:
        ligand_count = open(params["ligand"]).read().split("\n")
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


    logging.info('**************************************************************************\n')
    logging.info('Model Request Summary:\n')

    models = {}

    logging.info('XGBoost Model: Yes')
    xgb_path = os.path.join('utils','models','xgboost_models','495_models_58_booster.pkl')
    models['xgboost_model'] = pickle.load(open(xgb_path,'rb'))

    logging.info('Feedforward NN Model : Yes')
    logging.info(f'- Using Best {params["num_networks"]} Networks')
    models['ff_nn'] = os.path.join('utils','models','ff_nn_models')
    model_ranks = pickle.load(open(os.path.join(models['ff_nn'],'rankings.pkl'),'rb'))
    model_ranks = model_ranks[:params["num_networks"]]
    models['ff_nn'] = [os.path.join(models['ff_nn'], 'models',f'{model[1]}.hdf5') for model in model_ranks]

    logging.info('W&D NN Model : Yes')
    logging.info(f'- Using Best {params["num_networks"]} Networks')
    models['wd_nn'] = os.path.join('utils','models','wd_nn_models')
    model_ranks = pickle.load(open(os.path.join(models['wd_nn'],'rankings.pkl'),'rb'))
    model_ranks = model_ranks[:params["num_networks"]]
    models['wd_nn'] = [os.path.join(models['wd_nn'], 'models',f'{model[1]}.hdf5') for model in model_ranks]
    logging.info('\n')

    if params['pose_1']:

        logging.info('Calculating scores for first pose only in pdbqt file(s)\n')


    logging.info('**************************************************************************\n')

    return models

def score_ligand_batch(params, ligand_batch, model_binaries):

    with tqdm_joblib(tqdm(desc="Preparing features", total=len(ligand_batch))) as progress_bar:
        multi_pose_features = Parallel(n_jobs=params['threads'])(delayed(prepare_features)(ligand) for ligand in ligand_batch)

    multi_pose_features = pd.concat(multi_pose_features)

    multi_pose_features = scale_multipose_features(multi_pose_features)

    models = [(m[0], m[1], multi_pose_features) for m in model_binaries]

    model_results = list()

    logging.info('**************************************************************************\n')

    for model in models:
        model_results.append(score(model))
        logging.info('Done!')

    logging.info('**************************************************************************\n')

    merged_results = reduce(lambda x, y: pd.merge(x, y, on = ['Receptor','Ligand']), model_results)

    multi_models = ['xgboost_model',
                    'ff_nn_models_average',
                    'wd_nn_models_average']

    merged_results['SCORCH_pose_score'] = merged_results[multi_models].mean(axis=1)
    max_std = 0.4714 # result from [0, 0, 1] or [1, 1, 0]
    minimum_val = 1 - max_std
    
    merged_results['SCORCH_stdev'] = merged_results[multi_models].std(axis=1, ddof=0)
    merged_results['SCORCH_certainty'] = ((1-merged_results['SCORCH_stdev'])-minimum_val)/max_std

    merged_results = merged_results[['Receptor',
                                        'Ligand',
                                        'SCORCH_pose_score',
                                        'SCORCH_certainty']].copy()

    merged_results['Ligand_ID'] = merged_results['Ligand'].apply(get_ligand_id)
    merged_results['Pose_Number'] = merged_results['Ligand'].apply(lambda x: x.split('_pose_')[-1])

    del merged_results['Ligand']

    return merged_results

def create_final_results(params, ligand_scores):

    final_ligand_scores = pd.concat(ligand_scores)
    final_ligand_scores['SCORCH_score'] = final_ligand_scores.groupby(['Ligand_ID'])['SCORCH_pose_score'].transform('max')
    final_ligand_scores['best_pose'] = np.where(final_ligand_scores.SCORCH_score == final_ligand_scores.SCORCH_pose_score, 1, 0)

    if not params['return_pose_scores']:
        final_ligand_scores = final_ligand_scores.loc[final_ligand_scores.best_pose == 1]
        final_ligand_scores = final_ligand_scores[['Receptor',
                                                    'Ligand_ID',
                                                    'Pose_Number',
                                                    'SCORCH_score',
                                                    'SCORCH_certainty']].copy()
        final_ligand_scores = final_ligand_scores.sort_values(by='SCORCH_score', ascending=False)

    else:
        final_ligand_scores = final_ligand_scores.sort_values(by=['SCORCH_score','SCORCH_pose_score'], ascending=False)

    numerics = list(final_ligand_scores.select_dtypes(include=[np.number]))
    final_ligand_scores[numerics] = final_ligand_scores[numerics].round(5)

    return final_ligand_scores

def scoring(params):

    """
    Function: Score protein-ligand          
    complex(es)                             
                                            
    Inputs: User command line parameters    
    dictionary                              
                                            
    Output: Dataframe of scoring function   
    predictions                             
    """

    print_intro(params)

    if params['dock']:

       params, smi_dict = prepare_and_dock_inputs(params)

    from pprint import pprint

    pprint(params)

    input('')

    receptor_ligand_args = parse_ligand_poses(params)

    batches_needed = calculate_batches_needed(receptor_ligand_args)

    all_ligands_to_score = list_to_chunks(receptor_ligand_args, batches_needed)
    model_dict = prepare_models(params)
    model_binaries = list(model_dict.items())

    ligand_scores = list()

    logging.info('**************************************************************************\n')
    print('Scoring batches...')
    for batch_number, ligand_batch in enumerate(all_ligands_to_score):

        logging.info(f"Scoring ligand batch {batch_number + 1} of {batches_needed}")
        
        merged_results = score_ligand_batch(params, ligand_batch, model_binaries)

        ligand_scores.append(merged_results)

    final_ligand_scores = create_final_results(params, ligand_scores)

    if params['dock']:
        final_ligand_scores['Ligand_SMILE'] = final_ligand_scores['Ligand_ID'].map(smi_dict)

    return final_ligand_scores

#######################################################################
# Main script

if __name__ == "__main__":

    # parse user arguments
    print('Parsing arguments...')
    params = parse_args(sys.argv)
    print('Starting scoring...')

    # score complexes
    scoring_function_results = scoring(params)

    # output results
    if not params['out']:

        # send to stdout if no outfile given
        sys.stdout.write(scoring_function_results.to_csv(index=False))

    else:

        # otherwise save to user specified csv
        scoring_function_results.to_csv(params['out'], index=False)


# scorch.py end
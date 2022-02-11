#######################################################################
# SCORCH script version 1.0 - Run python scoring.py -h for help    #
#                                                                     #
# Script Authors:                                                     #
# @sammoneykyrle                                                      #
# @milesmcgibbon                                                      #
#                                                                     #
# School of Biological Sciences                                       #
# The University of Edinburgh                                         #
#######################################################################


# import all libraries and ignore tensorflow warnings
import textwrap
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
tf.get_logger().setLevel('ERROR')
from tensorflow.keras.models import load_model
from utils import binana
from sys import platform
from utils import kier
import logging
from utils.ecifs import *
from utils.dock_functions import *
import pandas as pd
import multiprocessing as mp
import sys
import joblib
from functools import reduce
import pickle
import numpy as np
import shutil
import json
import xgboost as xgb
from tqdm import tqdm
from warnings import filterwarnings
from itertools import product, chain
from functools import partialmethod
filterwarnings('ignore')

# get working directory where scoring function is being deployed
stem_path = os.getcwd()

def get_ligand_id(ligand):
    ###########################################
    # Function: Get ligand ID from pose scores#
    #                                         #
    # Inputs: Ligand pose name                #
    #                                         #
    # Output: Ligand ID base name             #
    ###########################################

    Ligand_ID = ligand.split('_pose')[0]
    return Ligand_ID

def run_binana(params, lig, rec):

    ###########################################
    # Function: Get BINANA descriptors for    #
    # protein-ligand complex                  #
    #                                         #
    # Inputs: BINANA parameters dictionary,   #
    # ligand as a pdbqt string block,         #
    # receptor pdbqt filepath                 #
    #                                         #
    # Output: BINANA protein-ligand complex   #
    # descriptor features as a DataFrame      #
    ###########################################

    if params.okay_to_proceed() == False:
        logging.critical(
            "Error: You need to specify the ligand and receptor PDBQT files to analyze using\nthe -receptor and -ligand tags from the command line.\n"
        )
        sys.exit(0)

    if params.error != "":
        logging.warning("Warning: The following command-line parameters were not recognized:")
        logging.warning(("   " + cmd_params.error + "\n"))

    output = binana.Binana(lig, rec, params).out

    return binana.parse(output, 0)

def kier_flexibility(lig):

    ###########################################
    # Function: Calculate Kier flexibility    #
    # for ligand                              #
    #                                         #
    # Inputs: ligand as a pdbqt string block  #
    #                                         #
    # Output: Kier flexibility                #
    ###########################################

    mol = kier.SmilePrep(lig)
    return kier.CalculateFlexibility(mol)

def calculate_ecifs(lig, rec):

    ###########################################
    # Function: Get ECIFs for protein-ligand  #
    # complex                                 #
    #                                         #
    # Inputs: ligand as a pdbqt string block, #
    # receptor pdbqt filepath                 #
    #                                         #
    # Output: ECIF protein-ligand complex     #
    # descriptor features as a DataFrame      #
    ###########################################

    ECIF_data = GetECIF(rec, lig, distance_cutoff=6.0)
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]
    ECIF_data = dict(zip(ECIFHeaders,ECIF_data))
    ECIF_df = pd.DataFrame(ECIF_data,index=[0])

    return ECIF_df

def extract(params):

    ###########################################
    # Function: Get all descriptor features   #
    # for protein-ligand complex              #
    #                                         #
    # Inputs: User defined params dictionary  #
    #                                         #
    # Output: All protein-ligand complex      #
    # descriptor features as a DataFrame      #
    ###########################################

    lig = params.params["ligand"]
    rec = params.params["receptor"]
    k = kier_flexibility(lig)
    bin = run_binana(params,lig,rec)
    ECIF = calculate_ecifs(lig, rec)
    df = pd.concat([ECIF,bin],axis=1)
    df['Kier Flexibility'] = k
    return df

def transform_df(df):

    ###########################################
    # Function: Condense and scale descriptor #
    # features for model input                #
    #                                         #
    # Inputs: Full Dataframe of               #
    # protein-ligand complex descriptors,     #
    # boolean for single pose model type,     #
    # boolean for further condensation with   #
    # principle component analysis            #
    #                                         #
    # Output: DataFrame of features for model #
    # input                                   #
    ###########################################

    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    scaler_58 = reference_headers.get('for_scaler_58')
    headers_58 = reference_headers.get('492_models_58')

    df = df[scaler_58]
    scaler = joblib.load(os.path.join('utils','params','58_maxabs_scaler_params.save'))
    scaled = scaler.transform(df)
    df[df.columns] = scaled
    df = df[headers_58]

    return df

def multiple_pose_check(lig, pose_1):

    ###########################################
    # Function: Transform ligand.pdbqt        #
    # poses/models into pdbqt string blocks   #
    #                                         #
    # Inputs: ligand.pdbqt filepath           #
    #                                         #
    # Output: List of model/pose pdbqt string #
    # blocks                                  #
    ###########################################

    pdbqt_pose_blocks = list()
    lig_text = open(lig, 'r').read()
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

    ###########################################
    # Function: Get prediction from MLP model #
    # for protein-ligand complex              #
    #                                         #
    # Inputs: Number of networks as integer,  #
    # condensed protein-ligand complex        #
    # features as DataFrame                   #
    #                                         #
    # Output: Float prediction                #
    ###########################################


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

    ###########################################
    # Function: Get prediction from XGB model #
    # for protein-ligand complex              #
    #                                         #
    # Inputs: Condensed protein-ligand        #
    # complex features as DataFrame,          #
    # boolean for single pose model           #
    #                                         #
    # Output: Float prediction                #
    ###########################################

    global xgboost_models
    dtest = xgb.DMatrix(df, feature_names=df.columns)
    prediction = xgboost_models.predict(dtest)
    return prediction

def test(params):

    ###########################################
    # Function: Wrapper to score              #
    # single protein-ligand complex/pose      #
    #                                         #
    # Inputs: User defined params dictionary  #
    #                                         #
    # Output: Float single prediction or      #
    # consensus mean prediction from all      #
    # three models                            #
    ###########################################

    cmd_params = binana.CommandLineParameters(params['binana_params'].copy())
    features = extract(cmd_params)
    results = []

    if params['ff_nn'] == True:
        df = transform_df(features)
        mlp_result = run_networks(params['num_networks'],df,'ff_nn')
        results.append(mlp_result)

    if params['wd_nn'] == True:
        df = transform_df(features)
        mlp_result = run_networks(params['num_networks'],df,'wd_nn')
        results.append(mlp_result)

    if params['xgboost_model'] == True:
        df = transform_df(features)
        xgb_result = run_xgboost_models(df)
        results.append(xgb_result)

    return np.mean(results)

def binary_concat(dfs, headers):

    ###########################################
    # Function: Concatenate list of           #
    # dataframes into a single dataframe by   #
    # sequentially writing to a single binary #
    # file (removes pd.concat bottleneck)     #
    #                                         #
    # Inputs: List of dataframes, dataframe   #
    # headers as a list                       #
    #                                         #
    # Output: Single combined dataframe       #
    ###########################################

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

    ###########################################
    # Function: Parse user arguments when     #
    # script is imported as a module          #
    #                                         #
    # Inputs: User arguments as a dictionary  #
    #                                         #
    # Output: Populated params dictionary     #
    ###########################################


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

    ###########################################
    # Function: Parse user defined command    #
    # line arguments                          #
    #                                         #
    # Inputs: Command line arguments          #
    #                                         #
    # Output: Populated params dictionary     #
    ###########################################

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
                    logging.info(line)
                elif  '-' in line:
                    logging.info(line)
                else:
                    logging.info(wrapper.fill(line))
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
            logging.basicConfig(level=logging.CRITICAL, format='%(message)s')

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

def prepare_features(receptor_ligand_args):

    ###########################################
    # Function: Wrapper to prepare            #
    # all requested protein-ligand            #
    # complexes/poses for scoring             #
    #                                         #
    # Inputs: User defined params dictionary  #
    # (as global),                            #
    # dictionary of paired ligand/receptor    #
    # filepaths                               #
    #                                         #
    # Output: Writes results as row to output #
    # csv file                                #
    ###########################################

    params = receptor_ligand_args[1]
    receptor_ligand_args = receptor_ligand_args[0]

    receptor = receptor_ligand_args[0]
    ligand = receptor_ligand_args[1]

    lig_block = receptor_ligand_args[2]
    head, name = os.path.split(ligand)
    pose_number = lig_block[0]
    lig_block = lig_block[1]

    ligand_name = name.replace('.pdbqt', pose_number)
    receptor_name = os.path.split(receptor)[-1]
    params['binana_params'][1] = receptor
    params['binana_params'][3] = lig_block
    cmd_params = binana.CommandLineParameters(params['binana_params'].copy())
    features = extract(cmd_params)

    return (receptor_name, ligand_name, features)

def score(models):

    ###########################################
    # Function: Score supplied ligands with   #
    # an individual model                     #
    #                                         #
    # Inputs: Tuple of (model_name,           #
    #                   model_binary_file,    #
    #                   feature dataframes)   #
    #                                         #
    # Output: Dataframe of model predictions  #
    ###########################################

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

def multiprocess_wrapper(function, items, threads):

    ###########################################
    # Function: Parallelise scoring           #
    # protein/ligand complexes                #
    #                                         #
    # Inputs: Function to parallelise         #
    # (def score),                            #
    # list of tuples as function input,       #
    # number of threads to parallelise across #
    #                                         #
    # Output: List of returned results        #
    ###########################################

    processes = min(threads, mp.cpu_count())
    with mp.Pool(processes) as p:
        results = list(tqdm(p.imap(function, items), total=len(items)))
        p.close()
        p.join()

    return results

def print_intro(params):

    ###########################################
    # Function: Prints chosen arguments to    #
    # stdout                                  #
    #                                         #
    # Inputs: User command line parameters    #
    # dictionary                              #
    #                                         #
    # Output: None                            #
    ###########################################


    logging.info('\n')
    logging.info('**************************************************************************')

    logging.info('SCORCH v1.0')
    logging.info('Miles McGibbon, Samuel Money-Kyrle, Vincent Blay & Douglas R. Houston\n')

    logging.info('**************************************************************************\n')

    if not params['dock']:
        logging.info(f'Parsed {len(params["ligand"])} ligand(s) for scoring against a single receptor...\n')

    else:
        ligand_count = len(open(params["ligand"]).read().split("\n"))
        logging.info(f'Parsed {ligand_count} ligand smiles for docking and scoring against a single receptor...\n')

        logging.info('**************************************************************************\n')

def prepare_models(params):

    ###########################################
    # Function: Loads machine-learning model  #
    # binaries                                #
    #                                         #
    # Inputs: User command line parameters    #
    # dictionary                              #
    #                                         #
    # Output: Dictionary of {model_name:      #
    #                        model_binary}    #
    ###########################################


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

def scoring(params):

    ###########################################
    # Function: Score protein-ligand          #
    # complex(es)                             #
    #                                         #
    # Inputs: User command line parameters    #
    # dictionary                              #
    #                                         #
    # Output: Dataframe of scoring function   #
    # predictions                             #
    ###########################################

    print_intro(params)

    if len(params['ligand']) == 1:
        binana_ligand = params['ligand'][0]
        binana_receptor = params['receptor']
    else:
        binana_ligand = os.path.dirname(params['ligand'][0]) + '/'
        binana_receptor = params['receptor'][0]

    params['binana_params'] = ['-receptor', binana_receptor, '-ligand', binana_ligand]

    if params['dock']:

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

        multiprocess_wrapper(make_pdbs_from_smiles, smi_dict.items(), params['threads'])

        pdbs = os.listdir(os.path.join('utils','temp','pdb_files',''))
        logging.info('Converting pdbs to pdbqts...')
        merged_pdb_args = merge_args(os.path.join('utils','MGLTools-1.5.6',''), pdbs)

        multiprocess_wrapper(autodock_convert, merged_pdb_args.items(), params['threads'])

        pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))

        if platform.lower() == 'darwin':
            os_name = 'mac'
        elif 'linux' in platform.lower():
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


    poses = list(map(lambda x: multiple_pose_check(x, params['pose_1']), params['ligand']))

    if params['pose_1']:
        poses = [pose[0] for pose in poses]
        receptor_ligand_args = list(zip(params['receptor'], params['ligand'], poses))
    else:
        receptor_ligand_args = list(map(lambda x,y,z: product([x],[y],z),params['receptor'],params['ligand'],poses))
        receptor_ligand_args = list(chain.from_iterable(receptor_ligand_args))

    receptor_ligand_args  = list(map(lambda x: (x, params), receptor_ligand_args))
    features = multiprocess_wrapper(prepare_features, receptor_ligand_args, params['threads'])
    feature_headers = list(features[0][2])
    features_df = binary_concat([i[2] for i in features], feature_headers)
    multi_pose_features = transform_df(features_df)
    multi_pose_features['Receptor'] = [i[0] for i in features]
    multi_pose_features['Ligand'] = [i[1] for i in features]

    models = prepare_models(params)
    models = list(models.items())
    models = [(m[0], m[1], multi_pose_features) for m in models]

    model_results = list()

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
    maxmum_val = 1
    merged_results['SCORCH_stdev'] = merged_results[multi_models].std(axis=1, ddof=0)
    merged_results['SCORCH_certainty'] = ((1-merged_results['SCORCH_stdev'])-minimum_val)/max_std

    merged_results = merged_results[['Receptor',
                                     'Ligand',
                                     'SCORCH_pose_score',
                                     'SCORCH_certainty']].copy()

    merged_results['Ligand_ID'] = merged_results['Ligand'].apply(get_ligand_id)
    merged_results['SCORCH_score'] = merged_results.groupby(['Ligand_ID'])['SCORCH_pose_score'].transform('max')
    merged_results['best_pose'] = np.where(merged_results.SCORCH_score == merged_results.SCORCH_pose_score, 1, 0)

    if not params['return_pose_scores']:
        merged_results = merged_results.loc[merged_results.best_pose == 1]
        merged_results = merged_results[['Receptor',
                                         'Ligand',
                                         'SCORCH_score',
                                         'SCORCH_certainty']].copy()

    numerics = list(merged_results.select_dtypes(include=[np.number]))
    merged_results[numerics] = merged_results[numerics].round(5)
    return merged_results

if __name__ == "__main__":

    params = parse_args(sys.argv)
    scoring_function_results = scoring(params)
    if not params['out']:
        sys.stdout.write(scoring_function_results.to_csv(index=False))
    else:
        scoring_function_results.to_csv(params['out'], index=False)

<<<<<<< HEAD
-return_pose_scores# SCORCH
=======
# SCORCH
>>>>>>> 2a86f1dd07e8b2b9c818443d97be01e659407d61

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.6.9-blue)](https://www.python.org/downloads/release/python-369/)
![](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

SCORCH is a fast scoring function based on a consensus of machine learning models. Scoring functions are used to evaluate poses of molecules obtained from molecular docking. SCORCH scores range from 0 to 1, with higher values indicating a higher probability of the molecule binding tightly to the receptor.

SCORCH uses `.pdbqt` files as input for the scoring, which is the format used by [Autodock](https://autodock.scripps.edu/), [Vina](https://dx.doi.org/10.1002/jcc.21334), and [GWOVina](https://cbbio.online/software/gwovina/index.html) docking software, among others. Additionally, this release contains an integrated pipeline to dock and score molecules in SMILES format using GWOVina.

SCORCH uses a variety of descriptors to characterize a docked pose, including [Binana 1.3](https://git.durrantlab.pitt.edu/jdurrant/binana/-/tree/1.3) and [ECIFs](https://github.com/DIFACQUIM/ECIF). The contributing models were trained on multiple docked poses for each ligand, labelled based on their RMSD to crystal structures. SCORCH has used over 54,000 poses in its training. As a result, SCORCH avoids biases and provides improved accuracy to identify true binder molecules in virtual screening. Read more in our [publication]().


![](https://raw.githubusercontent.com/miles-mcgibbon/miles-mcgibbon/main/.github/images/pose_labels.gif)

---

# Installation

Installation on linux and mac is achieved via [virtualenv](https://virtualenv.pypa.io/en/latest/). The installation of virtualenv, a local build of Python 3.6.9, scoring function dependencies and scoring function setup is all performed with the supplied setup bash script.

To install SCORCH:

```bash
# clone the GitHub repository
git clone https://github.com/miles-mcgibbon/ML-SCORCH.git

# ensure setup.sh is executable
cd ML-SCORCH
chmod 755 setup.sh

# execute the setup script
./setup.sh
```

# Receptor & Ligand Preparation

The scoring function accepts `.pdbqt` receptor and ligand files as inputs. These should be prepared with [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/) using the `prepare_receptor4.py` and `prepare_ligand4.py` Python scripts as follows:

```bash
# preparing a ligand
~/utils/MGLTools-1.5.6/bin/pythonsh \
~/utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l path/to/ligand.pdb \
-A hydrogens \
-o path/to/save/ligand.pdbqt
-U nphs

# preparing a receptor
~/utils/MGLTools-1.5.6/bin/pythonsh \
~/utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
-l path/to/receptor.pdb \
-A hydrogens \
-o path/to/save/receptor.pdbqt
-U nphs
```

`.sdf` and `.mol2` ligands can be converted to `.pdb` format using `rdkit` in Python:

```python
from rdkit import Chem

# for mol2 files
input_mol = Chem.MolFromMol2File("ligand.mol2")

# for sdf files
input_mol = Chem.SDMolSupplier("ligand.sdf")[0]

# save the output as pdb
Chem.MolToPDBFile(input_mol, "ligand.pdb")
```

# Usage

To use the scoring function, the virtual environment needs to be activated first:

```bash
source .scoring/bin/activate
```

<<<<<<< HEAD
The scoring function is supplied as the Python script `scorch.py`. Its main arguments are:
=======
The scoring function is supplied as the Python script `scoring.py`. Its main arguments are:
>>>>>>> 2a86f1dd07e8b2b9c818443d97be01e659407d61

|Argument     |Value                                                                                     |Importance                  |
|-------------|------------------------------------------------------------------------------------------|----------------------------|
|`-receptor`    |Filepath to receptor file (pdbqt)                                                       |Essential                   |
|`-ligand `     |Filepath to ligand(s) (pdbqt or SMILES)                                                 |Essential                   |
|`-ref_lig`    |Filepath to example ligand in receptor binding site (pdb or pdbqt)                       |Essential for SMILES ligands|
|`-out`         |Filepath for output csv (If not supplied, scores are written to stdout)                  |Optional (Default stdout)   |
<<<<<<< HEAD
|`-return_pose_scores` |If supplied, scoring values for individual poses in each ligand file are returned | Optional (Default False) |
=======
|`-return_scores_poses` |If supplied, scoring values for individual poses in each ligand file are returned | Optional (Default False) |
>>>>>>> 2a86f1dd07e8b2b9c818443d97be01e659407d61
|`-threads`     |Number of threads to use                                                                |Optional (Default 1)        |
|`-verbose`     |If supplied, progress bars and indicators displayed while scoring                  |Optional (Default False)    |

Additional options are explained in the function help.

### Scoring a Single PDBQT Ligand Against a Single PDBQT Receptor

For scoring a single ligand - `/home/user/ligands/ligand.pdbqt` - against a single receptor - `/home/user/receptors/receptor.pdbqt`:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt
```

### Scoring Multiple PDBQT Ligands Against a Single PDBQT Receptor

For scoring all `.pdbqt` ligands in the directory - `/home/user/ligands/` - against a single receptor - `/home/user/receptors/receptor.pdbqt` - just supply the directory path to the `-ligand` argument:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/  \
-out scoring_results.csv
# this writes the output scores to a file
-out scoring_results.csv \
# this parallelises the scoring over 6 threads
-threads 6 \
# this displays scoring progress
-verbose
```

The `-verbose` flag is here used in conjunction with the `-out` flag, otherwise the progress indicators will be written to the results file.


### Docking and Scoring Multiple SMILES Ligands Against a Single PDBQT Receptor


The module also includes a full pipeline to convert SMILES ligands to `.pdbqt` files using MGLTools 1.5.6, dock them using GWOVina, and score them with SCORCH:

```bash
python scorch.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/smiles/ligands.smi  \
-ref_lig /home/user/ligands/reference_ligand.pdb \
-out influenza_results.csv
```

The `-ligand` input should be supplied as `.smi`  or text file, with one SMILES ligand and an optional identifier per line, as in this example:

```bash
CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(=O)C)O)C(=O)OCC 49817880
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)[NH3+])C(=O)OCC 24848267
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)NOC 118722031
```

You also need to supply a `.pdb` or `.pdbqt` ligand on the active site of the protein receptor that you are screening using the `-ref_lig` argument. This is used to define the docking search space.

Docking settings can be changed by editing the `utils/params/dock_settings.json` file in the scoring function folder. See the function help for additional details.

### Importing SCORCH as a Python module

<<<<<<< HEAD
The main function from `scorch.py` can be imported and used in other Python scripts. It takes a dictionary of parameters as inputs and returns a pandas dataframe of model scores identical to the normal scoring function output:
=======
The main function from `scoring.py` can be imported and used in other Python scripts. It takes a dictionary of parameters as inputs and returns a pandas dataframe of model scores identical to the normal scoring function output:
>>>>>>> 2a86f1dd07e8b2b9c818443d97be01e659407d61

```python
from scorch import scoring

input_parameters = {'binana_params': ['-receptor',
                                       'path/to/receptor.pdbqt',
                                       '-ligand',
                                       'path/to/ligand.pdbqt'],
                  'dock': False,
<<<<<<< HEAD
		  '-return_pose_scores': False,
=======
		  '-return_poses_scores': False,
>>>>>>> 2a86f1dd07e8b2b9c818443d97be01e659407d61
                  'ligand': ['path/to/ligand.pdbqt'],
                  'out': 'output.csv',
                  'receptor': ['path/to/receptor.pdbqt'],
                  'ref_lig': None,
                  'threads': 4,
                  'verbose': False}

output = scoring(input_parameters)
```

# Output

Scores are output in `.csv` format. For example, scoring a single ligand pdbqt containing 10 docked poses against a single receptor file would yield the following output. Note that the output of SCORCH also includes a measure of the prediction certainty.

|Receptor      |Ligand        |SCORCH_score  |SCORCH_certainty|
|--------------|--------------|--------------|----------------|
<<<<<<< HEAD
|receptor.pdbqt|ligand_pose_1      |0.83521       |0.75148         |
|receptor.pdbqt|ligand_pose_2      |0.83781       |0.75926         |
|receptor.pdbqt|ligand_pose_3      |0.84241       |0.74976         |
|receptor.pdbqt|ligand_pose_4      |0.77493       |0.65857         |
|receptor.pdbqt|ligand_pose_5      |0.72339       |0.69700         |
|receptor.pdbqt|ligand_pose_6      |0.86471       |0.78335         |
|receptor.pdbqt|ligand_pose_7      |0.81688       |0.69923         |
|receptor.pdbqt|ligand_pose_8      |0.86669       |0.78808         |
|receptor.pdbqt|ligand_pose_9      |0.65023       |0.74499         |
|receptor.pdbqt|ligand_pose_10     |0.07947       |0.86123         |
=======
|receptor.pdbqt|ligand_1      |0.83521       |0.75148         |
|receptor.pdbqt|ligand_2      |0.83781       |0.75926         |
|receptor.pdbqt|ligand_3      |0.84241       |0.74976         |
|receptor.pdbqt|ligand_4      |0.77493       |0.65857         |
|receptor.pdbqt|ligand_5      |0.72339       |0.69700         |
|receptor.pdbqt|ligand_6      |0.86471       |0.78335         |
|receptor.pdbqt|ligand_7      |0.81688       |0.69923         |
|receptor.pdbqt|ligand_8      |0.86669       |0.78808         |
|receptor.pdbqt|ligand_9      |0.65023       |0.74499         |
|receptor.pdbqt|ligand_10     |0.07947       |0.86123         |
>>>>>>> 2a86f1dd07e8b2b9c818443d97be01e659407d61

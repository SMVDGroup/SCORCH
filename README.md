# ML-SCORCH

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.6.9-blue)](https://www.python.org/downloads/release/python-369/)
![](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

ML-SCORCH is a scoring function based on a consensus of machine learning models. Scoring functions are used to evaluate poses of molecules obtained from molecular docking. ML-SCORCH scores range from 0 to 1, with higher values indicating a higher probability of the molecule binding tightly to the receptor.

ML-SCORCH uses `.pdbqt` files as input for the scoring, which is the format used by [Autodock](https://autodock.scripps.edu/), [Vina](https://dx.doi.org/10.1002/jcc.21334), and [GWOVina](https://cbbio.online/software/gwovina/index.html) docking software, among others. Additionally, this release contains an integrated pipeline to dock and score molecules in SMILES format using GWOVina.

ML-SCORCH uses a variety of descriptors to characterize a docked pose, including [Binana 1.3](https://git.durrantlab.pitt.edu/jdurrant/binana/-/tree/1.3) and [ECIFs](https://github.com/DIFACQUIM/ECIF). The contributing models were trained on multiple docked poses for each ligand, labelled based on their RMSD to crystal structures. ML-SCORCH has used over 54,000 poses in its training. As a result, ML-SCORCH avoids biases and provides improved accuracy to identify true binder molecules in virtual screening. Read more in our [publication]().


![](https://raw.githubusercontent.com/miles-mcgibbon/miles-mcgibbon/main/.github/images/pose_labels.gif)

---

# Installation

Installation on linux and mac is achieved via [virtualenv](https://virtualenv.pypa.io/en/latest/). The installation of virtualenv, a local build of python 3.6.9, scoring function dependencies and scoring function setup is all performed with the supplied setup bash script.

To install ML-SCORCH:

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

The scoring function only accepts `.pdbqt` receptor and ligand files as inputs. These should be prepared with [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/) using the `prepare_receptor4.py` and `prepare_ligand4.py` python scripts as follows:

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

SDF and mol2 ligands can be converted to pdb format using `RDKit` in python:

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

The scoring function is supplied as the python script `scoring.py`. Its arguments are:

|Argument     |Value                                                                                     |Importance                  |
|-------------|------------------------------------------------------------------------------------------|----------------------------|
|`-receptor`    |Filepath to receptor file (pdbqt)                                                         |Essential                   |
|`-ligand `     |Filepath to ligand(s) (pdbqt or SMILES)                                                   |Essential                   |
|`-ref_lig`    |Filepath to example ligand in receptor binding site (pdb or pdbqt)                        |Essential for SMILES ligands|
|`-out`         |Filepath for output csv (If not supplied scores are written to stdout)                    |Optional (Default stdout)   |
|`-threads`     |Number of threads to use                                                                  |Optional (Default 1)        |
|`-num_networks`|Number of networks to average for mlpscore and widedeepscore model scores                 |Optional (Default 15)       |
|`-verbose`     |True if supplied, progress bars and indicators displayed while scoring                    |Optional (Default False)    |
|`-pose_1`      |True if supplied, only the first model in pdbqt ligand files is scored                    |Optional (Default False)    |
|`-detailed`    |True if supplied, outputs include scores from individual mlpscore and widedeepscore models|Optional (Default False)    |


### Scoring a Single PDBQT Ligand Against a Single PDBQT Receptor

For scoring a single ligand - `/home/user/ligands/ligand.pdbqt` - against a single receptor - `/home/user/receptors/receptor.pdbqt`:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt
```

Scoring can be sped up and monitored by using some optional flags:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt  \
# this writes the output scores to a file
-out scoring_results.csv \
# this parallelises the scoring over 6 threads
-threads 6 \
# this displays scoring progress
-verbose
```

Note - the `-verbose` flag must be used in conjunction with the `-out` flag, otherwise the progress indicators will be written to the results file.

### Scoring Multiple PDBQT Ligands Against a Single PDBQT Receptor

For scoring all ligands in the directory - `/home/user/ligands/` - against a single receptor - `/home/user/receptors/receptor.pdbqt` - just supply the directory path to the `-ligand` argument:


```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/  \
-out scoring_results.csv
```

### Docking and Scoring Multiple SMILES Ligands Against a Single PDBQT Receptor


The scoring function also includes a full pipeline to convert SMILES ligands to 3D pdbqt files using [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/), dock them using [GWOVina](https://doi.org/10.1111/cbdd.13764), and score them with ML-SCORCH.

SMILES inputs should be supplied as .smi or .txt files, with one SMILES ligand and optional SMILES identifier per line separated by a space. For example, you may have a `.smi` or `.txt` file of ligands which contains the following:

```bash
CCC(CC)O[C@@H]1C=C(C[C@H]([C@H]1NC(=O)C)O)C(=O)OCC 49817880
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)[NH3+])C(=O)OCC 24848267
CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)NOC 118722031
```

Alongside an .smi or .txt file, you need to supply a pdb or pdbqt ligand bound to the active site of the protein receptor file you are screening against as a reference for GWOVina docking using the `-ref_lig` argument:

```bash
python scoring.py \
-receptor /home/user/receptors/influenza_neuraminidase.pdbqt \
-ligand /home/user/smiles/influenza_inhibitors.smi  \
-ref_lig /home/user/ligands/oseltamivir.pdb
-out influenza_results.csv
```

GWOVina docking settings can be changed by editing the `utils/params/dock_settings.json` file in the scoring function folder. The additional `padding` variable is a value in angstroms which is added to the center coordinate of the reference ligand to determine the docking site limits.

### Importing the scoring function as a python module

The main function from `scoring.py` can be imported and used in other python scripts. It takes a dictionary of parameters as inputs and returns a pandas dataframe of model scores identical to the normal scoring function output:

```python
from scoring import scoring

input_parameters = {'binana_params': ['-receptor',
                                       'path/to/receptor.pdbqt',
                                       '-ligand',
                                       'path/to/ligand.pdbqt'],
                  'concise': True,
                  'dock': False,
                  'ligand': ['path/to/ligand.pdbqt'],
                  'num_networks': 15,
                  'out': 'output.csv',
                  'pose_1': False,
                  'receptor': ['path/to/receptor.pdbqt'],
                  'ref_lig': None,
                  'screen': False,
                  'single': True,
                  'threads': 4,
                  'verbose': False,
                  'wd_nn': True,
                  'ff_nn': True,
                  'xgboost_model': True}

output = scoring(input_parameters)
```

# Output

Scores are output in `.csv` format. For example, scoring a single ligand pdbqt containing 10 docked poses against a single receptor file would yield the following output. Note that the output of ML-SCORCH also includes a measure of the prediction certainty.

|Receptor      |Ligand        |xgboost_model|ff_nn_models_average|wd_nn_models_average|model_consensus|model_certainty|
|--------------|--------------|-------------|--------------------|--------------------|---------------|---------------|
|receptor.pdbqt|ligand_pose_1 |0.9974552    |0.783138            |0.7250332           |0.83520883     |0.7514803      |
|receptor.pdbqt|ligand_pose_2 |0.99668676   |0.77808344          |0.73867756          |0.83781594     |0.7592604      |
|receptor.pdbqt|ligand_pose_3 |0.99672925   |0.82014656          |0.7103614           |0.8424124      |0.74975705     |
|receptor.pdbqt|ligand_pose_4 |0.9974542    |0.7051348           |0.6221999           |0.7749297      |0.6585699      |
|receptor.pdbqt|ligand_pose_5 |0.92504144   |0.6327202           |0.6123963           |0.723386       |0.69700235     |
|receptor.pdbqt|ligand_pose_6 |0.9980556    |0.8461044           |0.74997765          |0.86471254     |0.7833506      |
|receptor.pdbqt|ligand_pose_7 |0.9973322    |0.80237544          |0.65094346          |0.81688374     |0.69922733     |
|receptor.pdbqt|ligand_pose_8 |0.9990158    |0.8433915           |0.757661            |0.8666894      |0.78807735     |
|receptor.pdbqt|ligand_pose_9 |0.48365158   |0.7629251           |0.70412195          |0.65023285     |0.74498904     |
|receptor.pdbqt|ligand_pose_10|0.004160531  |0.07061309          |0.1636547           |0.07947611     |0.8612344      |

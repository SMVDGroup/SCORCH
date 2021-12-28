# MLScore

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.6.9-blue)](https://www.python.org/downloads/release/python-369/)
[![OS](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)]()
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)](https://doi.org/10.1111/cbdd.13764)

---

MLScore is a scoring function based on a consensus of machine learning models. The contributing models were trained on multiple docked ligand poses labelled by RMSD, meaning only well docked strong binders are rewarded, and poorly docked strong binders and all weak binders are punished.

---

# Installation

Installation on linux and mac is achieved via [virtualenv](https://virtualenv.pypa.io/en/latest/). The installation of virtualenv, a local build of python 3.6.9, scoring function dependencies and scoring function setup is all performed with the supplied setup bash script.

To install MLScore:

```bash
# clone the GitHub repository
git clone https://github.com/miles-mcgibbon/MLSCORE.git

# ensure setup.sh is executable
cd MLSCORE
chmod 755 setup.sh

# execute the setup script
./setup.sh
```

# Usage

To use the scoring function, the virtual environment needs to be activated first:

```bash
source .scoring/bin/activate
```


## Overview

The scoring function itself is supplied as a python script `scoring.py`. Arguments are:

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

## Examples

#### Scoring a Single PDBQT Ligand Against a Single PDBQT Receptor

For scoring a single ligand - `/home/user/ligands/ligand.pdbqt` - against a single receptor - `/home/user/receptors/receptor.pdbqt`:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt
```

This would print the results to stdout; they could instead be redirected to a file:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt > scoring_results.csv
```

The `-out` argument could be used to achieve exactly the same result:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt  \
-out scoring_results.csv
```

Scoring could also be sped up and monitored by using some optional flags:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/ligand.pdbqt  \
-out scoring_results.csv \
# this parallelises the scoring over 6 threads
-threads 6 \
# this displays scoring progress
-verbose
```

#### Scoring Multiple PDBQT Ligands Against a Single PDBQT Receptor

For scoring all ligands in the directory - `/home/user/ligands/` - against a single receptor - `/home/user/receptors/receptor.pdbqt` - just supply the directory filepath to the `-ligand` argument:

```bash
python scoring.py \
-receptor /home/user/receptors/receptor.pdbqt \
-ligand /home/user/ligands/  \
-out scoring_results.csv
```

#### Scoring Multiple PDBQT Ligands Against Multiple PDBQT Receptors

The scoring function can also take a folder of folders as an argument to score multiple ligands against multiple receptors. For example, say you have a benchmarking dataset stored in the following way in the folder `dataset`:

```
dataset
└───TargetA
│   │   TargetA_receptor.pdbqt
│   │   TargetA_ligand.pdbqt
│
└───TargetB
    │   TargetB_receptor.pdbqt
    │   TargetB_ligand.pdbqt
```

As long as, in each subfolder, there is **one receptor file with the suffix receptor.pdbqt** and **one ligand file with the suffix ligand.pdbqt**, the scoring function can return scores for `TargetA_ligand.pdbqt` against `TargetA_receptor.pdbqt` and `TargetB_ligand.pdbqt` against `TargetB_receptor.pdbqt`. To do this, just supply the parent directory to both the `-ligand` and `-receptor` arguments:

```bash
python scoring.py \
-receptor /home/user/dataset/ \
-ligand /home/user/dataset/  \
-out scoring_results.csv
```

**Note - the `-verbose` flag must be used in conjunction with the `-out` flag, otherwise the progress indicators will be written to the results file**

#### Docking and Scoring Multiple SMILES Ligands Against a Single PDBQT Receptor

The scoring function also includes a full pipeline to convert SMILES ligands to 3D pdbqt files using [MGLTools 1.5.6](https://ccsb.scripps.edu/mgltools/1-5-6/), dock them using [GWOVina](https://doi.org/10.1111/cbdd.13764), and score them with MLScore.

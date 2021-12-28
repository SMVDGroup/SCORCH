# MLScore

[![License](https://img.shields.io/badge/License-MIT-blue)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.6.9-blue)](https://www.python.org/downloads/release/python-369/)
![OS](https://img.shields.io/badge/OS-linux%20%7C%20OS%20X-blueviolet)
[![Preprocessing](https://img.shields.io/badge/preprocessing-MGLTools%201.5.6-brightgreen)](https://ccsb.scripps.edu/mgltools/1-5-6/)]
[![Docking](https://img.shields.io/badge/docking-GWOVina%201.0-brightgreen)]( https://doi.org/10.1111/cbdd.13764)]
---

MLScore is a scoring function based on a consensus of machine learning models. The contributing models were trained on multiple docked ligand poses labelled by RMSD, meaning only well docked strong binders are rewarded, and poorly docked strong binders and all weak binders are punished.

---

# Installation

Installation on linux and mac is achieved via [virtualenv](https://virtualenv.pypa.io/en/latest/). The installation of virtualenv, a local build of python 3.6.9, scoring function dependencies and scoring function setup is all performed with the supplied setup bash script.

To install MLScore:

```bash
# clone the GitHub repository
git clone https://github.com/miles-mcgibbon/MLSCORE.git

# ensure setup.sh is executable and execute
cd MLSCORE
chmod 755 setup.sh
./setup.sh
```

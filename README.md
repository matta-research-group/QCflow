![logo](https://github.com/matta-research-group/QCflow/blob/main/QCflow_logo_narrow.jpg?raw=true)

# QCflow

[![QCflow stable](https://github.com/matta-research-group/QCflow/actions/workflows/run_test.yml/badge.svg?branch=qcflow-0.2)](https://github.com/matta-research-group/QCflow/actions/workflows/run_test.yml)

A cheminformatics -> quantum chemistry workflow toolkit leveraging [rdkit](https://github.com/rdkit/rdkit) and [cclib](https://github.com/cclib/cclib).

A typical workflow involves:

1. Generating a series of input molecules as `SMILES` - optionally, by combining different fragments into larger molecules/oligomers  

2. For each molecule:
    a. write `.com` Gaussian input files, `.sh` SLURM submission files 
    b. submit job (assuming you are working within a HPC)
    c. parse output file to submit further calculations or retrieve descriptors

3. Combine descriptors and results in a `pandas` dataframe format or similar for plotting / further analysis

## Supported QC codes 

Only Gaussian is supported at the moment, but we plan to add support for Psi4.


## Installation

```bash
git clone https://github.com/matta-research-group/QCflow/qcflow-0.2.git
cd QCflow
# install requirements into new environment
conda env create -f QCflow.yml
conda activate QCflow
# install the QCflow package
pip install .
```

## Usage Examples

The `run_calc` and `run_torsion` functions contain example workflows that submit gaussian calculations to the KCL CREATE HPC. 
Users external to Kings College London will need to alter the slurm.py file to match their HPC submission requirements.

## Calculation settings

QCflow can prepare and submit input files for the following jobs: 
- Single point calculation, neutral -> `sp`
- Single point calculation, anion → `sp_a`
- Single point calculation, cation → `sp_c`
- Single point calculation neutral charge, cationic geometry → `n_c_geo`
- Single point calculation neutral charge, anioinc geometry → `n_a_geo`
- Geometry optimisation, neutral -> `opt`
- Torsional scan, neutral → `tor`
- Optimisation anion → `opt_a`
- Optimisation cation → `opt_c`
- Optimisation neutral + Population analysis → `pop_opt_n`
- Single point Hirshfeld calculation → `sp_hirsh`


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

## Files

```bash
├── build
│   ├── bdist.macosx-10.13-x86_64
│   └── lib
│       └── QCflow
│           ├── energy_calculations.py
│           ├── find_torsion.py
│           ├── fragments.py
│           ├── __init__.py
│           ├── load_gaussian.py
│           ├── run_gaussian.py
│           ├── slurm.py
│           ├── testing_data.py
│           ├── torsion_parser.py
│           └── write_gaussian.py
├── LICENSE
├── QCflow
│   ├── energy_calculations.py
│   ├── find_torsion.py
│   ├── fragments.py
│   ├── future_functions
│   │   └── orbital_parse.py
│   ├── __init__.py
│   ├── load_gaussian.py
│   ├── __pycache__
│   │   ├── find_torsion.cpython-311.pyc
│   │   ├── find_torsion.cpython-39.pyc
│   │   ├── fragments.cpython-311.pyc
│   │   ├── fragments.cpython-38.pyc
│   │   ├── fragments.cpython-39.pyc
│   │   ├── __init__.cpython-311.pyc
│   │   ├── __init__.cpython-38.pyc
│   │   ├── __init__.cpython-39.pyc
│   │   ├── load_gaussian.cpython-311.pyc
│   │   ├── load_gaussian.cpython-39.pyc
│   │   ├── run_opt_neutral.cpython-311.pyc
│   │   ├── run_opt_neutral.cpython-39.pyc
│   │   ├── run_opt_set_dihedral.cpython-311.pyc
│   │   ├── slurm.cpython-311.pyc
│   │   ├── slurm.cpython-39.pyc
│   │   ├── torsion_parser.cpython-311.pyc
│   │   ├── torsion_parser.cpython-39.pyc
│   │   ├── torsion_run.cpython-311.pyc
│   │   ├── torsion_run.cpython-39.pyc
│   │   ├── write_gaussian.cpython-311.pyc
│   │   └── write_gaussian.cpython-39.pyc
│   ├── run_gaussian.py
│   ├── slurm.py
│   ├── testing_data.py
│   ├── torsion_parser.py
│   └── write_gaussian.py
├── qcflow.egg-info
│   ├── dependency_links.txt
│   ├── PKG-INFO
│   ├── SOURCES.txt
│   └── top_level.txt
├── QCflow_logo_narrow.jpg
├── QCflow.yml
├── README.md
├── setup.py
└── tests
    ├── b_18_v2
    │   ├── b_18_v2_n_a_geo.log
    │   ├── b_18_v2_n_c_geo.log
    │   ├── b_18_v2_opt_a.log
    │   ├── b_18_v2_opt_c.log
    │   ├── b_18_v2_opt.com
    │   ├── b_18_v2_opt.log
    │   ├── b_18_v2_opt.sh
    │   ├── b_18_v2_sp_a.log
    │   ├── b_18_v2_sp_c.log
    │   ├── b_18_v2_sp_hirsh.log
    │   └── b_18_v2_tor.log
    ├── example_dic.json
    ├── fake_dict.json
    ├── __pycache__
    │   ├── tests.cpython-312.pyc
    │   └── tests.cpython-312-pytest-7.4.4.pyc
    ├── test_dict.json
    ├── test.py
    └── torsion_1
        └── torsion_1_tor.log
```
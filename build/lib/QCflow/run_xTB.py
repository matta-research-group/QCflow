from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
from QCflow.run_gaussian import *
from QCflow.write_psi4 import *
from QCflow.write_xTB import *
import os
import json
import subprocess
import tempfile

def run_xTB(job_name, mol_name, mol_smile, time=2, cpus=10, functional='b3lyp', basis_set='6-31g*'):
    """
    Submits a xTB calculation to CREATE HPC.

    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'opt': Optimisation neutral
        - 'cation': Geometry optimisation cation (opt_c) and single of neutral charge, cation geometry (n_c_geo)
        - 'anion': Geometry optimisation anion (opt_a) and single of neutral charge, anion geometry (n_a_geo)
        - 'sp_c': Single point calculation of neutral geometry at cation charge
        - 'sp_a': Single point calculation of neutral geometry at anion charge
    mol_name : str
        Name of the oligomer.
    mol_smile : str
        SMILES string of the oligomer.
    time (int, optional): The time limit for the job in hours. Default is 24. (Max is 48)
    cpus (int, optional): The number of CPUs to allocate for the job. Default is 10.
    functional : str, optional
        Quantum chemistry functional to be used (default is 'B3LYP').
    basis_set : str, optional
        Basis set to be used (default is '6-31G*').

    
    Returns
    -------
    - For 'opt':
        - Converts the SMILES string to an RDKit molecule object.
        - Generates 3D coordinates for the molecule.
        - Predicts the conformer geometry.
        - Runs xTB optimisation using the generated coordinates.
        - Writes a psi4 input file with the optimised geometry.
        - Runs a single point calculation.
    Finally, the function writes a SLURM script, submits the job, and returns to the previous directory.
    """
    if os.path.exists(f'{mol_name}'):
        os.chdir(f'{mol_name}') #goes into directory
    else:
        os.mkdir(f'{mol_name}') #makes a directory for the molecule
        os.chdir(f'{mol_name}') #goes into directory

    if (job_name=='opt'):
        #turns smiles string into rdkit object
        mol = Chem.MolFromSmiles(mol_smile)
        #gets rdkit estimated coordinates of dimer
        mol3d = embed_molecule(mol)

        conf_geo = rdkit_predict_conf(mol_smile)
        #writes a guassian input file
        write_xTB_psi4(job_name, mol_name, mol_smile, functional, basis_set, mol=mol3d, conformer=conf_geo)
    
    #for reorganisation calcultions
    if (job_name=='cation') or (job_name=='anion') or (job_name=='sp_c') or (job_name=='sp_a'):

        write_xTB_psi4_reorg(job_name, mol_name, functional, basis_set)

    
    #writes the slurm file
    write_slurm_psi4(job_name, mol_name, time, cpus)
    #submits the slurm jon
    submit_slurm_job(job_name, mol_name)
    #goes back to previous dirctory
    os.chdir(os.path.dirname(os.getcwd()))

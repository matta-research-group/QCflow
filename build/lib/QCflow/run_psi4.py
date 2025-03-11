from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
from QCflow.run_gaussian import *
from QCflow.write_psi4 import *
import os
import json

def run_psi4(job_name, mol_name, mol_smile, time=24, cpus=10, functional='b3lyp', basis_set='6-31g*'):
    """
    Submits a psi4 calculation to CREATE HPC.

    Parameters
    ----------
    job_name (str): The type of job to run. Possible values:
        - 'sp': Single Point neutral
        - 'opt': Optimisation neutral
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
    - For 'opt' and 'sp':
        - Converts the SMILES string to an RDKit molecule object.
        - Generates 3D coordinates for the molecule.
        - Predicts the conformer geometry.
        - Writes a psi4 input file with the conformer geometry.
    Finally, the function writes a SLURM script, submits the job, and returns to the previous directory.
    """
    if os.path.exists(f'{mol_name}'):
        os.chdir(f'{mol_name}') #goes into directory
    else:
        os.mkdir(f'{mol_name}') #makes a directory for the molecule
        os.chdir(f'{mol_name}') #goes into directory

    if (job_name=='opt') or (job_name=='sp'):
        #turns smiles string into rdkit object
        mol = Chem.MolFromSmiles(mol_smile)
        #gets rdkit estimated coordinates of dimer
        mol3d = embed_molecule(mol)

        conf_geo = rdkit_predict_conf(mol_smile)
        #writes a guassian input file
        write_psi4(job_name, mol_name, mol_smile, functional, basis_set, mol=mol3d, conformer=conf_geo)

    
    #writes the slurm file
    write_slurm_psi4(job_name, mol_name, time, cpus)
    #submits the slurm jon
    submit_slurm_job(job_name, mol_name)
    #goes back to previous dirctory
    os.chdir(os.path.dirname(os.getcwd()))
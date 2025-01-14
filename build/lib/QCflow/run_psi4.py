from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
from QCflow.run_gaussian import *
from QCflow.write_psi4 import *
import os
import json

def run_psi4(job_name, mol_name, mol_smile, functional='B3LYP', basis_set='6-31G*'):
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
    functional : str, optional
        Quantum chemistry functional to be used (default is 'B3LYP').
    basis_set : str, optional
        Basis set to be used (default is '6-31G*').

    
    Returns
    -------
    - For 'sp_a', 'sp_c', 'opt_a', 'opt_c', 'n_a_geo', 'n_c_geo', and 'sp_hirsh':
        - Writes a Gaussian input file directly.
    - For 'opt', 'sp', and 'pop_opt_n':
        - Converts the SMILES string to an RDKit molecule object.
        - Generates 3D coordinates for the molecule.
        - Predicts the conformer geometry.
        - Writes a Gaussian input file with the conformer geometry.
    - For 'tor':
        - Converts the SMILES string to an RDKit molecule object.
        - Identifies the bond for torsional scan.
        - Determines the torsion angle.
        - Generates 3D coordinates for the molecule.
        - Writes a Gaussian input file with the torsion angle.
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
    write_slurm_psi4(job_name, mol_name)
    #submits the slurm jon
    submit_slurm_job(job_name, mol_name)
    #goes back to previous dirctory
    os.chdir(os.path.dirname(os.getcwd()))
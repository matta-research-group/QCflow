from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
import os
import json


def run_opt_failed(job_name, mol_name, mol_smile, mol_dic, functional='B3LYP', basis_set='6-31G*'):
    """
    Submits a neutral optimization of failed jobs

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    mol_dic : Dictionary the oligomer came from

    functional : preset is B3LYP

    basis_set : preset is 6-31G*
    """
    #goes into directory
    os.chdir(f'{mol_name}')
    #turns smiles string into rdkit object
    mol = Chem.MolFromSmiles(mol_smile)
    #gets rdkit estimated coordinates of dimer
    mol3d = embed_molecule(mol)
    #Saves the torsional scan as .csv and finds the lowest energy geometry
    #conf_geo = torsion_parser(mol_name, mol_dic)
    conf_geo = mol3d.GetConformer()
    #writes a guassian input file
    write_gaussian(job_name, mol_name, mol_smile, functional, basis_set, mol3d, 0, conf_geo)
    #writes the slurm file
    write_slurm(job_name, mol_name)
    #submits the slurm jon
    submit_slurm_job(job_name, mol_name)
    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))

def run_opt_failed_from_chk(job_name, mol_name, mol_smile, mol_dic, functional='B3LYP', basis_set='6-31G*'):
    """
    Submits a neutral optimization of failed jobs

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    mol_dic : Dictionary the oligomer came from

    functional : preset is B3LYP

    basis_set : preset is 6-31G*
    """
    #goes into directory
    os.chdir(f'{mol_name}')
    #turns smiles string into rdkit object
    mol = Chem.MolFromSmiles(mol_smile)
    #gets rdkit estimated coordinates of dimer
    mol3d = embed_molecule(mol)
    #Saves the torsional scan as .csv and finds the lowest energy geometry
    #conf_geo = torsion_parser(mol_name, mol_dic)
    conf_geo = mol3d.GetConformer()
    #writes a guassian input file
    write_gaussian(job_name, mol_name, mol_smile, functional, basis_set, mol3d, 0, conf_geo, old_chk=f'{mol_name}_{job_name}.chk')
    #writes the slurm file
    write_slurm(job_name, mol_name)
    #submits the slurm jon
    submit_slurm_job(job_name, mol_name)
    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))

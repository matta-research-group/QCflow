from fragments import *
from find_torsion import *
from write_gaussian import *
from slurm import *
from torsion_parser import *
import os
import json

def run_neu_an_geo(mol_name, mol_dic, functional='xB97XD', basis_set='cc-pVDZ'):
    """
    Submits a neutral optimization with population analysis when provided with the name of the mol

    mol_name : name of oligomer

    mol_dic : Dictionary the oligomer came from

    functional : preset is wB97XD

    basis_set : preset is cc-pVDZ
    """
    #goes into directory
    os.chdir(f'{mol_name}')
    #writes a guassian input file
    jobs = ['n_a_geo','n_c_geo']

    for j in jobs:
        write_gaussian(j, mol_name, mol_dic[mol_name], functional, basis_set)
        #writes the slurm file
        write_slurm(j, mol_name)
        #submits the slurm jon
        submit_slurm_job(j, mol_name)

    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))

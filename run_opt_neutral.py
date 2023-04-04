from fragments import *
from find_torsion import *
from write_gaussian import *
from slurm import *
from torsion_parser import *
import os
import json


def run_opt_neutral(mol_dic, functional='wB97XD', basis_set='cc-pVDZ'):
    """
    Submits a neutral optimization with population analysis when provided with the name of the mol

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    mol_dic : Dictionary the oligomer came from

    functional : preset is wB97XD

    basis_set : preset is cc-pVDZ
    """
    for mol_name in mol_dic.keys():
        #goes into directory
        os.chdir(mol_name)
        #turns smiles string into rdkit object
        mol = Chem.MolFromSmiles(mol_dic[mol_name])
        #gets rdkit estimated coordinates of dimer
        mol3d = embed_molecule(mol)
        #Saves the torsional scan as .csv and finds the lowest energy geometry
        conf_geo = torsion_parser(mol_name, mol_dic)
        #writes a guassian input file
        write_gaussian('pop_opt_n', mol_name, mol_dic[mol_name], functional, basis_set, mol3d, 0, conf_geo)
        #writes the slurm file
        write_slurm('pop_opt_n', mol_name)
        #submits the slurm jon
        submit_slurm_job('pop_opt_n', mol_name)
        #goes back to previous directory
        os.chdir(os.path.dirname(os.getcwd()))

from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
import os
import json

def opt_calculation_from_mol(mol_name, mol_file, functional='B3LYP', basis_set='6-31G*'):

    """
    Submits a geometry optimization without the need for a torsional scan. 
    The dihedral angles will be zero with respect to how you draw them in ChemDraw.
    Make sure ChemDraws were done correctly and the angles will be correct

    mol_name : name of oligomer

    mol_file : filename of structure, file should be in the form of a .mol

    functional : preset is B3LYP

    basis_set : preset is 6-31G*
    """
    #loads the mol file information
    mol3d = Chem.rdmolfiles.MolFromMolFile(mol_file, sanitize=True, removeHs=False)

    #makes a directory
    os.mkdir(f'{mol_name}')
    #goes into directory
    os.chdir(f'{mol_name}')

    #turns mol (rdkit object) into smiles string
    mol_smi = Chem.MolToSmiles(mol3d)
    #canconicalises smiles string
    con_smi = Chem.CanonSmiles(mol_smi)
    #gets conformer
    conf = mol3d.GetConformer()
    #writes a gaussian input file
    write_gaussian('opt', mol_name, con_smi, functional, basis_set, mol=mol3d, torsion=0, conformer=conf)
    #writes the slurm job
    write_slurm('opt', mol_name)
    #submits the slurm job
    submit_slurm_job('opt', mol_name)
    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))



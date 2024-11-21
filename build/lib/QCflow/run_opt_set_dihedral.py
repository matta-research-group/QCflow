from QCflow.fragments import *
from QCflow.find_torsion import *
from QCflow.write_gaussian import *
from QCflow.slurm import *
from QCflow.torsion_parser import *
import os
import json

def run_opt_neutral_set_dihedral(mol_name, mol_smile, functional='B3LYP', basis_set='6-31G*', dihedral_angle=0):
    """
    Submits a neutral optimization with population analysis when provided with the name of the mol

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    functional : preset is B3LYP

    basis_set : preset is 6-31G*

    dihedral_angle : Angle to set dihedral to (preset is 0)
    """
    #makes a directory for the molecule
    os.mkdir(f'{mol_name}')
    #goes into directory
    os.chdir(f'{mol_name}')
    #turns smiles string into rdkit object
    mol = Chem.MolFromSmiles(mol_smile)
    #gets rdkit estimated coordinates of dimer
    mol3d = embed_molecule(mol)
    #Saves the torsional scan as .csv and finds the lowest energy geometry
    conf_geo = setting_dihedral(mol_smile, dihedral_angle)
    #writes a guassian input file
    write_gaussian('opt', mol_name, mol_smile, functional, basis_set, mol=mol3d, torsion=0, conformer=conf_geo)
    #writes the slurm file
    write_slurm('opt', mol_name)
    #submits the slurm jon
    submit_slurm_job('opt', mol_name)
    #goes back to previous directory
    os.chdir(os.path.dirname(os.getcwd()))

def find_dihedral_from_mol(mol3d):
    conf = mol3d.GetConformer() #makes a conformer
    rot_bond = getBond(mol3d) #finds the rotatable bond
    tor_bonds = getTorsion(mol3d, rot_bond[0]) #finds the 4 atoms in the dihedral
    #angle_of_dihedral = rdMolTransforms.GetDihedralDeg(conf, tor_bonds[0], tor_bonds[1], tor_bonds[2], tor_bonds[3])
    #returns the dihedral angle of those 4 atoms


    return tor_bonds[0], tor_bonds[1], tor_bonds[2], tor_bonds[3]

def rdkit_predict_dihedral(mol_smiles, num_of_conformer=100, max_iter=500, min_energy_MMFF=10000, min_energy_index_MMFF=0):

    # Number of conformers to be generated
    #num_of_conformer=100
    #max_iter=500
    # Default values for min energy conformer
    #min_energy_MMFF=10000
    #min_energy_index_MMFF=0

    mol = Chem.MolFromSmiles(mol_smiles)
    mol_h_MMFF = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol_h_MMFF, numConfs=num_of_conformer,params=AllChem.ETKDG()) # Generate conformers (stored in side the mol object)
    results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(mol_h_MMFF,maxIters=max_iter)

    for index, result in enumerate(results_MMFF):
        if(min_energy_MMFF>result[1]):
            min_energy_MMFF=result[1]
            min_energy_index_MMFF=index

    lowest_mol = Chem.Mol(mol_h_MMFF,False,min_energy_index_MMFF)

    conf = lowest_mol.GetConformer()

    atom1, atom2, atom3, atom4 = find_dihedral_from_mol(lowest_mol)

    predictaed_angle = Chem.rdMolTransforms.GetDihedralDeg(conf, atom1, atom2, atom3, atom4)

    return predictaed_angle, conf

def rdkit_predict_conf(mol_smiles, num_of_conformer=100, max_iter=500, min_energy_MMFF=10000, min_energy_index_MMFF=0):
    """
    Generates conformers for a given molecule using RDKit and returns the lowest energy conformer.

    Parameters:
    mol_smiles (str): The SMILES representation of the molecule.
    num_of_conformer (int): The number of conformers to be generated (default: 100).
    max_iter (int): The maximum number of iterations for conformer optimization (default: 500).
    min_energy_MMFF (float): The minimum energy threshold for selecting the lowest energy conformer (default: 10000).
    min_energy_index_MMFF (int): The index of the lowest energy conformer (default: 0).

    Returns:
    conf (Chem.Conformer): The lowest energy conformer of the molecule.

    """
    # Number of conformers to be generated
    #num_of_conformer=100
    #max_iter=500
    # Default values for min energy conformer
    #min_energy_MMFF=10000
    #min_energy_index_MMFF=0
    
    mol = Chem.MolFromSmiles(mol_smiles)
    mol_h_MMFF = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol_h_MMFF, numConfs=num_of_conformer,params=AllChem.ETKDG()) # Generate conformers (stored in side the mol object)
    results_MMFF = AllChem.MMFFOptimizeMoleculeConfs(mol_h_MMFF,maxIters=max_iter)

    for index, result in enumerate(results_MMFF):
        if(min_energy_MMFF>result[1]):       
            min_energy_MMFF=result[1]
            min_energy_index_MMFF=index

    lowest_mol = Chem.Mol(mol_h_MMFF,False,min_energy_index_MMFF)

    conf = lowest_mol.GetConformer()

    return conf


def run_opt_neutral_rdkit(mol_name, mol_smile, functional='B3LYP', basis_set='6-31G*'):
    """
    Submits a neutral optimization with population analysis when provided with the name of the mol

    mol_name : name of oligomer

    mol_smile : SMILE string of oligomer

    functional : preset is B3LYP

    basis_set : preset is 6-31G*

    dihedral_angle : Angle to set dihedral to (preset is 0)
    """
    #makes a directory for the molecule
    os.mkdir(f'{mol_name}')
    #goes into directory
    os.chdir(f'{mol_name}')
    #turns smiles string into rdkit object
    mol = Chem.MolFromSmiles(mol_smile)
    #gets rdkit estimated coordinates of dimer
    mol3d = embed_molecule(mol)
    #Saves the torsional scan as .csv and finds the lowest energy geometry
    conf_geo = rdkit_predict_conf(mol_smile)
    #writes a guassian input file
    write_gaussian('opt', mol_name, mol_smile, functional, basis_set, mol=mol3d, torsion=0, conformer=conf_geo)
    #writes the slurm file
    write_slurm('opt', mol_name)
    #submits the slurm jon
    submit_slurm_job('opt', mol_name)
    #goes back to previous dirctory
    os.chdir(os.path.dirname(os.getcwd()))

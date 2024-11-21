import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
import numpy as np
import QCflow
from QCflow.load_gaussian import *
from QCflow.energy_calculations import *
from QCflow.torsion_parser import *
from QCflow.find_torsion import *

def getBondLinkers(mol, linker_type):
    '''
    From a rdkit molecule finds the rotatable bond

    mol : rdkit readable string of oligomer
    '''
    if linker_type == 'single':
        #rotatable bonds
        pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]')

    if linker_type == 'double':
        pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]C=C-!@[R!$(*#*)&!D1]')

    if linker_type == 'imine':
        pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]N=C-!@[R!$(*#*)&!D1]')

    if linker_type == 'thio':
        pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]')


    bonds = mol.GetSubstructMatches(pattern)
    
    return bonds

def getTorsion_one(mol, bond):
    '''
    Gets the torsion for the torsional scan

    mol : rdkit readable string of oligomer

    bond : rRtatable bond
    '''
    # get neighbors of first atom in bond
    for atom in mol.GetAtomWithIdx(bond[0]).GetNeighbors():
        idx = atom.GetIdx()
        mass = atom.GetMass()
        if idx!=bond[1]: #excludes the other atom in the bond
            if mass > 13: # N, S, O get priority
                first=idx
                break
            if mass > 12: # otherwise C
                first=idx

    return first, bond[0], bond[1], bond[2]     

def getTorsion_two(mol, bond):
    '''
    Gets the torsion for the torsional scan

    mol : rdkit readable string of oligomer

    bond : rRtatable bond
    '''
    for atom in mol.GetAtomWithIdx(bond[3]).GetNeighbors():
        idx=atom.GetIdx()
        mass = atom.GetMass()
        if idx!=bond[0]: #excludes the other atom in the bond
            if mass > 13: # N, S, O get priority
                last=idx
                break
            if mass > 12: # otherwise C
                last=idx
    
    return bond[1], bond[2], bond[3], last

def find_planarity(angle):
    
    planar = np.absolute(np.cos(np.deg2rad(angle)))
    return planar

def finding_multi_planairty(mol_name, mol_smiles, linker_type):
    
    data = cclib.io.ccread(f'{mol_name}/{mol_name}_opt.log')

    mol = Chem.MolFromSmiles(mol_smiles)

    mol_h = AllChem.AddHs(mol)
    mol_3d = embed_molecule(mol_h)
    #Embeds the molecule
    conf = mol_3d.GetConformer()

    for i in range(conf.GetNumAtoms()):
        correct_pos = conf.SetAtomPosition(i, data.converged_geometries[0][i])
    #Uses this min energy torsion to give the correct geometry

    #bond = getBondLinkers(mol_3d, linker_type)

    if linker_type == 'thio':
        bond = getBondLinkers(mol, 'thio')
        
        thio_bond_1 = getTorsion(mol, bond[0])
        thio_bond_2 = getTorsion(mol, bond[1])
        
        angle_1 = rdMolTransforms.GetDihedralDeg(conf, thio_bond_1[0], thio_bond_1[1], thio_bond_1[2], thio_bond_1[3])
        angle_2 = rdMolTransforms.GetDihedralDeg(conf, thio_bond_2[0], thio_bond_2[1], thio_bond_2[2], thio_bond_2[3])

    if linker_type == 'triple':
        angle_1 = 0
        angle_2 = 0

    if (linker_type == 'double') or (linker_type == 'imine'):
        bond = getBondLinkers(mol_3d, linker_type)
        
        torsion_one = getTorsion_one(mol, bond[0])
        torsion_two = getTorsion_two(mol, bond[0])

        angle_1 = rdMolTransforms.GetDihedralDeg(conf, torsion_one[0], torsion_one[1], torsion_one[2], torsion_one[3])
        angle_2 = rdMolTransforms.GetDihedralDeg(conf, torsion_two[0], torsion_two[1], torsion_two[2], torsion_two[3])

    planar_1 = find_planarity(angle_1)
    planar_2 = find_planarity(angle_2)

    average_planarity = (planar_1 + planar_2) / 2

    return average_planarity


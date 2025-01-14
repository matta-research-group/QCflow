from QCflow.fragments import *
from QCflow.load_gaussian import *
from QCflow.find_torsion import *
import json
import csv
import matplotlib.pyplot as plt
import cclib
import rdkit
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

def find_min_energy_index(data):
    """
    Finds the index of the minimum energy from the torsional scan data.

    Parameters
    ----------
    data (cclib.parser.data.ccData_optdone_bool): Loaded cclib data from the .log file containing the torsional information.

    Returns
    -------
    int: The index corresponding to the minimum energy value in the torsional scan data.
    """
    #finds opt energy of each 10 deg scan
    data.scanenergies = data.scfenergies[data.optstatus == 4]
    #finds min energy of all opt 36 scans
    mini = np.where(data.scanenergies==np.min(data.scanenergies))
    #gets index corresponding to minimum energy
    mini_value = mini[0][0]

    return mini_value

def min_angle(data):
    """
    Finds the angle of the minimum energy torsion.

    Parameters
    ----------
    data (cclib.parser.data.ccData): Loaded cclib data from the .log file containing the torsional information.

    Returns
    -------
    float: The dihedral angle corresponding to the minimum energy torsion.
    """
    data.scanenergies = data.scfenergies[data.optstatus == 4]
    min_ang_loc = np.where(data.scanenergies==np.min(data.scanenergies))
    #Location of the min angle
    min_ang = data.scanparm[0][min_ang_loc[0][0]]
    #Dihedral angle of the min energy
    return min_ang

def torsion_parser(mol_name, mol_smi):
    """
    Parses the torsional profile of a molecule and returns the optimized geometry at the lowest energy minimum.

    Parameters
    ----------
    mol_name (str): The name of the molecule (dimer or trimer).
    mol_smi (str): SMILES string of the molecule.

    Returns
    -------
    rdkit.Chem.rdchem.Conformer: The optimized geometry at the lowest energy minimum.

    The function performs the following steps:
    1. Converts the SMILES string into an RDKit molecule object.
    2. Adds hydrogen atoms to the molecule.
    3. Embeds the molecule in 3D space.
    4. Loads the torsional data from a .log file.
    5. Finds the minimum energy torsion.
    6. Sets the atom positions to the geometry corresponding to the minimum energy torsion.
    """
    mol = Chem.MolFromSmiles(mol_smi)
    # Converts the SMILES string into an RDKit molecule object
    mol_h = AllChem.AddHs(mol)
    # Adds hydrogen atoms
    mol_3d = embed_molecule(mol_h)
    # Embeds the molecule in 3D space
    conf = mol_3d.GetConformer()
    # Gets the conformer of the molecule
    data = cclib.io.ccread(f'{mol_name}/{mol_name}_tor.log')
    # Loads the torsional data from a .log file
    min_energy = find_min_energy_index(data)
    # Finds the minimum energy torsion

    for i in range(conf.GetNumAtoms()):
        conf.SetAtomPosition(i, data.converged_geometries[min_energy][i])
    # Sets the atom positions to the geometry corresponding to the minimum energy torsion

    return conf

def setting_dihedral(mol_smile, deg):
    """
    Sets the dihedral angle of the torsion for individual scans.

    Parameters
    ----------
    mol_smile (str): SMILES string of the dimer.
    deg (float): Desired dihedral angle (0 or 180 for planar).

    Returns
    -------
    rdkit.Chem.rdchem.Conformer: A conformer of the dimer with the specified dihedral angle.

    The function performs the following steps:
    1. Converts the SMILES string into an RDKit molecule object.
    2. Embeds the molecule to get estimated 3D coordinates.
    3. Retrieves the conformer of the molecule.
    4. Identifies the bond and torsion angles between fragments.
    5. Sets the specified dihedral angle for the torsion.
    6. Applies MMFF force field and minimizes the energy of the molecule.
    """

    mol = Chem.MolFromSmiles(mol_smile)
    #turns smiles string into rdkit object
    mol3d = embed_molecule(mol)
    #gets rdkit estimated coordinates of dimer
    conf = mol3d.GetConformer()
    #getting the conformer

    bond = getBond(mol)
    #finds the bond between the fragment
    torsion = getTorsion(mol, bond[0])
    #torsion of the bond between the fragment
    t1, t2, t3, t4 = [t for t in torsion]
    #getting atoms involved in torsion
    rdkit.Chem.rdMolTransforms.SetDihedralDeg(conf, t1, t2, t3, t4, deg)
    #setting the dihedral angle
    mp = AllChem.MMFFGetMoleculeProperties(mol3d)
    #mol properties
    ff = AllChem.MMFFGetMoleculeForceField(mol3d, mp)
    #force field
    for i in torsion:
        ff.MMFFAddPositionConstraint(i, 0, 1.e4)
    ff.Minimize(maxIts=10000)

    return conf

def finding_dihedral_opt(mol_smiles, log_data):
    """
    Calculates the dihedral angle of a molecule given its SMILES string and log data.

    Parameters
    ----------
        mol_smiles (str): The SMILES string of the molecule.
        log_data (object): An object containing log data, including converged geometries.

    Returns
    -------
        float: The dihedral angle in degrees.
    """

    data = log_data

    mol = Chem.MolFromSmiles(mol_smiles)
    #Converts the SMILE into rdkit readable string
    mol_h = AllChem.AddHs(mol)
    #Adds hydrogen atoms
    mol_3d = embed_molecule(mol_h)
    #Embeds the molecule
    conf = mol_3d.GetConformer()

    for i in range(conf.GetNumAtoms()):
        correct_pos = conf.SetAtomPosition(i, data.converged_geometries[0][i])
    #Uses this min energy torsion to give the correct geometry

    bond = getBond(mol)

    torsion_atoms = getTorsion(mol, bond[0])

    angle = rdMolTransforms.GetDihedralDeg(conf, torsion_atoms[0], torsion_atoms[1], torsion_atoms[2], torsion_atoms[3])

    return angle

def find_planarity(angle):
    """
    Calculate the planarity of a given angle.
    This function computes the planarity of an angle by taking the absolute value of the cosine of the angle converted to radians.
    
    Parameters
    ----------
    angle (float): The angle in degrees for which the planarity is to be calculated.
    
    Returns
    -------
    float: The planarity value.
    """
    
    planar = np.absolute(np.cos(np.deg2rad(angle)))
    return planar

def getBondLinkers(mol, linker_type):
    """
    From an RDKit molecule, finds the two atoms involved in specified type of rotatable bond.
    
    Parameters
    ----------
    mol (rdkit.Chem.Mol): RDKit molecule object.
    linker_type (str): Type of linker to search for. Options are 'single', 'double', 'imine', or 'thio'.
    
    Returns
    -------
    tuple: A tuple, where indices of the atoms involved in the matching bonds.
    """

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
    """
    Gets the first torsion of a multi torsion molecule. Works for triple, imine and double bonds.

    Parameters
    ----------
    mol (rdkit.Chem.Mol): RDKit molecule object representing the oligomer.
    bond (tuple): Tuple of atom indices representing the rotatable bond.

    Returns
    -------
    tuple: A tuple containing the indices of the four atoms defining the torsion angle.
    """

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
    """
    Gets the second torsion of a multi torsion molecule. Works for triple, imine and double bonds.

    Parameters
    ----------
    mol (rdkit.Chem.Mol): RDKit molecule object representing the oligomer.
    bond (tuple): Tuple of atom indices representing the rotatable bond.

    Returns
    -------
    tuple: A tuple containing the indices of the four atoms defining the torsion angle.
    """
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

def finding_multi_planairty(mol_name, mol_smiles, linker_type):
    """
    Determines the average planarity of a molecule based on its torsion angles.
    
    Parameters
    ----------
    mol_name (str): The name of the molecule, used to locate the optimization log file.
    mol_smiles (str): The SMILES representation of the molecule.
    linker_type (str): The type of linker in the molecule. Can be 'thio', 'triple', 'double', or 'imine'.
    
    Returns
    -------
    float: The average planarity of the molecule.
    """
    
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
        
        torsion_one = getTorsion_one(mol, bond[0]) #removed bond[0][0]
        torsion_two = getTorsion_two(mol, bond[0]) #removed bond[0][0]

        angle_1 = rdMolTransforms.GetDihedralDeg(conf, torsion_one[0], torsion_one[1], torsion_one[2], torsion_one[3])
        angle_2 = rdMolTransforms.GetDihedralDeg(conf, torsion_two[0], torsion_two[1], torsion_two[2], torsion_two[3])

    planar_1 = find_planarity(angle_1)
    planar_2 = find_planarity(angle_2)

    average_planarity = (planar_1 + planar_2) / 2

    return average_planarity

def update_conformer_from_xyz(mol, xyz_file):
    """
    Updates the conformer of an RDKit molecule using coordinates from an XYZ file.
    
    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule.
        xyz_file (str): Path to the XYZ file containing new coordinates.
    
    Returns:
        rdkit.Chem.Mol: The molecule with updated conformer.
    """
    # Read the XYZ file
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    
    # Parse number of atoms from the first line
    num_atoms = int(lines[0].strip())
    
    # Parse the coordinates
    coordinates = []
    atom_symbols = []
    for line in lines[2:2 + num_atoms]:  # Skip the first two lines
        parts = line.split()
        atom_symbols.append(parts[0])  # Atom symbol
        coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])  # x, y, z
    
    # Ensure atom count matches
    if len(coordinates) != mol.GetNumAtoms():
        raise ValueError("Number of atoms in the XYZ file does not match the RDKit molecule.")
    
    # Create or update a conformer
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, coord in enumerate(coordinates):
        conf.SetAtomPosition(i, coord)
    
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    
    return mol


def update_conformer_from_xyz(mol, xyz_file):
    """
    Updates the conformer of an RDKit molecule using coordinates from an XYZ file.
    
    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule.
        xyz_file (str): Path to the XYZ file containing new coordinates.
    
    Returns:
        rdkit.Chem.Mol: The molecule with updated conformer.
    """
    # Read the XYZ file
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    
    # Parse number of atoms from the first line
    num_atoms = int(lines[0].strip())
    
    # Parse the coordinates
    coordinates = []
    atom_symbols = []
    for line in lines[2:2 + num_atoms]:  # Skip the first two lines
        parts = line.split()
        atom_symbols.append(parts[0])  # Atom symbol
        coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])  # x, y, z
    
    # Ensure atom count matches
    if len(coordinates) != mol.GetNumAtoms():
        raise ValueError("Number of atoms in the XYZ file does not match the RDKit molecule.")
    
    # Create or update a conformer
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, coord in enumerate(coordinates):
        conf.SetAtomPosition(i, coord)
    
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    
    return mol

def finding_planairty_psi4(mol_name, mol_smiles, linker_type, job_name):
    """
    Determines the average planarity of a molecule based on its torsion angles.
    
    Parameters
    ----------
    mol_name (str): The name of the molecule, used to locate the optimization xyz file.
    mol_smiles (str): The SMILES representation of the molecule.
    linker_type (str): The type of linker in the molecule. Can be 'single', 'thio', 'triple', 'double', or 'imine'.
    job_name (str): The type of job to run. Possible values:
        - 'sp': Single Point neutral
        - 'opt': Optimisation neutral
    
    Returns
    -------
    float: The average planarity of the molecule.
    """
    
    #data = cclib.io.ccread(f'{mol_name}/{mol_name}_opt.xyz')
    data = f'{mol_name}/{mol_name}_{job_name}.xyz'
    

    mol = Chem.MolFromSmiles(mol_smiles)

    mol_h = AllChem.AddHs(mol)
    #Embeds the molecule
    mol_update = update_conformer_from_xyz(mol_h, data)
    conf = mol_update.GetConformer()

    if linker_type == 'single':

        bond = getBond(mol)
        torsion_atoms = getTorsion(mol, bond[0])

        angle = rdMolTransforms.GetDihedralDeg(conf, torsion_atoms[0], torsion_atoms[1], torsion_atoms[2], torsion_atoms[3])
        single_planarity = find_planarity(angle)


    if linker_type == 'thio':
        bond = getBondLinkers(mol_update, 'thio')
        
        thio_bond_1 = getTorsion(mol, bond[0])
        thio_bond_2 = getTorsion(mol, bond[1])
        
        angle_1 = rdMolTransforms.GetDihedralDeg(conf, thio_bond_1[0], thio_bond_1[1], thio_bond_1[2], thio_bond_1[3])
        angle_2 = rdMolTransforms.GetDihedralDeg(conf, thio_bond_2[0], thio_bond_2[1], thio_bond_2[2], thio_bond_2[3])

    if linker_type == 'triple':
        angle_1 = 0
        angle_2 = 0

    if (linker_type == 'double') or (linker_type == 'imine'):
        bond = getBondLinkers(mol_update, linker_type)
        
        torsion_one = getTorsion_one(mol_update, bond[0])
        torsion_two = getTorsion_two(mol_update, bond[0])

        angle_1 = rdMolTransforms.GetDihedralDeg(conf, torsion_one[0], torsion_one[1], torsion_one[2], torsion_one[3])
        angle_2 = rdMolTransforms.GetDihedralDeg(conf, torsion_two[0], torsion_two[1], torsion_two[2], torsion_two[3])

    if (linker_type == 'double') or (linker_type == 'imine') or (linker_type == 'thio') or (linker_type == 'triple'):

        planar_1 = find_planarity(angle_1)
        planar_2 = find_planarity(angle_2)

        planarity = (planar_1 + planar_2) / 2

    if (linker_type == 'single'):
        planarity = single_planarity

    return planarity
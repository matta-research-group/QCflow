import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def getBond(mol):
    """
    Finds the rotatable bonds in a given RDKit molecule.

    Parameters
    ----------
    mol (rdkit.Chem.Mol): An RDKit molecule object.

    Returns
    -------
    tuple: A tuple of tuples, where each inner tuple contains the indices of atoms that form a rotatable bond.
    """
    #rotatable bonds
    pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]')
    bonds = mol.GetSubstructMatches(pattern)
    return bonds

def getTorsion(mol,bond):
    """
    Gets the torsion for the torsional scan.

    Parameters
    ----------
    mol (rdkit.Chem.Mol): RDKit molecule object representing the oligomer.
    bond (tuple): Tuple of two integers representing the indices of the rotatable bond.

    Returns
    -------
    tuple: A tuple of four integers representing the indices of the atoms involved in the torsion.
            The format is (first_atom, bond_atom1, bond_atom2, last_atom).
            'first_atom' is the neighbor of 'bond_atom1' with the highest priority (N, S, O > C).
            'last_atom' is the neighbor of 'bond_atom2' with the highest priority (N, S, O > C).
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
    # get neighbors of second atom in bond
    for atom in mol.GetAtomWithIdx(bond[1]).GetNeighbors():
        idx=atom.GetIdx()
        mass = atom.GetMass()
        if idx!=bond[0]: #excludes the other atom in the bond
            if mass > 13: # N, S, O get priority
                last=idx
                break
            if mass > 12: # otherwise C
                last=idx

    return (first, bond[0], bond[1], last)

def embed_molecule(mol):
    """
    Generates the 3D rdkit.Chem.Mol object of the given RDKit molecule and adds hydrogen atoms.

    Parameters
    ----------
    mol (rdkit.Chem.Mol): An RDKit molecule object.

    Returns
    -------
    rdkit.Chem.Mol: The RDKit molecule object with embedded 3D coordinates and added hydrogen atoms.
    """

    addhs = Chem.AddHs(mol)
    AllChem.EmbedMolecule(addhs) # the embedded molecule
    return addhs

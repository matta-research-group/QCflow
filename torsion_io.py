import numpy as np
import matplotlib.pyplot as plt
import rdkit
import cclib
import itertools
import os
import os.path
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests
from itertools import combinations
from typing import List

def getBond(mol):
    #rotatable bonds
    pattern = Chem.MolFromSmarts('[R!$(*#*)&!D1]-!@[R!$(*#*)&!D1]')
    bonds = mol.GetSubstructMatches(pattern)
    return bonds

#3 14 21 22 S 35 10.0
def getTorsion(mol,bond):
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
    addhs = Chem.AddHs(mol)
    AllChem.EmbedMolecule(addhs) # the embedded molecule
    return addhs

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


def read_fragments(name_of_smi_file):
    """
    reads fragment SMILES from file
    returns a dictionary of fragments
    """
    smile_list = np.genfromtxt(name_of_smi_file, dtype='str')
    #creates a dictionary of fragments
    fragment_dic = { i : d for i, d in enumerate(smile_list)}
    return fragment_dic

def make_dimers_dic(fragment_dic):
    """
    Reads fragment dictionary
    generates SMILES for all possible dimers
    returns dictionary of dimers as SMILES
    """
    dimer_smiles = (f'{v1}.{v2}'.replace('(*)', '9')
                    for v1, v2 in itertools.combinations(fragment_dic.values(), r=2))
    #All the combinations of dimer names from fragment dictionary
    dimer_names = [(f'{k1}_{k2}')
                    for k1, k2 in itertools.combinations(fragment_dic.keys(), r=2)]

    # Convert all SMILES to mol objects
    #mols = [Chem.MolFromSmiles(smi) for smi in dimer_smiles]

    #relating the dimer names and dimer smiles back to one another
    dimer_dic = { k : v for k, v in zip(dimer_names, dimer_smiles) }

    return dimer_dic

def smiles_to_mols(SMILES):
    '''
    Provide a dimer SMILE string which is turned into a rdkit molecule
    '''
    mol = Chem.MolFromSmiles(SMILES)
    return mol

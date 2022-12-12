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


def fragment_dictionary(name_of_smi_file):
    #loads the smile file
    smile_list = np.genfromtxt(f'{name_of_smi_file}.smi',dtype='str')
    #creates a dictionary of fragments
    fragment_dic = { i : d for i, d in enumerate(smile_list)}
    return fragment_dic

def combined_fragments_dic(dic):
    # Get all SMILES for all possible combinations (of 2 fragments) of the fragment list
    #All the combinations of dimers smiles from fragment dictionary
    combo_smiles = (f'{v1}.{v2}'.replace('(*)', '9')
                    for v1, v2 in itertools.combinations(dic.values(), r=2))
    #All the combinations of dimer names from fragment dictionary
    combo_names = [(f'{k1}-{k2}')
                    for k1, k2 in itertools.combinations(dic.keys(), r=2)]

        # Convert all SMILES to mol objects
    #mols = [Chem.MolFromSmiles(smi) for smi in combo_smiles]

    #relating the dimer names and dimer smiles back to one another
    combo_dic = { k : v for k, v in zip(combo_names, combo_smiles) }

    return combo_dic
    #Creates a list of all the possible combinations as RDKit info

def combined_fragments_dic_mol(dic):
    # Get all SMILES for all possible combinations (of 2 fragments) of the fragment list
    #All the combinations of dimers smiles from fragment dictionary
    combo_smiles = (f'{v1}.{v2}'.replace('(*)', '9')
                    for v1, v2 in itertools.combinations(dic.values(), r=2))
    #All the combinations of dimer names from fragment dictionary
    combo_names = [(f'{k1}-{k2}')
                    for k1, k2 in itertools.combinations(dic.keys(), r=2)]

        # Convert all SMILES to mol objects
    mols = [Chem.MolFromSmiles(smi) for smi in combo_smiles]

    #relating the dimer names and dimer smiles back to one another
    combo_dic = { k : v for k, v in zip(combo_names, mols) }

    return combo_dic
    #Creates a list of all the possible combinations as RDKit info

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

def make_dimer_smiles(a, b):
    return a.format('', '8') + '.' + b.format('8', '')

def make_trimer_smiles(a, b):
    aba = a.format('', '8') + '.' + b.format('8', '9') + '.' + a.format('9', '')
    bab = b.format('', '8') + '.' + a.format('8', '9') + '.' + b.format('9', '')
    return aba, bab

def read_fragments(name_of_smi_file):
    """
    reads fragment SMILES from file
    returns a dictionary of fragments
    """
    smile_list = np.genfromtxt(name_of_smi_file, dtype='str')

    smi_mol = [str(smi).replace("(*)", "{}") for smi in smile_list]

    #creates a dictionary of fragments
    fragment_dic = { i : d for i, d in enumerate(smi_mol)}
    return fragment_dic

def make_dimer_dic(fragment_dic):

    """
    Reads fragment dictionary
    generates SMILES for all possible dimers
    returns dictionary of dimers as SMILES
    """

    mol_smiles = [make_dimer_smiles(v1, v2)
                    for v1, v2 in itertools.combinations(fragment_dic.values(), r=2)]

    mol_names = [(f'{k1}_{k2}')
                    for k1, k2 in itertools.combinations(fragment_dic.keys(), r=2)]

    mol_dic = { k : v for k, v in zip(mol_names, mol_smiles) }

    return mol_dic

def make_trimer_dic(fragment_dic):

    """
    Reads fragment dictionary
    generates SMILES for all possible trimers
    returns dictionary of trimers as SMILES
    """

    mol_smiles = [make_trimer_smiles(v1, v2)
                    for v1, v2 in itertools.combinations(fragment_dic.values(), r=2)]
    #Makes all the possible trimers aba and bab
    k1k2k1 = []
    k2k1k2 = []
    for seperate in mol_smiles:
        aba = seperate[0]
        bab = seperate[1]

        k1k2k1.append(aba)
        k2k1k2.append(bab)
    #Makes the two lists of aba and bab

    aba = [(f'{k1}_{k2}_{k1}')
                    for k1, k2 in itertools.combinations(fragment_dic.keys(), r=2)]
    #combinations of keys for aba
    bab = [(f'{k2}_{k1}_{k2}')
                    for k1, k2 in itertools.combinations(fragment_dic.keys(), r=2)]
    #combinations of keys for bab
    aba_dic = { k : v for k, v in zip(aba, k1k2k1) }
    #Makes dictionary of all the aba keys with the aba trimers
    bab_dic = { k : v for k, v in zip(bab, k2k1k2) }
    #Makes dicitionary of all the bab keys with the bab trimers
    trimer_dic = aba_dic | bab_dic
    #Combines these two dictionaries together
    return trimer_dic

import numpy as np
import rdkit
import itertools
import os
import os.path
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import requests
from itertools import combinations
from typing import List

def make_dimer_smiles(a, b):
    """
    Combines two fragments together to make a dimer

    a : SMILE string of a fragment with attachment points

    b : SMILE string of a fragment with attachment points
    """
    #todo: should be also returning the dimer with
    # a.format('','8') + b.format('', '8')
    # and potentially all combinations (?)
    # to take into account asymmetric/functionalised rings
    return a.format('', '8') + '.' + b.format('8', '')


def make_trimer_smiles(a, b):
    """
    Combines two fragments to make a trimer with both the aba and bab
    configuration produced

    a : SMILE string of a fragment with attachment points

    b : SMILE string of a fragment with attachment points
    """
    aba = a.format('', '8') + '.' + b.format('8', '9') + '.' + a.format('9', '')
    bab = b.format('', '8') + '.' + a.format('8', '9') + '.' + b.format('9', '')
    return aba, bab

def read_fragments(name_of_smi_file):
    """
    Reads fragment SMILES from a .smi file and returns a
    dictionary of fragments

    name_of_smi_file : name of the .smi file where the fragments are
    """
    smile_list = np.genfromtxt(name_of_smi_file, dtype='str')
    #The format of the attachment points are changed as a space is required for combination
    smi_mol = [str(smi).replace("(*)", "{}") for smi in smile_list]

    #creates a dictionary of fragments
    fragment_dic = { i : d for i, d in enumerate(smi_mol)}
    return fragment_dic

def make_dimer_dic(fragment_dic):

    """
    When provided the fragment dictionary, generates SMILES for all possible
    dimers and returns dictionary of dimers.

    Where the key is the name of the dimer, if fragment 0 is combined with
    fragment 1 then the name is 0_1.

    The value is the SMILE string of the dimer

    fragment_dic : The dictionary of all the fragments

    """

    mol_smiles = [make_dimer_smiles(v1, v2)
                    for v1, v2 in itertools.combinations(fragment_dic.values(), r=2)]

    mol_names = [(f'{k1}_{k2}')
                    for k1, k2 in itertools.combinations(fragment_dic.keys(), r=2)]

    mol_dic = { k : v for k, v in zip(mol_names, mol_smiles) }

    return mol_dic

def make_dimer_dic_from_2_dic(fragment_dic_1,fragment_dic_2):

    """
    When provided 2 different fragment dictionaries, generates SMILES for all possible
    dimers and returns dictionary of dimers.

    Where the key is the name of the dimer, if fragment 0 is combined with
    fragment 1 then the name is 0_1.

    The value is the SMILE string of the dimer

    mol_dic : The dictionary of all the fragments

    """

    mol_smiles = [make_dimer_smiles(v1, v2)
                    for v1, v2 in itertools.product(list(fragment_dic_1.values()),list(fragment_dic_2.values()))]

    mol_names = [(f'{k1}_{k2}')
                    for k1, k2 in itertools.product(fragment_dic_1.keys(),fragment_dic_2.keys())]

    mol_dic = { k : v for k, v in zip(mol_names, mol_smiles) }

    return mol_dic

def make_trimer_dic(fragment_dic):

    """
    When provided the fragment dictionary, generates SMILES for all possible
    trimers (ABA & BAB) and returns dictionary of trimers.

    Where the key is the name of the trimer, if fragment 0 is combined with
    fragment 1 then the name is 0_1_0 & 1_0_1.

    The value is the SMILE string of the trimer

    fragment_dic : The dictionary of all the fragments
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

def make_trimer_dic_from_2_dic(frag_dic_1,frag_dic_2):

    """
    When provided 2 fragment dictionaries, generates SMILES for all possible
    trimers (ABA & BAB) and returns dictionary of trimers.

    Where the key is the name of the trimer, if fragment 0 is combined with
    fragment 1 then the name is 0_1_0 & 1_0_1.

    The value is the SMILE string of the trimer

    fragment_dic : The dictionary of all the fragments
    """

    mol_smiles = [make_trimer_smiles(v1, v2)
                    for v1, v2 in itertools.product(list(frag_dic_1.values()),list(frag_dic_2.values()))]
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
                    for k1, k2 in itertools.product(frag_dic_1.keys(),frag_dic_2.keys())]
    #combinations of keys for aba
    bab = [(f'{k2}_{k1}_{k2}')
                    for k1, k2 in itertools.product(frag_dic_1.keys(),frag_dic_2.keys())]
    #combinations of keys for bab
    aba_dic = { k : v for k, v in zip(aba, k1k2k1) }
    #Makes dictionary of all the aba keys with the aba trimers
    bab_dic = { k : v for k, v in zip(bab, k2k1k2) }
    #Makes dicitionary of all the bab keys with the bab trimers
    trimer_dic = aba_dic | bab_dic
    #Combines these two dictionaries together
    return trimer_dic

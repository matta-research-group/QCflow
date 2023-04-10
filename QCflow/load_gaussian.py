from QCflow.fragments import *
from QCflow.torsion_run import *
from QCflow.write_gaussian import *
import json
import csv
import cclib

def load_data(mol_name, job_name):
    """
    Opens and reads the .log data file of a torsional scan

    mol_name : number name of dimer or trimer

    job_name : The type of job run. Possible runs:
                Torsional scan neutral → tor
                Optimisation neutral/Population analysis → pop_opt_n
                Vertical anion → ver_a
                Vertical cation → ver_c
                Optimisation anion → opt_a
                Optimisation cation → opt_c
    """

    data = cclib.io.ccread(f'{mol_name}_{job_name}.log')

    return data

def save_dictionary(mol_dic, name_of_dic):
    """
    Saves the created dictionary as a .json file

    mol_dic : dictionary of oligomers. Where the key is the number of the
    oligomer and the value is the SMILE string

    name_of_dic : desired name of the json file i.e. name_of_dic.json
    """
    with open((name_of_dic), "w") as fp:
        json.dump(mol_dic,fp)

def open_dictionary(dictionary_file):
    """
    Opens a saved dictionary file

    dictionary_file : name of the .json file in which the dictionary is saved in
    """
    with open((dictionary_file), 'r') as f:
      mol_dic = json.load(f)

    return mol_dic

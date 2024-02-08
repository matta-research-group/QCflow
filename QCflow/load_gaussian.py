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
                Vertical anion → sp_a
                Vertical cation → sp_c
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

def data_dic(mol_dic, job_type):
    """
    Takes a dictionary of oligomers and the calculation performed.
    Returns a dictionary of parsed rdkit objects with the key being the name of the oligomer and the
    value being the cclib object.

    mol_dic : dictionary of oligomers and their SMILE strings

    job_type : Type of calculation performed

    """
    parsed_dic = {} #empty dictionary
    for k in mol_dic.keys():
        data = load_data(k, job_type) #parses log file
        parsed_dic[f'{k}'] = data
    return parsed_dic

def clean_directory(mol_dic, file_type, dir_loc):
    """
    Removes unwanted files

    mol_dic : Dictionary of oligomers in which all the calculations in which
                the user wants to clean up the directories

    file_type : File extension, as string ('.txt', '.err', '.out')

    dir_loc : Location of the main directory in which all the oligomers directories
                are located. Example('/scratch/prj/mime/tor_input')
    """
    
    for k, v in mol_dic.items():
        dir_name = f'{dir_loc}/{k}'
        #The directory of the calculations
        location = os.listdir(dir_name)
        for item in location:
            if item.endswith(str(file_type)):
                os.remove(os.path.join(dir_name, item))

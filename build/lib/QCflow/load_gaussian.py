import json
import csv
import cclib

def load_data(mol_name, job_name):
    """
    Opens and reads the .log data file.

    Parameters:
        mol_name (str): Name of molecule.
        job_name (str): The type of job run. Possible runs:
            - Single Point neutral → sp
            - Optimisation neutral → opt
            - Torsional scan neutral → tor
            - Optimisation neutral + Population analysis → pop_opt_n
            - Single point anion → sp_a
            - Single point cation → sp_c
            - Optimisation anion → opt_a
            - Optimisation cation → opt_c
            - neutral charge, optimised anion geometry → n_a_geo
            - neutral charge, optimised cation geometry → n_c_geo
            - Single Point Hirshfeld → sp_hirsh

    Returns:
        cclib.parser.data.ccData: Parsed data from the .log file.
    """

    data = cclib.io.ccread(f'{mol_name}/{mol_name}_{job_name}.log')

    return data

def save_dictionary(mol_dic, name_of_dic):
    """
    Saves the created dictionary as a .json file.

    Parameters:
    mol_dic (dict): Dictionary of oligomers where the key is the number of the oligomer and the value is the SMILE string.
    name_of_dic (str): Desired name of the json file (without the .json extension).

    Returns:
    None
    """
    with open((name_of_dic), "w") as fp:
        json.dump(mol_dic,fp)

def open_dictionary(dictionary_file):
    """
    Opens a saved dictionary file.

    Parameters:
    dictionary_file (str): The name of the .json file in which the dictionary is saved.

    Returns:
    dict: The dictionary loaded from the specified .json file.
    """

    with open((dictionary_file), 'r') as f:
      mol_dic = json.load(f)

    return mol_dic

def data_dic(mol_dic, job_type):
    """
    Takes a dictionary of oligomers and the calculation performed, and returns a dictionary of parsed cclib objects.

    Parameters:
    mol_dic (dict): A dictionary where keys are the names of oligomers and values are their SMILE strings.
    job_type (str): The type of calculation performed.

    Returns:
    dict: A dictionary where keys are the names of the oligomers and values are the parsed cclib objects.
    """

    parsed_dic = {} #empty dictionary
    for k in mol_dic.keys():
        data = load_data(k, job_type) #parses log file
        parsed_dic[f'{k}'] = data
    return parsed_dic

import json
import cclib
from QCflow.load_gaussian import *

def success_test(parsed_dic):
    """
    Tests if the calculation has been successful for each cclib parsed object.

    Parameters
    ----------
    parsed_dic (dict): A dictionary where the key is the name of the oligomer and the value is the cclib object.

    Returns
    -------
    tuple: Two lists - the first list contains the names of oligomers that passed, and the second list contains the names of oligomers that failed.
    """
    passed = []
    failed = []
    for k, v in parsed_dic.items():
        if v.metadata['success'] == True:
            passed.append(k)
        if v.metadata['success'] == False:
            failed.append(k)
    return passed, failed

def functional_test(parsed_dic, functional):
    """
    Tests if the calculations in the parsed cclib objects were done using the specified functional.

    Parameters
    ----------
    parsed_dic (dict): A dictionary where the key is the name of the oligomer and the value is the cclib object.
    functional (str): The desired functional of the calculation.

    Returns
    -------
    tuple: Two lists - the first list contains the names of oligomers that passed the test, 
            and the second list contains the names of oligomers that failed the test.
    """

    passed = []
    failed = []
    for k, v in parsed_dic.items():
        if v.metadata['functional'] == f'{functional}':
            passed.append(k)
        if v.metadata['functional'] != f'{functional}':
            failed.append(k)
    return passed, failed

def basis_set_test(parsed_dic, basis_set):
    """
    Tests if the calculations in the parsed dictionary have been done using the correct basis set.

    Parameters
    ----------
    parsed_dic (dict): A dictionary of parsed rdkit objects with the key being the name of the oligomer 
                       and the value being the cclib object.
    basis_set (str): The desired basis set for the calculation.

    Returns
    -------
    tuple: Two lists - the first list contains the names of oligomers that passed the basis set check,
           and the second list contains the names of oligomers that failed the basis set check.
    """

    passed = []
    failed = []
    for k, v in parsed_dic.items():
        if v.metadata['basis_set'] == f'{basis_set}':
            passed.append(k)
        if v.metadata['basis_set'] != f'{basis_set}':
            failed.append(k)
    return passed, failed

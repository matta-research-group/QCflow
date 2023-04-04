import json
import cclib

def test_success(parsed_dic):
    """
    When given rdkit parsed object it tests if the calculation has been successful.
    Returns two lists (passed and failed)

    parsed_dic : dictionary of parsed rdkit ojects with key being name of oligomer and
                value being the cclib object
    """
    passed = []
    failed = []
    for k, v in parsed_dic.items():
        if v.metadata['success'] == True:
            passed.append(k)
        if v.metadata['success'] == False:
            failed.append(k)
    return passed, failed

def test_functional(parsed_dic, functional):
    """
    When given rdkit parsed object it tests if the calculation has been done at correct functional.
    Returns two lists (passed and failed)

    parsed_dic : dictionary of parsed rdkit ojects with key being name of oligomer and
                value being the cclib object

    functional : desired functional of calculation
    """
    passed = []
    failed = []
    for k, v in parsed_dic.items():
        if v.metadata['functional'] == f'{functional}':
            passed.append(k)
        if v.metadata['functional'] != f'{functional}':
            failed.append(k)
    return passed, failed

def test_basis_set(parsed_dic, basis_set):
    """
    When given rdkit parsed object it tests if the calculation has been done at correct basis set.
    Returns two lists (passed and failed)

    parsed_dic : dictionary of parsed rdkit ojects with key being name of oligomer and
                value being the cclib object

    basis_set : desired basis_set of calculation
    """
    passed = []
    failed = []
    for k, v in parsed_dic.items():
        if v.metadata['functional'] == f'{basis_set}':
            passed.append(k)
        if v.metadata['functional'] != f'{basis_set}':
            failed.append(k)
    return passed, failed

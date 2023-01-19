from fragments_io import *
from torsion_io import *
from write_input_io import *
from file_io import *
import json
import csv

def load_torsional_data(mol_name, job_type):
    """
    Opens and reads the .log data file of a torsional scan
    """

    data = cclib.io.ccread(f'{mol_name}_{job_type}.log')

    return data

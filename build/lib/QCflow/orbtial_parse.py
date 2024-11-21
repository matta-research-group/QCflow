import QCflow
from QCflow.load_gaussian import *

def open_gaussian_output(job_name, mol_name):
    """
    Opens the Gaussian output file for a given job and molecule name.

    Parameters:
    job_name (str): The name of the job.
    mol_name (str): The name of the molecule.

    Returns:
    list: A list containing the lines of the log file.
    """
    gauss_out = f'{mol_name}/{mol_name}_{job_name}.log'
    log_file_lines = [line for line in open(gauss_out, 'r')]

    return log_file_lines


def homo_lumo_parse(opened_file):
    """
    Parses the given log file to extract the highest occupied molecular orbital (HOMO),
    lowest unoccupied molecular orbital (LUMO), and the energy gap between them.

    Args:
        opened_file (list): A list of strings representing the lines of the log file.

    Returns:
        tuple: A tuple containing the HOMO energy (in eV), LUMO energy (in eV), and the energy gap (in eV).

    This function was inspired by https://github.com/t-young31/g09_scripts/blob/master/g09extract.py
    """

    last_occ_eigenvalue = None
    first_virt_eigenvalue = None
    hartree_to_eV = 27.21138505

    homo = []
    lumo = []

    log_file_lines = opened_file

    for line in reversed(log_file_lines):

            if 'Alpha  occ. eigenvalues' in line and last_occ_eigenvalue is None:
                last_occ_eigenvalue = line.split()[-1]
                e_homo = last_occ_eigenvalue
                homo_float = (float(e_homo))
                homo.append(homo_float)

            if 'Alpha virt. eigenvalues' in line and last_occ_eigenvalue is None:
                first_virt_eigenvalue = line.split()[4]
                e_lumo = first_virt_eigenvalue
                lumo_float = (float(e_lumo))
                lumo.append(lumo_float)

    homo_eV = homo[-1] * hartree_to_eV
    lumo_eV = lumo[-1] * hartree_to_eV
    energy_gap = lumo_eV - homo_eV

    return homo_eV, lumo_eV, energy_gap


def orbtials_calc(job_name, mol_name):
    """
    Calculates the orbital properties for a given job and molecule.

    Parameters:
    job_name (str): The name of the job.
    mol_name (str): The name of the molecule.

    Returns:
    tuple: A tuple containing the passed molecules, homo energies, lumo energies, and energy gaps.
    """
    log_file_lines = open_gaussian_output(job_name, mol_name)
   
    passed_molecules = []
    homo = []
    lumo = []
    energy_gap = []

    for line in log_file_lines:
        if line.startswith(" Normal termination of Gaussian"):
                passed_molecules.append(mol_name)

                homo_eV, lumo_eV, energy_gap_eV = homo_lumo_parse(log_file_lines)
                
                homo.append(homo_eV)
                lumo.append(lumo_eV)
                energy_gap.append(energy_gap_eV)

    
    return passed_molecules, homo, lumo, energy_gap

def parse_orbtials(job_name, mol_dic):

    homo = {}
    lumo = {}
    energy_gap = {}
    successful_mol = []

    for k, v in mol_dic.items():
        passed_molecules, homo_eV, lumo_eV, energy_gap_eV = orbtials_calc(job_name, k)

        homo[k] = homo_eV[0]
        lumo[k] = lumo_eV[0]
        energy_gap[k] = energy_gap_eV[0]
        successful_mol.append(passed_molecules)

    return homo, lumo, energy_gap, successful_mol
from QCflow.load_gaussian import *
import cclib
import numpy as np

def extract_data_from_txt(file_path):
    """
    Extracts data from a .txt file with energy, HOMO, LUMO, and energy gap information.

    Parameters
    ----------
    file_path (str): Path to the .txt file.

    Returns
    -------
    dict: A dictionary with extracted values.
    """
    data = {}

    with open(file_path, 'r') as file:
        for line in file:
            if "Optimized energy:" in line:
                data["optimized_energy"] = float(line.split(":")[1].split()[0])  # Extract energy in eV
            elif "Singl Point energy:" in line:
                data["sp_energy"] = float(line.split(":")[1].split()[0]) #Extract single point energy in eV
            elif "HOMO:" in line:
                data["homo"] = float(line.split(":")[1].split()[0])  # Extract HOMO in eV
            elif "LUMO:" in line:
                data["lumo"] = float(line.split(":")[1].split()[0])  # Extract LUMO in eV
            elif "Energy gap" in line:
                data["energy_gap"] = float(line.split(":")[1].split()[0])  # Extract energy gap in eV

    return data

def cal_reorg(opt_n,sp_c,opt_c,n_c_geo, calculation_software='Gaussian'):
    """
    Calculate the reorganization energy.

    Parameters
    ----------

    calculation_software (str): The computational chemistry software used to perform the calculations. Default is 'Gaussian'. Can be 'Gaussian' or 'Psi4'.

    For Gaussian:

    opt_n (cclib.io.ccread): cclib object for the neutral population optimization analysis.
    sp_c (cclib.io.ccread): cclib object for the vertical anion or cation.
    opt_c (cclib.io.ccread): cclib object for the optimized anion or cation.
    n_c_geo (cclib.io.ccread): cclib object for the neutral ion at anion or cation geometry.

    For Psi4:

    opt_n (dict): Dictionary containing the optimized energy of the neutral population optimization analysis.
    sp_c (dict): Dictionary containing the single point energy of the vertical anion or cation.
    opt_c (dict): Dictionary containing the optimized energy of the anion or cation.
    n_c_geo (dict): Dictionary containing the single point energy of the neutral ion at anion or cation geometry.

    Returns
    -------
    float: The reorganization energy in eV.

    The reorganization energy is calculated using the following formula:
    reorg_en = (EcN - EnN) + (EnC - EcC)
        - EnN is the SCF energy of the neutral population optimization.
        - EcN is the SCF energy of the vertical anion or cation.
        - EcC is the SCF energy of the optimized anion or cation.
        - EnC is the SCF energy of the neutral ion charge at anion or cation geometry.
    """
    if calculation_software == 'Gaussian':

        EnN = opt_n.scfenergies[opt_n.optstatus==4][0]
        EcN = sp_c.scfenergies[0]
        EcC = opt_c.scfenergies[opt_c.optstatus==4][0]
        EnC = n_c_geo.scfenergies[0]

        reorg_en = (EcN-EnN)+(EnC-EcC)
    
    if calculation_software == 'Psi4':
        EnN = opt_n['optimized_energy']
        EcN = sp_c['sp_energy']
        EcC = opt_c['optimized_energy']
        EnC = n_c_geo['sp_energy']

        reorg_en = (EcN-EnN)+(EnC-EcC)

    return reorg_en

def cal_HOMO(opt):
    """
    Calculate the Highest Occupied Molecular Orbital (HOMO) energy.

    This function takes a cclib object representing the optimized neutral population and returns the HOMO energy in electron volts (eV).

    Parameters
    ----------
    opt (cclib.parser.data.ccData): A cclib object parsed from a .log file containing the optimization analysis of the neutral population.

    Returns
    -------
    float: The HOMO energy in electron volts (eV).
    """
    HOMO = opt.moenergies[0][opt.homos[0]]

    return HOMO

def cal_LUMO(opt):
    """
    Calculate the LUMO energy from the optimized neutral population.

    Parameters
    ----------
    opt (cclib.parser.data.ccData_optdone): The cclib object containing the parsed .log file for neutral population optimization analysis.

    Returns
    -------
    float: The LUMO energy in electron volts (eV).
    """
    LUMO = opt.moenergies[0][opt.homos[0]+1]

    return LUMO

def cal_gap(opt):
    """
    Calculate the HOMO-LUMO gap for an optimized neutral population.

    This function computes the energy difference between the Highest Occupied Molecular Orbital (HOMO) 
    and the Lowest Unoccupied Molecular Orbital (LUMO) in electron volts (eV).

    Parameters
    ----------
    opt (cclib.parser.ccData): A cclib object representing the parsed .log file for the 
                                     neutral population optimization analysis.

    Returns
    -------
    float: The energy difference between the HOMO and LUMO (HOMO-LUMO gap) in eV.
    """

    gap = cal_LUMO(opt) - cal_HOMO(opt)

    return gap

def cal_IP(opt_n, cation, IP_type):

    """
    Calculate the ionization potential (IP) given the cclib objects for the optimized neutral and cation populations.

    Parameters
    ----------
    opt_n (cclib.parser.data.ccData_optdone): cclib object for the optimized neutral population.
    cation (cclib.parser.data.ccData_optdone): cclib object for the optimized or vertical cation.
    IP_type (str): Type of ionization potential to calculate. Can be 'adiabatic' or 'vertical'.

    Returns
    -------
    float: The ionization potential (IP) calculated as the difference between the SCF energies of the optimized cation and neutral populations (eV).
    """

    EnN = opt_n.scfenergies[opt_n.optstatus==4][0]

    if IP_type == 'adiabatic':

        EcC = cation.scfenergies[cation.optstatus==4][0]

        IP = EcC - EnN

    if IP_type == 'vertical':

        EcN = cation.scfenergies[0]

        IP = EcN - EnN

    return IP

def cal_EA(opt_n, anion, EA_type):
    """
    Calculate the Electron Affinity (EA) given the optimized neutral population and optimized anion.
    
    Parameters
    ----------
    opt_n (cclib.io.ccread): Parsed cclib object for neutral population optimization analysis.
    anion (cclib.io.ccread): Parsed cclib object for optimized or vertical anion.
    EA_type (str): Type of electron affinity to calculate. Can be 'adiabatic' or 'vertical'.
    
    Returns
    -------
    float: The calculated Electron Affinity (EA).
    
    Notes
    ------
        - The function assumes that the SCF energies are available in the `scfenergies` attribute of the cclib objects.
        - The function uses the first SCF energy where the optimization status is 4 (indicating convergence).
    """

    EnN = opt_n.scfenergies[opt_n.optstatus==4][0]

    if EA_type == 'adiabatic':

        EaA = anion.scfenergies[anion.optstatus==4][0]

        EA = EnN - EaA

    if EA_type == 'vertical':
        
        EaN = anion.scfenergies[0]
    
        EA = EnN - EaN

    return EA

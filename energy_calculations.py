import cclib
from load_gaussian import *

def log_data(name,job):
    '''
    when providede the name of the dimer and the job name, the function fetches and parses the relevent files,
    and returns it as a cclib object

    name: str, name of the dimer (e.g.'a_2')
    job: str, the job type (e.g. 'ver_c')
    '''
    data = cclib.ccopen(f'{name}/{name}_{job}.log').parse()
    return data

def cal_reorg(pop_opt_n,ver,opt,neutral_ion_geo):
    """
    when provided with the cclib object for the corresponding files (using log_data() function),
    the function return the reorganisation energy.

    pop_opt_n: cclib obj, cclib parsed .log file for neutral population optimisation analaysis
    ver: cclib obj, cclib parsed .log file for vertical anion or cation
    opt: cclib obj, cclib parsed .log file for optimised anion or cation
    neutral_ion_geo: cclib obj, cclib parsed .log file for neutral ion geometry

    """
    EnN = pop_opt_n.scfenergies[pop_opt_n.optstatus==4][0]
    EcN = ver.scfenergies[0]
    EcC = opt.scfenergies[opt.optstatus==4][0]
    EnC = neutral_ion_geo.scfenergies[0]

    reorg_en = (EcN-EnN)+(EnC-EcC)

    return reorg_en

def get_HOMO(data):
    '''
    when provided a cclib object, returns the HOMO energy (eV)

    data: cclib object  
    '''
    HOMO = pop_opt_n.moenergies[0][pop_opt_n.homos[0]]
    return HOMO

def get_LUMO(data):
    '''
    when provided a cclib object, returns the LUMO energy (eV)

    data: cclib object 
    '''
    LUMO = pop_opt_n.moenergies[0][pop_opt_n.homos[0]+1]
    return LUMO

def calc_HL_gap(data):
    '''
    returns the energy difference between HOMO and LUMO (i.e. HOMO-LUMO gap) in eV.

    data: cclib object

    '''
    HL_gap = get_LUMO(data) - get_HOMO(data)

    return HL_gap

def calc_IP(pop_opt_n,ver_c,opt_c, kind):
    '''
    when provided the cclb object for the optimised neutral population,vertical cation, optimised cation
    and the type of the IP, function returns the ionisation energy.

    pop_opt_n: cclib obj, cclib parsed .log file for neutral population optimisation analaysis
    ver_c: cclib obj, cclib parsed .log file for vertical cation
    opt_c: cclib obj, cclib parsed .log file for optimised cation

    kind: str, type of IP:
        vertical = 'ver'
        adiabatic = 'adi'
    '''
    EnN = pop_opt_n.scfenergies[pop_opt_n.optstatus==4][0]

    if kind == 'ver':
        EcN = ver_c.scfenergies[0]
        IP = EcN - EnN

    if kind == 'adi':
        EcC = opt_c.scfenergies[opt_c.optstatus==4][0]
        IP = EcC - EnN

    return IP

def calc_EA(pop_opt_n, opt_a, an_geo_n, kind):
    '''
    when provided the optimised neutral population,optimised anion, neutral anion geometry and
    the type of the EA, function returns the energy for Electron Affinity.

    pop_opt_n: cclib obj, cclib parsed .log file for neutral population optimisation analaysis
    opt_a: cclib obj, cclib parsed .log file for optimised anion
    an_geo_n: cclib obj, cclib parsed .log file for neutral anion geometry

    kind: str, type of EA:
        vertical = 'ver'
        adiabatic = 'adi'
    '''
    EaA = opt_a.scfenergies[opt_a.optstatus==4][0]

    if kind == 'ver':
        EnA = an_geo_n.scfenergies[0]
        EA = EaA - EnA

    if kind == 'adi':
        EnN = pop_opt_n.scfenergies[pop_opt_n.optstatus==4][0]
        EA = EaA - EnN

    return EA

def calc_CSI_dimer(data, dimer_mol, a_mol, b_mol):
    '''
    when provided the dimer name, the optimised neutral population energy calculation, dimer dictionary,
    dictionary containing all melanin fragments and the dictionary of all organic electronic fragements,
    function returns the the charge seperation index between the fragments within the oligomer in eV.

    
    data: cclib obj, cclib parsed .log file for neutral population optimisation analaysis
    
    #dimer_name: str, name of the dimer in the dimer dictionary (e.g.'a_2')
    #dimer_dic: str, name of the .json file containing all dimers
    #mel_dic: str, name of the .json file containing all the melanin fragments
    #all_frag_dic: str, name of the .json file containing all the organic electronic fragements
    

    '''
    #dimer =  open_dictionary(dimer_dic)[dimer_name]   #when this is operational, the dictionary file should be mel_frag.json
    #dimer_mol = rdkit.Chem.MolFromSmiles(dimer.replace('{}',''))

    #a = open_dictionary(mel_dic)[dimer_name.split('_')[0]]  #this only works for frag-frag dimer, ideally 'a' should be changed to mel_frag2p.json
    #a_mol = rdkit.Chem.MolFromSmiles(a.replace('{}',''))

    #b = open_dictionary(all_frag_dic)[dimer_name.split('_')[1]]
    #b_mol = rdkit.Chem.MolFromSmiles(b.replace('{}',''))

    mulliken = data.atomcharges['mulliken']

    charge_a = mulliken[list(dimer_mol.GetSubstructMatch(a_mol))].sum()

    charge_b = mulliken[list(dimer_mol.GetSubstructMatch(b_mol))].sum()


    CSI = ((charge_a - charge_b)**2)**0.5

    return CSI

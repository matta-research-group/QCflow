import sys
sys.path.append('/Users/hongdingzhang/Year4/Fragment_library_code/github_main/QCflow_old/')
from torsion_parser import *
from energy_cal import*
import math
pi = math.pi

#to use these functions, make sure to first used the cclib_dic() function in cclib_parsing.py to create the relevant dictionaries.

def make_planarity_dic(tor_cclib_dic):
    """
    Function takes a cclib-parsed dimer dictionary of the torsion job and produces a dictionary where the keys are
    the dimer names and the values are the planarity values.

    tor_cclib_dic: dictionary
    """
    angles = []
    for dimer in tor_cclib_dic.values():
        angles.append(min_angle(dimer))

    angles_in_radian = []
    for angle in angles:
        angles_in_radian.append(abs(math.sin(angle * pi/180)))

    planarity_dic = {k:v for k,v in zip(dimer_dic.keys(),angles_in_radian)}

    return planarity_dic

def make_HOMO_dic(pop_opt_n_cclib_dic):
    """
    Function takes a cclib-parsed dimer dictionary of the population analysis job and produces a dictionary where
    the keys are the name of the dimers and the values are the HOMO values.

    pop_opt_n_cclib_dic: dictionary
    """
    HOMO_ls = []
    for obj in pop_opt_n_cclib_dic.vaues():
        HOMO_ls.append(cal_HOMO(obj))
    HOMO_dic = {k:v for k,v in zip(pop_opt_n_cclib_dic.keys(),HOMO_ls)}

    return HOMO_dic

def make_LUMO_dic(pop_opt_n_cclib_dic):
    """
    Function takes a cclib-parsed dimer dictionary of the population analysis job and produces a dictionary where
    the keys are the name of the dimers and the values are the LUMO values.

    pop_opt_n_cclib_dic: dictionary
    """
    LUMO_ls = []
    for obj in pop_opt_n_cclib_dic.values:
        LUMO_ls.append(cal_LUMO(obj))

    LUMO_dic = {k:v for k,v in zip(pop_opt_n_cclib_dic.keys(),LUMO_ls)}

    return LUMO_dic

def make_IP_dic(pop_opt_n_cclib_dic,ver_c_cclib_dic,opt_c_cclib_dic):
    """
    Function takes cclib-parsed dimer dictionaries of the population analysis job, ver_c job and opt_c job to produces
    a dictionary where the keys are the name of the dimers and the values are the LUMO values.

    pop_opt_n_cclib_dic: dictionary
    ver_c_cclib_dic: dictionary
    opt_c_cclib_dic: dictionary
    """
    IP_ls = []
    for pop_opt_n,ver_c,opt_c in zip(pop_opt_n_cclib_dic.values(),ver_c_cclib_dic.values(),opt_c_cclib_dic.values()):
        IP_ls.append(cal_IP(pop_opt_n,ver_c,opt_c,'adi'))

    IP_dic = {k:v for k,v in zip(pop_opt_n_cclib_dic.keys(),IP_ls)}

    return IP_dic

def make_EA_dic(pop_opt_n_cclib_dic,opt_a_cclib_dic,n_a_geo_cclib_dic):
    """
    Function takes cclib-parsed dimer dictionaries of the population analysis job, ver_c job and opt_c job to produces
    a dictionary where the keys are the name of the dimers and the values are the LUMO values.

    pop_opt_n_cclib_dic: dictionary
    ver_c_cclib_dic: dictionary
    opt_c_cclib_dic: dictionary
    """
    EA_ls = []
    for pop_opt_n,opt_a,an_geo_n in zip(pop_opt_n_cclib_dic.values(),opt_a_cclib_dic.values(),n_a_geo_cclib_dic.values()):
        EA.append(cal_EA(pop_opt_n,opt_a,an_geo_n,'adi'))

    EA_dic = {k:v for k,v in zip(dimer_dic.keys(),EA_ls)}

    return EA_dic

def make_cat_reorg_dic(pop_opt_n_cclib_dic,ver_c_cclib_dic,opt_c_cclib_dic,n_c_geo_cclib_dic):
    """
    Function takes cclib-parsed dimer dictionaries that are relavent to calculate cationic reorganisation energy
    and produces a dictionary where the keys are name of the dimers and the values are the cationic reorganisation values.

    pop_opt_n_cclib_dic: dictionary
    ver_c_cclib_dic: dictionary
    opt_c_cclib_dic: dictionary
    n_c_geo_cclib_dic: dictionary
    """
    cat_reorg = []

    for pop_opt_n,ver,opt,neutral_ion_geo in zip(pop_opt_n_cclib_dic.values(),ver_c_cclib_dic.values(),
                                             opt_c_cclib_dic.values(),n_c_geo_cclib_dic.values()):
        cat_reorg.append(cal_reorg(pop_opt_n,ver,opt,neutral_ion_geo))

    cat_reorg_dic = {k:v for k,v in zip(dimer_dic.keys(),cat_reorg)}

    return cat_reorg_dic

def make_an_reorg_dic(pop_opt_n_cclib_dic,ver_a_cclib_dic,opt_a_cclib_dic,n_a_geo_cclib_dic):
    """
    Function takes cclib-parsed dimer dictionaries that are relavent to calculate cationic reorganisation energy
    and produces a dictionary where the keys are name of the dimers and the values are the cationic reorganisation values.

    pop_opt_n_cclib_dic: dictionary
    ver_a_cclib_dic: dictionary
    opt_a_cclib_dic: dictionary
    n_a_geo_cclib_dic: dictionary
    """
    an_reorg = []
    for pop_opt_n,ver,opt,neutral_ion_geo in zip(pop_opt_n_cclib_dic.values(),ver_a_cclib_dic.values(),
                                             opt_a_cclib_dic.values(),n_a_geo_cclib_dic.values()):
        an_reorg.append(cal_reorg(pop_opt_n,ver,opt,neutral_ion_geo))

    an_reorg_dic = {k:v for k,v in zip(dimer_dic.keys(),an_reorg)}

    return an_reorg_dic

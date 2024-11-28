# This python script is a example script that was used to extract data on a large scale (100 molecule or more)
# This python script was submited on the KCL CREATE HPC and extracted the energy gap and planarity data
# The json files for the dictionaries are not provided but the example placeholder names are shown

import QCflow
from QCflow.load_gaussian import *
from QCflow.energy_calculations import *
from QCflow.torsion_parser import *

dyads = open_dictionary('low_predicted_energy_gap_dyads.json') #load dictionary of the calculated dyads, key is the name and value is the SMILES string

opt_data = {}
for k, v in dyads.items():
    opt = load_data(k, 'opt')
    opt_data[k] = opt

failed = {}
for k, v in opt_data.items():
    if v.metadata['success'] == False:
        failed[k] = dyads[k]

save_dictionary(failed, 'failed_dyads.json') #stores the failed systems so they can be manually checked

passed_dyads = {k: v for k, v in opt_data.items() if k not in failed} #filters out the failed dyads and creates a new dictionary with only the passed dyads

passed_dyads_smi = {}
for k, v in passed_dyads.items():
    passed_dyads_smi[k] = dyads[k]

save_dictionary(passed_dyads_smi, 'passed_dyads_smi.json') #stores the passed systems

linker_type = {} #creates a dictionary to store the linker type
for k, v in passed_dyads_smi.items():
    linker_name = k.split('_')[2]
    linker_type[k] = linker_name

planarity_score = {} #calculates planarity based on the linker used in the dyad
for k, v in passed_dyads.items():
    if linker_type[k] == 'single':
        di_angle = finding_dihedral_opt(k, passed_dyads_smi[k])
        planarity_score[k] = find_planarity(di_angle)

    if linker_type[k] != 'single':
        planrity_multi = finding_multi_planairty(k, passed_dyads_smi[k], linker_type[k])
        planarity_score[k] = planrity_multi

save_dictionary(planarity_score, 'planarity_score_dyads.json') #stores the planarity score

eg_data = {} #calculates the energy gap of the dyads
for k, v in opt_data.items():
    eg = cal_gap(v)
    eg_data[k] = eg

save_dictionary(eg_data, 'eg_data.json') #stores the energy gap
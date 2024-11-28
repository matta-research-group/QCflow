# Calculations are submitted in bulk to the KCL CREATE HPC
# QCflow functionality alows for this to be done in a single python script and can submit large quantities of calculations
# The json files for the dictionaries are not provided but the example placeholder names are shown

import QCflow
from QCflow.load_gaussian import *
from QCflow.run_gaussian import *

dyads = open_dictionary('low_predicted_energy_gap_dyads.json') #key is the name, value is the SMILES string

for k, v in dyads.items():
    run_calc('opt', k, v) #optimisation calculation

#if I wanted to run the initial stages of a reorganisation energy calculation it would look like this:

for k, v in dyads.items():
    run_calc('opt_c', k, v) #optimisation calculation at cationic charge
    run_calc('opt_a', k, v) #optimisation calculation at anioinc charge
    run_calc('sp_a', k, v) #single point calculation at anioinc charge, neutral geometry
    run_calc('sp_c', k, v) #single point calculation at cationic charge, neutral geometry
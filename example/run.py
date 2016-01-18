import numpy as np
import cobra.io

model = cobra.io.load_json_model('ecoli_core_kegg.json')

external_indicies = [
    12, # Biomass
    27, # Glucose
    # 35, # Oxygen
    22, # CO2
    34, # Ammonium
    # 36, # Phosphate
]

y0 = [
    0.1,  # Biomass
    40.,  # Glucose
    # 100., # Oxygen
    0,    # CO2
    50.,   # Ammonium
    # 50.,   # Phosphate
]

parameters = np.array([5, 2])

from dFBA import DFBA_Simulator
from rbf import GlucoseUptake

dfbasimulator = DFBA_Simulator(model, 
                               np.array(external_indicies, dtype=np.int32),
                               np.array(y0),
                               GlucoseUptake(
                                   np.array(parameters, dtype=np.float)),
                               collect_fluxes=True)

# This command should generate a python warning, and issue a CVodes Error
# indicating the function evaluation failed. This is due to the fact that
# biomass growth is not possible once glucose concetation falls to a critical
# threshold
dfbasimulator.integrate(t_end=10., res=20)

print dfbasimulator.ys[:,:2]

dfbasimulator.print_final_stats()



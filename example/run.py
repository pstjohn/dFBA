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

parameters = np.array([10, 0.5])

from dFBA import DFBA_Simulator
from rbf import GlucoseUptake

dfbasimulator = DFBA_Simulator(model, 
                               np.array(external_indicies, dtype=np.int32),
                               np.array(y0),
                               GlucoseUptake(parameters))

dfbasimulator.integrate(10., 20)

print dfbasimulator.ys[:,:2]

dfbasimulator.print_final_stats()



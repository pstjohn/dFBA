import numpy as np
from Models import ecoli_core_kegg

model = ecoli_core_kegg()

# external_indicies = np.array(
#     [model.reactions.index('Biomass_Ecoli_core_w_GAM')] +
#     [model.reactions.index(r.id) for r in
#      model.reactions.iquery('system_boundary', 'boundary')], dtype=np.int32)

# {r.id : model.reactions.index(r.id) for r in
#  model.reactions.iquery('system_boundary', 'boundary')}

external_indicies = [
    12, # Biomass
    27, # Glucose
    35, # Oxygen
    22, # CO2
    34, # Ammonium
    36, # Phosphate
]

y0 = [
    0.1,  # Biomass
    40.,  # Glucose
    100., # Oxygen
    0,    # CO2
    50.,   # Ammonium
    50.,   # Phosphate
]


import dfbasim
from rbf import GlucoseUptake
# cvode.main(model, np.array(external_indicies, dtype=np.int32), np.array(y0))
dfbasimulator = dfbasim.DFBA_Simulator(model, 
                               np.array(external_indicies, dtype=np.int32),
                               np.array(y0),
                               GlucoseUptake())

dfbasimulator.integrate(10., 20)

print dfbasimulator.ys[:,:2]

dfbasimulator.print_final_stats()



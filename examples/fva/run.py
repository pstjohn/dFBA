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
    0.05,  # Biomass
    40.,  # Glucose
    # 100., # Oxygen
    0,    # CO2
    50.,   # Ammonium
    # 50.,   # Phosphate
]

parameters = np.array([5, 2], dtype=np.float)

from rbf import GlucoseUptake
from dFBA.sim.dfvasim import DFVA_Simulator

dfva = DFVA_Simulator(model, external_indicies, y0, GlucoseUptake(parameters),
                      .99)


dfva.update_bounds_func_parameters(parameters)
dfva.set_death_rate(0.01)
dfva.integrate(t_end=10., res=20)

print dfva.ys_max[:,-1] - dfva.ys_min[:,-1]
print dfva.ys_dot_max[:,-1] - dfva.ys_dot_min[:,-1]

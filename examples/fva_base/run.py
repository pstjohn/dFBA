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

c_inner = np.array([r.objective_coefficient for r in model.reactions])

rxn_ind = model.reactions.index('EX_nh4_e')

c_outer = np.zeros(len(model.reactions))
c_outer[rxn_ind] = +1


from dFBA.sim.dfvasim_base import DFVA_Simulator_base
from rbf import GlucoseUptake

dfvasimulator_ub = DFVA_Simulator_base(
    model, np.array(external_indicies, dtype=np.int32),
    np.array(y0), GlucoseUptake(np.array(parameters, dtype=np.float)),
    c_inner, c_outer, .99)

dfvasimulator_lb = DFVA_Simulator_base(
    model, np.array(external_indicies, dtype=np.int32),
    np.array(y0), GlucoseUptake(np.array(parameters, dtype=np.float)),
    c_inner, -c_outer, .99)

# This command should generate a python warning, and issue a CVodes Error
# indicating the function evaluation failed. This is due to the fact that
# biomass growth is not possible once glucose concetation falls to a critical
# threshold
dfvasimulator_ub.integrate(t_end=10., res=20)
dfvasimulator_lb.integrate(t_end=10., res=20)

print dfvasimulator_ub.ys[:,3] - dfvasimulator_lb.ys[:,3]

dfvasimulator_lb.print_final_stats()



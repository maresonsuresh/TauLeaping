import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create clarke1 model
def create_clarke1(parameter_values=None):
    model = gillespy2.Model(name="clarke1")
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    Pi = gillespy2.Species(name="Pi", initial_value=0.0, mode="discrete")
    R_smad_P_cyt = gillespy2.Species(name="R_smad_P_cyt", initial_value=0.0, mode="discrete")
    R_smad_P_nuc = gillespy2.Species(name="R_smad_P_nuc", initial_value=0.0, mode="discrete")
    R_smad_P_smad4_cyt = gillespy2.Species(name="R_smad_P_smad4_cyt", initial_value=0.0, mode="discrete")
    R_smad_P_smad4_nuc = gillespy2.Species(name="R_smad_P_smad4_nuc", initial_value=0.0, mode="discrete")
    R_smad_cyt = gillespy2.Species(name="R_smad_cyt", initial_value=162000.0, mode="discrete")
    R_smad_nuc = gillespy2.Species(name="R_smad_nuc", initial_value=18000.0, mode="discrete")
    receptor = gillespy2.Species(name="receptor", initial_value=10000.0, mode="discrete")
    smad4_cyt = gillespy2.Species(name="smad4_cyt", initial_value=120000.0, mode="discrete")
    smad4_nuc = gillespy2.Species(name="smad4_nuc", initial_value=30000.0, mode="discrete")
    model.add_species([Pi, R_smad_P_cyt, R_smad_P_nuc, R_smad_P_smad4_cyt, R_smad_P_smad4_nuc, R_smad_cyt, R_smad_nuc, receptor, smad4_cyt, smad4_nuc])

    # Rate parameters
    KCAT = gillespy2.Parameter(name="KCAT", expression=3.51)  # min^-1
    K1 = gillespy2.Parameter(name="K1", expression=289000.0)  # substance
    k5nc = gillespy2.Parameter(name="k5nc", expression=5.63)  # min^-1
    k5cn = gillespy2.Parameter(name="k5cn", expression=0.563)  # min^-1
    k4nc = gillespy2.Parameter(name="k4nc", expression=0.783)  # min^-1
    k4cn = gillespy2.Parameter(name="k4cn", expression=0.00497)  # min^-1
    k2a = gillespy2.Parameter(name="k2a", expression=0.000065)  # dimensionless
    k2d = gillespy2.Parameter(name="k2d", expression=0.0399)  # min^-1
    k3 = gillespy2.Parameter(name="k3", expression=16.6)  # min^-1
    k6d = gillespy2.Parameter(name="k6d", expression=0.0492)  # min^-1
    k6a = gillespy2.Parameter(name="k6a", expression=0.000144)  # dimensionless
    Vmax7 = gillespy2.Parameter(name="Vmax7", expression=17100.0)  # items/min
    K7 = gillespy2.Parameter(name="K7", expression=8950.0)  # substance
    sum_R_smad_cyt = gillespy2.Parameter(name="sum_R_smad_cyt", expression=0.0)
    sum_R_smad_nuc = gillespy2.Parameter(name="sum_R_smad_nuc", expression=0.0)
    sum_smad4_cyt = gillespy2.Parameter(name="sum_smad4_cyt", expression=0.0)
    sum_smad4_nuc = gillespy2.Parameter(name="sum_smad4_nuc", expression=0.0)
    
    # Add Parameters to Model
    model.add_parameter([
        KCAT, K1, k5nc, k5cn, k4nc, k4cn, k2a, k2d, k3, k6d, k6a, Vmax7, K7,
        sum_R_smad_cyt, sum_R_smad_nuc, sum_smad4_cyt, sum_smad4_nuc
    ])  

    # Receptor degradation
    reaction_0 = gillespy2.Reaction(
        name="reaction_0",
        reactants={'receptor': 1},
        products={},
        propensity_function="50"
    )

    # Phosphorylation
    reaction_1 = gillespy2.Reaction(
        name="reaction_1",
        reactants={'R_smad_cyt': 1},
        products={'R_smad_P_cyt': 1},
        propensity_function="KCAT * receptor * R_smad_cyt / (K1 + R_smad_cyt)"
    )

    # Complex formation (reversible)
    reaction_2 = gillespy2.Reaction(
        name="reaction_2",
        reactants={'R_smad_P_cyt': 1, 'smad4_cyt': 1},
        products={'R_smad_P_smad4_cyt': 1},
        propensity_function="k2a * R_smad_P_cyt * smad4_cyt - k2d * R_smad_P_smad4_cyt"
    )

    # Complex translocation
    reaction_3 = gillespy2.Reaction(
        name="reaction_3",
        reactants={'R_smad_P_smad4_cyt': 1},
        products={'R_smad_P_smad4_nuc': 1},
        propensity_function="k3 * R_smad_P_smad4_cyt"
    )

    # Smad4 translocation (reversible)
    reaction_4 = gillespy2.Reaction(
        name="reaction_4",
        reactants={'smad4_nuc': 1},
        products={'smad4_cyt': 1},
        propensity_function="k4nc * smad4_nuc - k4cn * smad4_cyt"
    )

    # R-Smad translocation (reversible)
    reaction_5 = gillespy2.Reaction(
        name="reaction_5",
        reactants={'R_smad_nuc': 1},
        products={'R_smad_cyt': 1},
        propensity_function="k5nc * R_smad_nuc - k5cn * R_smad_cyt"
    )

    # Complex dissociation in nucleus (reversible)
    reaction_6 = gillespy2.Reaction(
        name="reaction_6",
        reactants={'R_smad_P_smad4_nuc': 1},
        products={'smad4_nuc': 1, 'R_smad_P_nuc': 1},
        propensity_function="k6d * R_smad_P_smad4_nuc - k6a * smad4_nuc * R_smad_P_nuc"
    )

    # Dephosphorylation
    reaction_7 = gillespy2.Reaction(
        name="reaction_7",
        reactants={'R_smad_P_nuc': 1},
        products={'R_smad_nuc': 1, 'Pi': 1},
        propensity_function="Vmax7 * R_smad_P_nuc / (K7 + R_smad_P_nuc)"
    )

    # Add Reactions to Model
    model.add_reaction([reaction_0, reaction_1, reaction_2, reaction_3, reaction_4, reaction_5, reaction_6, reaction_7])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=100)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_clarke1()

start_time = time.time()
results = model.run(algorithm="SSA")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

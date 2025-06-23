import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create hatakeyama1 model
def create_hatakeyama1(parameter_values=None):
    model = gillespy2.Model(name="hatakeyama1")
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    Akt = gillespy2.Species(name="Akt", initial_value=10.0, mode="discrete")
    AktPIP = gillespy2.Species(name="AktPIP", initial_value=0.0, mode="discrete")
    AktPIP3 = gillespy2.Species(name="AktPIP3", initial_value=0.0, mode="discrete")
    AktPIPP = gillespy2.Species(name="AktPIPP", initial_value=0.0, mode="discrete")
    E = gillespy2.Species(name="E", initial_value=7.0, mode="discrete")
    ERK = gillespy2.Species(name="ERK", initial_value=1000.0, mode="discrete")
    ERKP = gillespy2.Species(name="ERKP", initial_value=0.0, mode="discrete")
    ERKPP = gillespy2.Species(name="ERKPP", initial_value=0.0, mode="discrete")
    GS = gillespy2.Species(name="GS", initial_value=10.0, mode="discrete")
    HRG = gillespy2.Species(name="HRG", initial_value=330.0, mode="discrete")
    MEK = gillespy2.Species(name="MEK", initial_value=120.0, mode="discrete")
    MEKP = gillespy2.Species(name="MEKP", initial_value=0.0, mode="discrete")
    MEKPP = gillespy2.Species(name="MEKPP", initial_value=0.0, mode="discrete")
    MKP3 = gillespy2.Species(name="MKP3", initial_value=2.0, mode="discrete")
    PI3K = gillespy2.Species(name="PI3K", initial_value=10.0, mode="discrete")
    PI3Kstar = gillespy2.Species(name="PI3Kstar", initial_value=0.0, mode="discrete")
    PIP3 = gillespy2.Species(name="PIP3", initial_value=0.0, mode="discrete")
    PP2A = gillespy2.Species(name="PP2A", initial_value=11.0, mode="discrete")
    P_I = gillespy2.Species(name="P_I", initial_value=800.0, mode="discrete")
    R = gillespy2.Species(name="R", initial_value=80.0, mode="discrete")
    RHRG = gillespy2.Species(name="RHRG", initial_value=0.0, mode="discrete")
    RHRG2 = gillespy2.Species(name="RHRG2", initial_value=0.0, mode="discrete")
    RP = gillespy2.Species(name="RP", initial_value=0.0, mode="discrete")
    RPI3K = gillespy2.Species(name="RPI3K", initial_value=0.0, mode="discrete")
    RPI3Kstar = gillespy2.Species(name="RPI3Kstar", initial_value=0.0, mode="discrete")
    RShGS = gillespy2.Species(name="RShGS", initial_value=0.0, mode="discrete")
    RShP = gillespy2.Species(name="RShP", initial_value=0.0, mode="discrete")
    RShc = gillespy2.Species(name="RShc", initial_value=0.0, mode="discrete")
    Raf = gillespy2.Species(name="Raf", initial_value=100.0, mode="discrete")
    Rafstar = gillespy2.Species(name="Rafstar", initial_value=0.0, mode="discrete")
    RasGDP = gillespy2.Species(name="RasGDP", initial_value=120.0, mode="discrete")
    RasGTP = gillespy2.Species(name="RasGTP", initial_value=0.0, mode="discrete")
    ShGS = gillespy2.Species(name="ShGS", initial_value=0.0, mode="discrete")
    ShP = gillespy2.Species(name="ShP", initial_value=0.0, mode="discrete")
    Shc = gillespy2.Species(name="Shc", initial_value=1000.0, mode="discrete")
    internalization = gillespy2.Species(name="internalization", initial_value=0.0, mode="discrete")

    # Add all species to the model
    model.add_species([
        Akt, AktPIP, AktPIP3, AktPIPP, E, ERK, ERKP, ERKPP, GS, HRG,
        MEK, MEKP, MEKPP, MKP3, PI3K, PI3Kstar, PIP3, PP2A, P_I, R,
        RHRG, RHRG2, RP, RPI3K, RPI3Kstar, RShGS, RShP, RShc, Raf, Rafstar,
        RasGDP, RasGTP, ShGS, ShP, Shc, internalization
    ])

    # Define Global Parameters
    AktPP_percent = gillespy2.Parameter(name="AktPP_percent", expression=0.0)
    ERKPP_percent = gillespy2.Parameter(name="ERKPP_percent", expression=0.0)
    K10 = gillespy2.Parameter(name="K10", expression=340.0)
    K11 = gillespy2.Parameter(name="K11", expression=0.181)
    K12 = gillespy2.Parameter(name="K12", expression=0.0571)
    K13 = gillespy2.Parameter(name="K13", expression=11.7)
    K14 = gillespy2.Parameter(name="K14", expression=8.07)
    K15 = gillespy2.Parameter(name="K15", expression=317.0)
    K16 = gillespy2.Parameter(name="K16", expression=2200.0)
    K17 = gillespy2.Parameter(name="K17", expression=317.0)
    K18 = gillespy2.Parameter(name="K18", expression=60.0)
    K19 = gillespy2.Parameter(name="K19", expression=146000.0)
    K20 = gillespy2.Parameter(name="K20", expression=160.0)
    K21 = gillespy2.Parameter(name="K21", expression=146000.0)
    K22 = gillespy2.Parameter(name="K22", expression=60.0)
    K26 = gillespy2.Parameter(name="K26", expression=3680.0)
    K27 = gillespy2.Parameter(name="K27", expression=39.1)
    K28 = gillespy2.Parameter(name="K28", expression=9.02)
    K30 = gillespy2.Parameter(name="K30", expression=80000.0)
    K31 = gillespy2.Parameter(name="K31", expression=4.35)
    K32 = gillespy2.Parameter(name="K32", expression=80000.0)
    K33 = gillespy2.Parameter(name="K33", expression=12.0)
    K4 = gillespy2.Parameter(name="K4", expression=50.0)
    MEKPP_percent = gillespy2.Parameter(name="MEKPP_percent", expression=0.0)
    PI3Kstar_percent = gillespy2.Parameter(name="PI3Kstar_percent", expression=0.0)
    RP_percent = gillespy2.Parameter(name="RP_percent", expression=0.0)
    Rafstar_percent = gillespy2.Parameter(name="Rafstar_percent", expression=0.0)
    ShP_percent = gillespy2.Parameter(name="ShP_percent", expression=0.0)
    V10 = gillespy2.Parameter(name="V10", expression=0.0154)
    V12 = gillespy2.Parameter(name="V12", expression=0.289)
    V26 = gillespy2.Parameter(name="V26", expression=2620.0)
    V28 = gillespy2.Parameter(name="V28", expression=17000.0)
    V30 = gillespy2.Parameter(name="V30", expression=20000.0)
    V32 = gillespy2.Parameter(name="V32", expression=20000.0)
    V4 = gillespy2.Parameter(name="V4", expression=62.5)
    k1 = gillespy2.Parameter(name="k1", expression=0.0012)
    k11 = gillespy2.Parameter(name="k11", expression=0.222)
    k13 = gillespy2.Parameter(name="k13", expression=1.53)
    k14 = gillespy2.Parameter(name="k14", expression=0.00673)
    k15 = gillespy2.Parameter(name="k15", expression=3.5)
    k16 = gillespy2.Parameter(name="k16", expression=0.058)
    k17 = gillespy2.Parameter(name="k17", expression=2.9)
    k18 = gillespy2.Parameter(name="k18", expression=0.058)
    k19 = gillespy2.Parameter(name="k19", expression=9.5)
    k2 = gillespy2.Parameter(name="k2", expression=0.01)
    k20 = gillespy2.Parameter(name="k20", expression=0.3)
    k21 = gillespy2.Parameter(name="k21", expression=16.0)
    k22 = gillespy2.Parameter(name="k22", expression=0.27)
    k23 = gillespy2.Parameter(name="k23", expression=0.1)
    k24 = gillespy2.Parameter(name="k24", expression=9.85)
    k25 = gillespy2.Parameter(name="k25", expression=45.8)
    k27 = gillespy2.Parameter(name="k27", expression=16.9)
    k29 = gillespy2.Parameter(name="k29", expression=507.0)
    k3 = gillespy2.Parameter(name="k3", expression=1.0)
    k31 = gillespy2.Parameter(name="k31", expression=0.107)
    k33 = gillespy2.Parameter(name="k33", expression=0.211)
    k34 = gillespy2.Parameter(name="k34", expression=0.001)
    k5 = gillespy2.Parameter(name="k5", expression=0.1)
    k6 = gillespy2.Parameter(name="k6", expression=20.0)
    k7 = gillespy2.Parameter(name="k7", expression=60.0)
    k8 = gillespy2.Parameter(name="k8", expression=2040.0)
    k9 = gillespy2.Parameter(name="k9", expression=40.8)
    k_1 = gillespy2.Parameter(name="k_1", expression=0.00076)
    k_2 = gillespy2.Parameter(name="k_2", expression=0.1)
    k_23 = gillespy2.Parameter(name="k_23", expression=2.0)
    k_24 = gillespy2.Parameter(name="k_24", expression=0.0985)
    k_25 = gillespy2.Parameter(name="k_25", expression=0.047)
    k_29 = gillespy2.Parameter(name="k_29", expression=234.0)
    k_3 = gillespy2.Parameter(name="k_3", expression=0.01)
    k_34 = gillespy2.Parameter(name="k_34", expression=0.0)
    k_5 = gillespy2.Parameter(name="k_5", expression=1.0)
    k_6 = gillespy2.Parameter(name="k_6", expression=5.0)
    k_7 = gillespy2.Parameter(name="k_7", expression=546.0)
    k_8 = gillespy2.Parameter(name="k_8", expression=15700.0)
    k_9 = gillespy2.Parameter(name="k_9", expression=0.0)

    # Add all parameters to the model
    model.add_parameter([
        AktPP_percent, ERKPP_percent, K10, K11, K12, K13, K14, K15, K16, K17,
        K18, K19, K20, K21, K22, K26, K27, K28, K30, K31,
        K32, K33, K4, MEKPP_percent, PI3Kstar_percent, RP_percent, Rafstar_percent, ShP_percent,
        V10, V12, V26, V28, V30, V32, V4, k1, k11, k13,
        k14, k15, k16, k17, k18, k19, k2, k20, k21, k22,
        k23, k24, k25, k27, k29, k3, k31, k33, k34, k5,
        k6, k7, k8, k9, k_1, k_2, k_23, k_24, k_25, k_29,
        k_3, k_34, k_5, k_6, k_7, k_8, k_9
    ])

    # Define Reactions
    reaction_1 = gillespy2.Reaction(
        name="reaction_1",
        reactants={'R': 1, 'HRG': 1},
        products={'RHRG': 1},
        propensity_function="(k1 * R * HRG - k_1 * RHRG)"
    )

    reaction_2 = gillespy2.Reaction(
        name="reaction_2",
        reactants={'RHRG': 2},
        products={'RHRG2': 1},
        propensity_function="(k2 * pow(RHRG, 2) - k_2 * RHRG2)"
    )

    reaction_3 = gillespy2.Reaction(
        name="reaction_3",
        reactants={'RHRG2': 1},
        products={'RP': 1},
        propensity_function="(k3 * RHRG2 - k_3 * RP)"
    )

    reaction_4 = gillespy2.Reaction(
        name="reaction_4",
        reactants={'RP': 1},
        products={'RHRG2': 1},
        propensity_function="V4 * RP / (K4 + RP)"
    )

    reaction_5 = gillespy2.Reaction(
        name="reaction_5",
        reactants={'RP': 1, 'Shc': 1},
        products={'RShc': 1},
        propensity_function="(k5 * RP * Shc - k_5 * RShc)"
    )

    reaction_6 = gillespy2.Reaction(
        name="reaction_6",
        reactants={'RShc': 1},
        products={'RShP': 1},
        propensity_function="(k6 * RShc - k_6 * RShP)"
    )

    reaction_7 = gillespy2.Reaction(
        name="reaction_7",
        reactants={'RShP': 1, 'GS': 1},
        products={'RShGS': 1},
        propensity_function="(k7 * RShP * GS - k_7 * RShGS)"
    )

    reaction_8 = gillespy2.Reaction(
        name="reaction_8",
        reactants={'RShGS': 1},
        products={'ShGS': 1, 'RP': 1},
        propensity_function="(k8 * RShGS - k_8 * ShGS * RP)"
    )

    reaction_9 = gillespy2.Reaction(
        name="reaction_9",
        reactants={'ShGS': 1},
        products={'GS': 1, 'ShP': 1},
        propensity_function="(k9 * ShGS - k_9 * GS * ShP)"
    )

    reaction_10 = gillespy2.Reaction(
        name="reaction_10",
        reactants={'ShP': 1},
        products={'Shc': 1},
        propensity_function="V10 * ShP / (K10 + ShP)"
    )

    reaction_11 = gillespy2.Reaction(
        name="reaction_11",
        reactants={'RasGDP': 1},
        products={'RasGTP': 1},
        propensity_function="k11 * ShGS * RasGDP / (K11 + RasGDP)"
    )

    reaction_12 = gillespy2.Reaction(
        name="reaction_12",
        reactants={'RasGTP': 1},
        products={'RasGDP': 1},
        propensity_function="V12 * RasGTP / (K12 + RasGTP)"
    )

    reaction_13 = gillespy2.Reaction(
        name="reaction_13",
        reactants={'Raf': 1},
        products={'Rafstar': 1},
        propensity_function="k13 * RasGTP * Raf / (K13 + Raf)"
    )

    reaction_14 = gillespy2.Reaction(
        name="reaction_14",
        reactants={'Rafstar': 1},
        products={'Raf': 1},
        propensity_function="k14 * (AktPIPP + E) * Rafstar / (K14 + Rafstar)"
    )

    reaction_15 = gillespy2.Reaction(
        name="reaction_15",
        reactants={'MEK': 1},
        products={'MEKP': 1},
        propensity_function="k15 * Rafstar * MEK / (K15 * (1 + MEKP / K17) + MEK)"
    )

    reaction_16 = gillespy2.Reaction(
        name="reaction_16",
        reactants={'MEKP': 1},
        products={'MEK': 1},
        propensity_function="k16 * PP2A * MEKP / (K16 * (1 + MEKPP / K18 + AktPIP / K31 + AktPIPP / K33) + MEKP)"
    )

    reaction_17 = gillespy2.Reaction(
        name="reaction_17",
        reactants={'MEKP': 1},
        products={'MEKPP': 1},
        propensity_function="k17 * Rafstar * MEKP / (K17 * (1 + MEK / K15) + MEKP)"
    )

    reaction_18 = gillespy2.Reaction(
        name="reaction_18",
        reactants={'MEKPP': 1},
        products={'MEKP': 1},
        propensity_function="k18 * PP2A * MEKPP / (K18 * (1 + MEKP / K16 + AktPIPP / K31 + AktPIPP / K33) + MEKPP)"
    )

    reaction_19 = gillespy2.Reaction(
        name="reaction_19",
        reactants={'ERK': 1},
        products={'ERKP': 1},
        propensity_function="k19 * MEKPP * ERK / (K19 * (1 + ERKP / K21) + ERK)"
    )

    reaction_20 = gillespy2.Reaction(
        name="reaction_20",
        reactants={'ERKP': 1},
        products={'ERK': 1},
        propensity_function="k20 * MKP3 * ERKP / (K20 * (1 + ERKPP / K22) + ERKP)"
    )

    reaction_21 = gillespy2.Reaction(
        name="reaction_21",
        reactants={'ERKP': 1},
        products={'ERKPP': 1},
        propensity_function="k21 * MEKPP * ERKP / (K21 * (1 + ERK / K19) + ERKP)"
    )

    reaction_22 = gillespy2.Reaction(
        name="reaction_22",
        reactants={'ERKPP': 1},
        products={'ERKP': 1},
        propensity_function="k22 * MKP3 * ERKPP / (K22 * (1 + ERKP / K20) + ERKPP)"
    )

    reaction_23 = gillespy2.Reaction(
        name="reaction_23",
        reactants={'RP': 1, 'PI3K': 1},
        products={'RPI3K': 1},
        propensity_function="(k23 * RP * PI3K - k_23 * RPI3K)"
    )

    reaction_24 = gillespy2.Reaction(
        name="reaction_24",
        reactants={'RPI3K': 1},
        products={'RPI3Kstar': 1},
        propensity_function="(k24 * RPI3K - k_24 * RPI3Kstar)"
    )

    reaction_25 = gillespy2.Reaction(
        name="reaction_25",
        reactants={'RPI3Kstar': 1},
        products={'RP': 1, 'PI3Kstar': 1},
        propensity_function="(k25 * RPI3Kstar - k_25 * RP * PI3Kstar)"
    )

    reaction_26 = gillespy2.Reaction(
        name="reaction_26",
        reactants={'PI3Kstar': 1},
        products={'PI3K': 1},
        propensity_function="V26 * PI3Kstar / (K26 + PI3Kstar)"
    )

    reaction_27 = gillespy2.Reaction(
        name="reaction_27",
        reactants={'P_I': 1},
        products={'PIP3': 1},
        propensity_function="k27 * PI3Kstar * P_I / (K27 + P_I)"
    )

    reaction_28 = gillespy2.Reaction(
        name="reaction_28",
        reactants={'PIP3': 1},
        products={'P_I': 1},
        propensity_function="V28 * PIP3 / (K28 + PIP3)"
    )

    reaction_29 = gillespy2.Reaction(
        name="reaction_29",
        reactants={'PIP3': 1, 'Akt': 1},
        products={'AktPIP3': 1},
        propensity_function="(k29 * PIP3 * Akt - k_29 * AktPIP3)"
    )

    reaction_30 = gillespy2.Reaction(
        name="reaction_30",
        reactants={'AktPIP3': 1},
        products={'AktPIP': 1},
        propensity_function="V30 * AktPIP3 / (K30 * (1 + AktPIP / K32) + AktPIP3)"
    )

    reaction_31 = gillespy2.Reaction(
        name="reaction_31",
        reactants={'AktPIP': 1},
        products={'AktPIP3': 1},
        propensity_function="k31 * PP2A * AktPIP / (K31 * (1 + MEKP / K16 + MEKPP / K18 + AktPIPP / K33) + AktPIP)"
    )

    reaction_32 = gillespy2.Reaction(
        name="reaction_32",
        reactants={'AktPIP': 1},
        products={'AktPIPP': 1},
        propensity_function="V32 * AktPIP / (K32 * (1 + AktPIP3 / K30) + AktPIP)"
    )

    reaction_33 = gillespy2.Reaction(
        name="reaction_33",
        reactants={'AktPIPP': 1},
        products={'AktPIP': 1},
        propensity_function="k33 * PP2A * AktPIPP / (K33 * (1 + MEKP / K16 + MEKPP / K18 + AktPIP / K31) + AktPIPP)"
    )

    reaction_34 = gillespy2.Reaction(
        name="reaction_34",
        reactants={'RP': 1},
        products={'internalization': 1},
        propensity_function="(k34 * RP - k_34 * internalization)"
    )

    # Add all reactions to model
    model.add_reaction([
        reaction_1, reaction_2, reaction_3, reaction_4, reaction_5, reaction_6,
        reaction_7, reaction_8, reaction_9, reaction_10, reaction_11, reaction_12,
        reaction_13, reaction_14, reaction_15, reaction_16, reaction_17, reaction_18,
        reaction_19, reaction_20, reaction_21, reaction_22, reaction_23, reaction_24,
        reaction_25, reaction_26, reaction_27, reaction_28, reaction_29, reaction_30,
        reaction_31, reaction_32, reaction_33, reaction_34
    ])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=400)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_hatakeyama1()

start_time = time.time()
results = model.run(algorithm="ODE")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

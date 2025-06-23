import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create bachmann model
def create_bachmann(parameter_values=None):
    model = gillespy2.Model(name="bachmann")

    # Define Variables (GillesPy2.Species)
    CIS = gillespy2.Species(name="CIS", initial_value=0.0, mode="continuous")
    CISRNA = gillespy2.Species(name="CISRNA", initial_value=0.0, mode="continuous")
    Epo = gillespy2.Species(name="Epo", initial_value=0.000000124997, mode="continuous")
    EpoRJAK2 = gillespy2.Species(name="EpoRJAK2", initial_value=3.97622, mode="continuous")
    EpoRJAK2_CIS = gillespy2.Species(name="EpoRJAK2_CIS", initial_value=0.0, mode="continuous")
    EpoRpJAK2 = gillespy2.Species(name="EpoRpJAK2", initial_value=0.0, mode="continuous")
    SHP1 = gillespy2.Species(name="SHP1", initial_value=26.7251, mode="continuous")
    SHP1Act = gillespy2.Species(name="SHP1Act", initial_value=0.0, mode="continuous")
    SOCS3 = gillespy2.Species(name="SOCS3", initial_value=0.0, mode="continuous")
    SOCS3RNA = gillespy2.Species(name="SOCS3RNA", initial_value=0.0, mode="continuous")
    STAT5 = gillespy2.Species(name="STAT5", initial_value=79.7535, mode="continuous")
    p12EpoRpJAK2 = gillespy2.Species(name="p12EpoRpJAK2", initial_value=0.0, mode="continuous")
    p1EpoRpJAK2 = gillespy2.Species(name="p1EpoRpJAK2", initial_value=0.0, mode="continuous")
    p2EpoRpJAK2 = gillespy2.Species(name="p2EpoRpJAK2", initial_value=0.0, mode="continuous")
    pSTAT5 = gillespy2.Species(name="pSTAT5", initial_value=0.0, mode="continuous")
    CISnRNA1 = gillespy2.Species(name="CISnRNA1", initial_value=0.0, mode="continuous")
    CISnRNA2 = gillespy2.Species(name="CISnRNA2", initial_value=0.0, mode="continuous")
    CISnRNA3 = gillespy2.Species(name="CISnRNA3", initial_value=0.0, mode="continuous")
    CISnRNA4 = gillespy2.Species(name="CISnRNA4", initial_value=0.0, mode="continuous")
    CISnRNA5 = gillespy2.Species(name="CISnRNA5", initial_value=0.0, mode="continuous")
    SOCS3nRNA1 = gillespy2.Species(name="SOCS3nRNA1", initial_value=0.0, mode="continuous")
    SOCS3nRNA2 = gillespy2.Species(name="SOCS3nRNA2", initial_value=0.0, mode="continuous")
    SOCS3nRNA3 = gillespy2.Species(name="SOCS3nRNA3", initial_value=0.0, mode="continuous")
    SOCS3nRNA4 = gillespy2.Species(name="SOCS3nRNA4", initial_value=0.0, mode="continuous")
    SOCS3nRNA5 = gillespy2.Species(name="SOCS3nRNA5", initial_value=0.0, mode="continuous")
    npSTAT5 = gillespy2.Species(name="npSTAT5", initial_value=0.0, mode="continuous")
    # Add Vairables to Model
    model.add_species([CIS, CISRNA, Epo, EpoRJAK2, EpoRJAK2_CIS, EpoRpJAK2, SHP1, SHP1Act, SOCS3, SOCS3RNA, STAT5, p12EpoRpJAK2, p1EpoRpJAK2, p2EpoRpJAK2, pSTAT5,
                       CISnRNA1, CISnRNA2, CISnRNA3, CISnRNA4, CISnRNA5, SOCS3nRNA1, SOCS3nRNA2, SOCS3nRNA3, SOCS3nRNA4, SOCS3nRNA5, npSTAT5])

    # Define Parameters
    ActD = gillespy2.Parameter(name="ActD", expression=0.0)
    CISEqc = gillespy2.Parameter(name="CISEqc", expression=432.871)
    CISEqcOE = gillespy2.Parameter(name="CISEqcOE", expression=0.530261)
    CISInh = gillespy2.Parameter(name="CISInh", expression=784653000.0)
    CISRNADelay = gillespy2.Parameter(name="CISRNADelay", expression=0.144775)
    CISRNAEqc = gillespy2.Parameter(name="CISRNAEqc", expression=1.0)
    CISRNATurn = gillespy2.Parameter(name="CISRNATurn", expression=1000.0)
    CISTurn = gillespy2.Parameter(name="CISTurn", expression=0.00839842)
    CISoe = gillespy2.Parameter(name="CISoe", expression=0.0)
    EpoRActJAK2 = gillespy2.Parameter(name="EpoRActJAK2", expression=0.267308)
    EpoRCISInh = gillespy2.Parameter(name="EpoRCISInh", expression=1000000.0)
    EpoRCISRemove = gillespy2.Parameter(name="EpoRCISRemove", expression=5.42932)
    JAK2ActEpo = gillespy2.Parameter(name="JAK2ActEpo", expression=633253.0)
    JAK2EpoRDeaSHP1 = gillespy2.Parameter(name="JAK2EpoRDeaSHP1", expression=142.722)
    SHP1ActEpoR = gillespy2.Parameter(name="SHP1ActEpoR", expression=0.001)
    SHP1Dea = gillespy2.Parameter(name="SHP1Dea", expression=0.00816391)
    SOCS3Eqc = gillespy2.Parameter(name="SOCS3Eqc", expression=173.653)
    SOCS3EqcOE = gillespy2.Parameter(name="SOCS3EqcOE", expression=0.679157)
    SOCS3Inh = gillespy2.Parameter(name="SOCS3Inh", expression=10.408)
    SOCS3RNADelay = gillespy2.Parameter(name="SOCS3RNADelay", expression=1.06465)
    SOCS3RNAEqc = gillespy2.Parameter(name="SOCS3RNAEqc", expression=1.0)
    SOCS3RNATurn = gillespy2.Parameter(name="SOCS3RNATurn", expression=0.00830844)
    SOCS3Turn = gillespy2.Parameter(name="SOCS3Turn", expression=10000.0)
    SOCS3oe = gillespy2.Parameter(name="SOCS3oe", expression=0.0)
    STAT5ActEpoR = gillespy2.Parameter(name="STAT5ActEpoR", expression=38.9757)
    STAT5ActJAK2 = gillespy2.Parameter(name="STAT5ActJAK2", expression=0.0780965)
    STAT5Exp = gillespy2.Parameter(name="STAT5Exp", expression=0.0745155)
    STAT5Imp = gillespy2.Parameter(name="STAT5Imp", expression=0.0268889)
    epo_level = gillespy2.Parameter(name="epo_level", expression=0.000000124997)
    init_EpoRJAK2 = gillespy2.Parameter(name="init_EpoRJAK2", expression=3.97622)
    init_SHP1 = gillespy2.Parameter(name="init_SHP1", expression=26.7251)
    init_STAT5 = gillespy2.Parameter(name="init_STAT5", expression=79.7535)
    cyt = gillespy2.Parameter(name="cyt", expression=0.4)
    nuc = gillespy2.Parameter(name="nuc", expression=0.275)
    
    # Add Parameters to Model
    model.add_parameter([
    ActD, CISEqc, CISEqcOE, CISInh, CISRNADelay, CISRNAEqc, CISRNATurn, CISTurn, CISoe, EpoRActJAK2, EpoRCISInh, EpoRCISRemove, JAK2ActEpo,
    JAK2EpoRDeaSHP1, SHP1ActEpoR, SHP1Dea, SOCS3Eqc, SOCS3EqcOE, SOCS3Inh, SOCS3RNADelay, SOCS3RNAEqc, SOCS3RNATurn, SOCS3Turn, SOCS3oe,
    STAT5ActEpoR, STAT5ActJAK2, STAT5Exp, STAT5Imp, epo_level, init_EpoRJAK2, init_SHP1, init_STAT5, cyt, nuc
    ])

    # Define Reactions
    reaction_1 = gillespy2.Reaction(
        name="reaction_1",
        reactants={'EpoRJAK2': 1},
        products={'EpoRpJAK2': 1},
        propensity_function="JAK2ActEpo * Epo * EpoRJAK2 / (SOCS3Inh * SOCS3 / SOCS3Eqc + 1) * cyt"
    )

    reaction_2 = gillespy2.Reaction(
        name="reaction_2",
        reactants={'EpoRpJAK2': 1},
        products={'EpoRJAK2': 1},
        propensity_function="JAK2EpoRDeaSHP1 * SHP1Act * EpoRpJAK2 / init_SHP1 * cyt"
    )

    reaction_3 = gillespy2.Reaction(
        name="reaction_3",
        reactants={'EpoRpJAK2': 1},
        products={'p1EpoRpJAK2': 1},
        propensity_function="EpoRActJAK2 * EpoRpJAK2 / (SOCS3Inh * SOCS3 / SOCS3Eqc + 1) * cyt"
    )

    reaction_4 = gillespy2.Reaction(
        name="reaction_4",
        reactants={'EpoRpJAK2': 1},
        products={'p2EpoRpJAK2': 1},
        propensity_function="3 * EpoRActJAK2 * EpoRpJAK2 / ((SOCS3Inh * SOCS3 / SOCS3Eqc + 1) * (EpoRCISInh * EpoRJAK2_CIS + 1)) * cyt"
    )

    reaction_5 = gillespy2.Reaction(
        name="reaction_5",
        reactants={'p1EpoRpJAK2': 1},
        products={'p12EpoRpJAK2': 1},
        propensity_function="3 * EpoRActJAK2 * p1EpoRpJAK2 / ((SOCS3Inh * SOCS3 / SOCS3Eqc + 1) * (EpoRCISInh * EpoRJAK2_CIS + 1)) * cyt"
    )

    reaction_6 = gillespy2.Reaction(
        name="reaction_6",
        reactants={'p2EpoRpJAK2': 1},
        products={'p12EpoRpJAK2': 1},
        propensity_function="EpoRActJAK2 * p2EpoRpJAK2 / (SOCS3Inh * SOCS3 / SOCS3Eqc + 1) * cyt"
    )

    reaction_7 = gillespy2.Reaction(
        name="reaction_7",
        reactants={'p1EpoRpJAK2': 1},
        products={'EpoRJAK2': 1},
        propensity_function="JAK2EpoRDeaSHP1 * SHP1Act * p1EpoRpJAK2 / init_SHP1 * cyt"
    )

    reaction_8 = gillespy2.Reaction(
        name="reaction_8",
        reactants={'p2EpoRpJAK2': 1},
        products={'EpoRJAK2': 1},
        propensity_function="JAK2EpoRDeaSHP1 * SHP1Act * p2EpoRpJAK2 / init_SHP1 * cyt"
    )

    reaction_9 = gillespy2.Reaction(
        name="reaction_9",
        reactants={'p12EpoRpJAK2': 1},
        products={'EpoRJAK2': 1},
        propensity_function="JAK2EpoRDeaSHP1 * SHP1Act * p12EpoRpJAK2 / init_SHP1 * cyt"
    )

    reaction_10 = gillespy2.Reaction(
        name="reaction_10",
        reactants={'EpoRJAK2_CIS': 1},
        products={},
        propensity_function="EpoRCISRemove * EpoRJAK2_CIS * (p12EpoRpJAK2 + p1EpoRpJAK2) / init_EpoRJAK2 * cyt"
    )

    reaction_11 = gillespy2.Reaction(
        name="reaction_11",
        reactants={'SHP1': 1},
        products={'SHP1Act': 1},
        propensity_function="SHP1ActEpoR * SHP1 * (EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2) / init_EpoRJAK2 * cyt"
    )

    reaction_12 = gillespy2.Reaction(
        name="reaction_12",
        reactants={'SHP1Act': 1},
        products={'SHP1': 1},
        propensity_function="SHP1Dea * SHP1Act * cyt"
    )

    reaction_13 = gillespy2.Reaction(
        name="reaction_13",
        reactants={'STAT5': 1},
        products={'pSTAT5': 1},
        propensity_function="STAT5ActJAK2 * STAT5 * (EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2) / (init_EpoRJAK2 * (SOCS3Inh * SOCS3 / SOCS3Eqc + 1)) * cyt"
    )

    reaction_14 = gillespy2.Reaction(
        name="reaction_14",
        reactants={'STAT5': 1},
        products={'pSTAT5': 1},
        propensity_function="STAT5ActEpoR * STAT5 * pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2) / (pow(init_EpoRJAK2, 2) * (CISInh * CIS / CISEqc + 1) * (SOCS3Inh * SOCS3 / SOCS3Eqc + 1)) * cyt"
    )

    reaction_15 = gillespy2.Reaction(
        name="reaction_15",
        reactants={'pSTAT5': 1},
        products={'npSTAT5': 1},
        propensity_function="STAT5Imp * pSTAT5 * cyt"
    )

    reaction_16 = gillespy2.Reaction(
        name="reaction_16",
        reactants={'npSTAT5': 1},
        products={'STAT5': 1},
        propensity_function="STAT5Exp * npSTAT5 * nuc"
    )   

    reaction_17 = gillespy2.Reaction(
        name="reaction_17",
        reactants={},
        products={'CISnRNA1': 1},
        propensity_function="CISRNAEqc * CISRNATurn * npSTAT5 * (1 - ActD) / init_STAT5 * nuc"
    )

    reaction_18 = gillespy2.Reaction(
        name="reaction_18",
        reactants={'CISnRNA1': 1},
        products={'CISnRNA2': 1},
        propensity_function="CISRNADelay * CISnRNA1 * nuc"
    )

    reaction_19 = gillespy2.Reaction(
        name="reaction_19",
        reactants={'CISnRNA2': 1},
        products={'CISnRNA3': 1},
        propensity_function="CISRNADelay * CISnRNA2 * nuc"
    )

    reaction_20 = gillespy2.Reaction(
        name="reaction_20",
        reactants={'CISnRNA3': 1},
        products={'CISnRNA4': 1},
        propensity_function="CISRNADelay * CISnRNA3 * nuc"
    )

    reaction_21 = gillespy2.Reaction(
        name="reaction_21",
        reactants={'CISnRNA4': 1},
        products={'CISnRNA5': 1},
        propensity_function="CISRNADelay * CISnRNA4 * nuc"
    )

    reaction_22 = gillespy2.Reaction(
        name="reaction_22",
        reactants={'CISnRNA5': 1},
        products={'CISRNA': 1},
        propensity_function="CISRNADelay * CISnRNA5 * nuc"
    )

    # Cytoplasmic RNA and protein reactions
    reaction_23 = gillespy2.Reaction(
        name="reaction_23",
        reactants={'CISRNA': 1},
        products={},
        propensity_function="CISRNATurn * CISRNA * cyt"
    )

    reaction_24 = gillespy2.Reaction(
        name="reaction_24",
        reactants={},
        products={'CIS': 1},
        propensity_function="CISEqc * CISTurn * CISRNA / CISRNAEqc * cyt"
    )

    reaction_25 = gillespy2.Reaction(
        name="reaction_25",
        reactants={'CIS': 1},
        products={},
        propensity_function="CISTurn * CIS * cyt"
    )

    reaction_26 = gillespy2.Reaction(
        name="reaction_26",
        reactants={},
        products={'CIS': 1},
        propensity_function="CISoe * CISEqc * CISTurn * CISEqcOE * cyt"
    )

    reaction_27 = gillespy2.Reaction(
        name="reaction_27",
        reactants={},
        products={'SOCS3nRNA1': 1},
        propensity_function="SOCS3RNAEqc * SOCS3RNATurn * npSTAT5 * (1 - ActD) / init_STAT5 * nuc"
    )

    reaction_28 = gillespy2.Reaction(
        name="reaction_28",
        reactants={'SOCS3nRNA1': 1},
        products={'SOCS3nRNA2': 1},
        propensity_function="SOCS3RNADelay * SOCS3nRNA1 * nuc"
    )

    reaction_29 = gillespy2.Reaction(
        name="reaction_29",
        reactants={'SOCS3nRNA2': 1},
        products={'SOCS3nRNA3': 1},
        propensity_function="SOCS3RNADelay * SOCS3nRNA2 * nuc"
    )

    reaction_30 = gillespy2.Reaction(
        name="reaction_30",
        reactants={'SOCS3nRNA3': 1},
        products={'SOCS3nRNA4': 1},
        propensity_function="SOCS3RNADelay * SOCS3nRNA3 * nuc"
    )

    reaction_31 = gillespy2.Reaction(
        name="reaction_31",
        reactants={'SOCS3nRNA4': 1},
        products={'SOCS3nRNA5': 1},
        propensity_function="SOCS3RNADelay * SOCS3nRNA4 * nuc"
    )

    reaction_32 = gillespy2.Reaction(
        name="reaction_32",
        reactants={'SOCS3nRNA5': 1},
        products={'SOCS3RNA': 1},
        propensity_function="SOCS3RNADelay * SOCS3nRNA5 * nuc"
    )

    reaction_33 = gillespy2.Reaction(
        name="reaction_33",
        reactants={'SOCS3RNA': 1},
        products={},
        propensity_function="SOCS3RNATurn * SOCS3RNA * cyt"
    )

    reaction_34 = gillespy2.Reaction(
        name="reaction_34",
        reactants={},
        products={'SOCS3': 1},
        propensity_function="SOCS3Eqc * SOCS3Turn * SOCS3RNA / SOCS3RNAEqc * cyt"
    )

    reaction_35 = gillespy2.Reaction(
        name="reaction_35",
        reactants={'SOCS3': 1},
        products={},
        propensity_function="SOCS3Turn * SOCS3 * cyt"
    )

    reaction_36 = gillespy2.Reaction(
        name="reaction_36",
        reactants={},
        products={'SOCS3': 1},
        propensity_function="SOCS3oe * SOCS3Eqc * SOCS3Turn * SOCS3EqcOE * cyt"
    )

    # Add reactions to model
    model.add_reaction([
        reaction_1, reaction_2, reaction_3, reaction_4, reaction_5, reaction_6,
        reaction_7, reaction_8, reaction_9, reaction_10, reaction_11, reaction_12,
        reaction_13, reaction_14, reaction_15, reaction_16,
        reaction_17, reaction_18, reaction_19, reaction_20, reaction_21, reaction_22,
        reaction_23, reaction_24, reaction_25, reaction_26,
        reaction_27, reaction_28, reaction_29, reaction_30, reaction_31, reaction_32,
        reaction_33, reaction_34, reaction_35, reaction_36
    ])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=400, num_points=401)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_bachmann()

start_time = time.time()
results = model.run(algorithm="ODE")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

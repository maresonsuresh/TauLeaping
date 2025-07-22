import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create heatshock model
def create_heatshock(parameter_values=None):
    model = gillespy2.Model(name="heatshock")
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    NatP = gillespy2.Species(name="NatP", initial_value=6000000, mode="discrete")  
    HCom = gillespy2.Species(name="HCom", initial_value=5900, mode="discrete")  
    Hsp90 = gillespy2.Species(name="Hsp90", initial_value=300000, mode="discrete")  
    HSF1 = gillespy2.Species(name="HSF1", initial_value=100, mode="discrete")  
    ROS = gillespy2.Species(name="ROS", initial_value=100, mode="discrete")  
    ATP = gillespy2.Species(name="ATP", initial_value=10000, mode="discrete")  
    ADP = gillespy2.Species(name="ADP", initial_value=1000, mode="discrete")  
    HSE = gillespy2.Species(name="HSE", initial_value=1, mode="discrete")  
    DiH = gillespy2.Species(name="DiH", initial_value=0, mode="discrete")  
    TriH = gillespy2.Species(name="TriH", initial_value=0, mode="discrete")  
    HSETriH = gillespy2.Species(name="HSETriH", initial_value=0, mode="discrete")  
    MisP = gillespy2.Species(name="MisP", initial_value=0, mode="discrete")
    MCom = gillespy2.Species(name="MCom", initial_value=0, mode="discrete")
    AggP = gillespy2.Species(name="AggP", initial_value=0, mode="discrete")
    # Add Variables to Model
    model.add_species([NatP, HCom, Hsp90, HSF1, ROS, ATP, ADP, HSE, DiH, TriH, HSETriH, MisP, MCom, AggP])

    # Protein folding and heat shock response parameters
    k1 = gillespy2.Parameter(name="k1", expression=10.0)  # Protein synthesis (mol/s)
    k2 = gillespy2.Parameter(name="k2", expression=0.00002)  # Misfolding (mol⁻¹s⁻¹)
    k3 = gillespy2.Parameter(name="k3", expression=50.0)  # Hsp90 binds misfolded protein (mol⁻¹s⁻¹)
    k4 = gillespy2.Parameter(name="k4", expression=0.00001)  # Misfolded complex dissociation (s⁻¹)
    k5 = gillespy2.Parameter(name="k5", expression=4.0e6)  # Refolding (mol⁻¹s⁻¹)
    k6 = gillespy2.Parameter(name="k6", expression=6.0e7)  # Protein degradation (mol⁻¹s⁻¹)
    k7 = gillespy2.Parameter(name="k7", expression=1.0e7)  # Aggregation (mol⁻¹s⁻¹)
    k8 = gillespy2.Parameter(name="k8", expression=500.0)  # Hsp90 binds HSF1 (mol⁻¹s⁻¹)
    k9 = gillespy2.Parameter(name="k9", expression=1.0)  # HSF1 complex dissociation (s⁻¹)
    k10 = gillespy2.Parameter(name="k10", expression=0.01)  # HSF1 dimerization (mol⁻¹s⁻¹)
    k11 = gillespy2.Parameter(name="k11", expression=100.0)  # HSF1 trimerization (mol⁻¹s⁻¹)
    k12 = gillespy2.Parameter(name="k12", expression=0.5)  # Trimer dissociation (s⁻¹)
    k13 = gillespy2.Parameter(name="k13", expression=0.5)  # Dimer dissociation (s⁻¹)
    k14 = gillespy2.Parameter(name="k14", expression=0.05)  # HSE binding (mol⁻¹s⁻¹)
    k15 = gillespy2.Parameter(name="k15", expression=0.08)  # HSE dissociation (s⁻¹)
    k16 = gillespy2.Parameter(name="k16", expression=1000.0)  # Hsp90 transcription (s⁻¹)
    k17 = gillespy2.Parameter(name="k17", expression=8.02e9)  # Hsp90 degradation (s⁻¹)
    k18 = gillespy2.Parameter(name="k18", expression=12.0)  # ATP formation (s⁻¹)
    k19 = gillespy2.Parameter(name="k19", expression=0.02)  # ADP formation (s⁻¹)
    k20 = gillespy2.Parameter(name="k20", expression=0.1)  # ROS production (mol/s)
    k21 = gillespy2.Parameter(name="k21", expression=0.001)  # ROS removal (s⁻¹)
    # Add Parameters to Model
    model.add_parameter([
        k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
        k11, k12, k13, k14, k15, k16, k17, k18, k19,
        k20, k21
    ])

    # Define Reactions
    r1 = gillespy2.Reaction(name="r1", reactants={}, products={'NatP': 1}, rate="k1")
    # Protein misfolding (NatP + ROS → MisP + ROS)
    r2 = gillespy2.Reaction(
        name="r2", 
        reactants={'NatP':1, 'ROS':1}, 
        products={'MisP':1, 'ROS':1}, 
        rate="k2"
    )

    # Hsp90 binding to misfolded protein (MisP + Hsp90 → MCom)
    r3 = gillespy2.Reaction(
        name="r3",
        reactants={'MisP':1, 'Hsp90':1},
        products={'MCom':1},
        rate="k3"
    )

    # Misfolded complex dissociation (MCom → MisP + Hsp90)
    r4 = gillespy2.Reaction(
        name="r4",
        reactants={'MCom':1},
        products={'MisP':1, 'Hsp90':1},
        rate="k4"
    )

    # Protein refolding (MCom + ATP → NatP + Hsp90 + ADP)
    r5 = gillespy2.Reaction(
        name="r5",
        reactants={'MCom':1, 'ATP':1},
        products={'NatP':1, 'Hsp90':1, 'ADP':1},
        rate="k5"
    )

    # Protein degradation (MisP + ATP → ADP)
    r6 = gillespy2.Reaction(
        name="r6",
        reactants={'MisP':1, 'ATP':1},
        products={'ADP':1},
        rate="k6"
    )

    # Protein aggregation (2MisP → AggP)
    r7 = gillespy2.Reaction(
        name="r7",
        reactants={'MisP':2},
        products={'AggP':1},
        rate="k7"
    )

    # Hsp90 binding to HSF1 (Hsp90 + HSF1 → HCom)
    r8 = gillespy2.Reaction(
        name="r8",
        reactants={'Hsp90':1, 'HSF1':1},
        products={'HCom':1},
        rate="k8"
    )

    # HSF1 complex dissociation (HCom → Hsp90 + HSF1)
    r9 = gillespy2.Reaction(
        name="r9",
        reactants={'HCom':1},
        products={'Hsp90':1, 'HSF1':1},
        rate="k9"
    )

    # HSF1 dimerization (2HSF1 → DiH)
    r10 = gillespy2.Reaction(
        name="r10",
        reactants={'HSF1':2},
        products={'DiH':1},
        rate="k10"
    )

    # HSF1 trimerization (HSF1 + DiH → TriH)
    r11 = gillespy2.Reaction(
        name="r11",
        reactants={'HSF1':1, 'DiH':1},
        products={'TriH':1},
        rate="k11"
    )

    # HSF1 trimer dissociation (TriH → HSF1 + DiH)
    r12 = gillespy2.Reaction(
        name="r12",
        reactants={'TriH':1},
        products={'HSF1':1, 'DiH':1},
        rate="k12"
    )

    # HSF1 dimer dissociation (DiH → 2HSF1)
    r13 = gillespy2.Reaction(
        name="r13",
        reactants={'DiH':1},
        products={'HSF1':2},
        rate="k13"
    )

    # HSE binding (TriH + HSE → HSETriH)
    r14 = gillespy2.Reaction(
        name="r14",
        reactants={'TriH':1, 'HSE':1},
        products={'HSETriH':1},
        rate="k14"
    )

    # HSE dissociation (HSETriH → TriH + HSE)
    r15 = gillespy2.Reaction(
        name="r15",
        reactants={'HSETriH':1},
        products={'TriH':1, 'HSE':1},
        rate="k15"
    )

    # Hsp90 transcription (HSETriH → HSETriH + Hsp90)
    r16 = gillespy2.Reaction(
        name="r16",
        reactants={'HSETriH':1},
        products={'HSETriH':1, 'Hsp90':1},
        rate="k16"
    )

    # Hsp90 degradation (Hsp90 + ATP → ADP)
    r17 = gillespy2.Reaction(
        name="r17",
        reactants={'Hsp90':1, 'ATP':1},
        products={'ADP':1},
        rate="k17"
    )

    # ATP generation (ADP → ATP)
    r18 = gillespy2.Reaction(
        name="r18",
        reactants={'ADP':1},
        products={'ATP':1},
        rate="k18"
    )

    # ATP consumption (ATP → ADP)
    r19 = gillespy2.Reaction(
        name="r19",
        reactants={'ATP':1},
        products={'ADP':1},
        rate="k19"
    )

    # ROS production (∅ → ROS)
    r20 = gillespy2.Reaction(
        name="r20",
        reactants={},
        products={'ROS':1},
        rate="k20"
    )

    # ROS removal (ROS → ∅)
    r21 = gillespy2.Reaction(
        name="r21",
        reactants={'ROS':1},
        products={},
        rate="k21"
    )


    # Add Reactions to Model
    model.add_reaction([
        r2, r3, r4, r5, r6, r7, r8, r9, r10,
        r11, r12, r13, r14, r15, r16, r17,
        r18, r19, r20, r21
    ])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=1, num_points=2)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_heatshock()

start_time = time.time()
results = model.run(algorithm="Tau-Leaping")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

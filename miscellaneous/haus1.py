import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
from math import tanh
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create haus1 model
def create_haus1(parameter_values=None):
    model = gillespy2.Model(name="haus1")
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    A = gillespy2.Species(name="A", initial_value=43, mode="discrete")
    AC = gillespy2.Species(name="AC", initial_value=0, mode="discrete")
    Aa = gillespy2.Species(name="Aa", initial_value=0, mode="discrete")
    AaC = gillespy2.Species(name="AaC", initial_value=0, mode="discrete")
    Ad = gillespy2.Species(name="Ad", initial_value=0, mode="discrete")
    Ah = gillespy2.Species(name="Ah", initial_value=0, mode="discrete")
    An = gillespy2.Species(name="An", initial_value=3, mode="discrete")
    B = gillespy2.Species(name="B", initial_value=56, mode="discrete")
    BC = gillespy2.Species(name="BC", initial_value=0, mode="discrete")
    Bn = gillespy2.Species(name="Bn", initial_value=0, mode="discrete")
    Cf = gillespy2.Species(name="Cf", initial_value=0, mode="discrete")
    En = gillespy2.Species(name="En", initial_value=6, mode="discrete")

    # Add Vairables to Model
    model.add_species([A, AC, Aa, AaC, Ad, Ah, An, B, BC, Bn, Cf, En])

    # Define Parameters
    Chemostat = gillespy2.Parameter(name="Chemostat", expression=1.0)
    Dvar = gillespy2.Parameter(name="Dvar", expression=0.075)
    G = gillespy2.Parameter(name="G", expression=40.0)
    K1 = gillespy2.Parameter(name="K1", expression=0.00158)
    K10 = gillespy2.Parameter(name="K10", expression=0.000014)
    K2 = gillespy2.Parameter(name="K2", expression=0.00181)
    K4 = gillespy2.Parameter(name="K4", expression=1.87)
    K8 = gillespy2.Parameter(name="K8", expression=0.00000792)
    V1 = gillespy2.Parameter(name="V1", expression=4.94)
    V10 = gillespy2.Parameter(name="V10", expression=4.75)
    V2 = gillespy2.Parameter(name="V2", expression=2.92)
    V4 = gillespy2.Parameter(name="V4", expression=22.8)
    V8 = gillespy2.Parameter(name="V8", expression=64.8)
    alpha3 = gillespy2.Parameter(name="alpha3", expression=0.00517)
    alpha5 = gillespy2.Parameter(name="alpha5", expression=0.014)
    alpha6 = gillespy2.Parameter(name="alpha6", expression=0.00537)
    alpha7 = gillespy2.Parameter(name="alpha7", expression=4790.0)
    alpha9 = gillespy2.Parameter(name="alpha9", expression=347000.0)
    compartment = gillespy2.Parameter(name="compartment", expression=1.0)
    k0 = gillespy2.Parameter(name="k0", expression=0.0)
    n = gillespy2.Parameter(name="n", expression=485.0)
    pstar = gillespy2.Parameter(name="pstar", expression=4.5)
    rAd = gillespy2.Parameter(name="rAd", expression=0.00547)
    rAdplus = gillespy2.Parameter(name="rAdplus", expression=0.1037)
    rAh = gillespy2.Parameter(name="rAh", expression=0.289)
    rAhplus = gillespy2.Parameter(name="rAhplus", expression=2.5594)
    rCf = gillespy2.Parameter(name="rCf", expression=0.000324)
    rCfplus = gillespy2.Parameter(name="rCfplus", expression=1.063)

    # Add Parameters to Model
    model.add_parameter([
        Chemostat, Dvar, G, K1, K10, K2, K4, K8,
        V1, V10, V2, V4, V8,
        alpha3, alpha5, alpha6, alpha7, alpha9,
        compartment, k0, n, pstar,
        rAd, rAdplus, rAh, rAhplus, rCf, rCfplus
    ])

    # Define Reactions
    v_1 = gillespy2.Reaction(
        name="v_1",
        reactants={},
        products={'AC': 1},
        propensity_function="(2*V1*G/(K1+G))"
    )

    v_2 = gillespy2.Reaction(
        name="v_2",
        reactants={'AC': 1},
        products={'A': 1},
        propensity_function="V2*AC/(K2+AC)"
    )

    v_3 = gillespy2.Reaction(
        name="v_3",
        reactants={'A': 1, 'AaC': 1},
        products={'Aa': 1, 'AC': 1},
        propensity_function="alpha3 * A * AaC * Cf"
    )

    v_4 = gillespy2.Reaction(
        name="v_4",
        reactants={'AC': 1},
        products={'AaC': 1},
        propensity_function="V4*AC/(K4+AC)"
    )

    v_5 = gillespy2.Reaction(
        name="v_5",
        reactants={'AC': 1},
        products={'En': 1},
        propensity_function="alpha5 * AC * Ah"
    )

    v_6 = gillespy2.Reaction(
        name="v_6",
        reactants={'AaC': 1, 'B': 1},
        products={'Aa': 1, 'BC': 1},
        propensity_function="alpha6 * B * AaC * Cf"
    )

    v_7 = gillespy2.Reaction(
        name="v_7",
        reactants={'Aa': 1},
        products={'An': 1},
        propensity_function="alpha7 * Aa * Ad"
    )

    v_8 = gillespy2.Reaction(
        name="v_8",
        reactants={'BC': 1},
        products={'B': 1},
        propensity_function="V8*BC/(K8+BC)"
    )

    v_9 = gillespy2.Reaction(
        name="v_9",
        reactants={'BC': 1},
        products={'Bn': 1},
        propensity_function="alpha9 * BC * Ah"
    )

    v_10 = gillespy2.Reaction(
        name="v_10",
        reactants={},
        products={'Ad': 1},
        propensity_function="(rAd + rAdplus * (1 - tanh(n * (5.7 - 0.6 * tanh(141.9436 * 0.1846) + 0.6 * tanh(-0.1846 * (time - 141.9436)) - pstar)))"
    )

    v_11 = gillespy2.Reaction(
        name="v_11",
        reactants={},
        products={'Cf': 1},
        propensity_function="(rCf + rCfplus * (1 - tanh(n * (5.7 - 0.6 * tanh(141.9436 * 0.1846) + 0.6 * tanh(-0.1846 * (time - 141.9436)) - pstar)))"
    )

    v_12 = gillespy2.Reaction(
        name="v_12",
        reactants={},
        products={'Ah': 1},
        propensity_function="(rAh + rAhplus * (1 - tanh(n * (5.7 - 0.6 * tanh(141.9436 * 0.1846) + 0.6 * tanh(-0.1846 * (time - 141.9436)) - pstar)))"
    )

    v_13 = gillespy2.Reaction(
        name="v_13",
        reactants={'AaC': 1},
        products={'BC': 1},
        propensity_function="V10*AaC/(K10+AaC)"
    )

    v_14 = gillespy2.Reaction(
        name="v_14",
        reactants={'AC': 1},
        products={},
        propensity_function="Dvar*AC"
    )

    v_15 = gillespy2.Reaction(
        name="v_15",
        reactants={'A': 1},
        products={},
        propensity_function="Dvar*A"
    )

    v_16 = gillespy2.Reaction(
        name="v_16",
        reactants={'En': 1},
        products={},
        propensity_function="Dvar*En"
    )

    v_17 = gillespy2.Reaction(
        name="v_17",
        reactants={'AaC': 1},
        products={},
        propensity_function="Dvar*AaC"
    )

    v_18 = gillespy2.Reaction(
        name="v_18",
        reactants={'Aa': 1},
        products={},
        propensity_function="Dvar*Aa"
    )

    v_19 = gillespy2.Reaction(
        name="v_19",
        reactants={'BC': 1},
        products={},
        propensity_function="Dvar*BC"
    )

    v_20 = gillespy2.Reaction(
        name="v_20",
        reactants={'B': 1},
        products={},
        propensity_function="Dvar*B"
    )

    v_21 = gillespy2.Reaction(
        name="v_21",
        reactants={'An': 1},
        products={},
        propensity_function="Dvar*An"
    )

    v_22 = gillespy2.Reaction(
        name="v_22",
        reactants={'Bn': 1},
        products={},
        propensity_function="Dvar*Bn"
    )

    # Add all reactions to the model
    model.add_reaction([
        v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9,
        v_10, v_11, v_12, v_13, v_14, v_15, v_16,
        v_17, v_18, v_19, v_20, v_21, v_22
    ])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=400)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_haus1()

start_time = time.time()
results = model.run(algorithm="ODE")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

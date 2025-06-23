import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create zhao1 model
def create_zhao1(parameter_values=None):
    model = gillespy2.Model(name="zhao1")
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    I = gillespy2.Species(name="I", initial_value=149400, mode="discrete")
    Ip = gillespy2.Species(name="Ip", initial_value=16600, mode="discrete")
    S = gillespy2.Species(name="S", initial_value=667200, mode="discrete")
    Sp = gillespy2.Species(name="Sp", initial_value=166800.0, mode="discrete")
    # Add Vairables to Model
    model.add_species([I, Ip, S, Sp])

    # Define Parameters
    Lambda = gillespy2.Parameter(name="Lambda", expression=38094.0)
    alphai = gillespy2.Parameter(name="alphai", expression=0.5)
    alphas = gillespy2.Parameter(name="alphas", expression=0.5)
    ba = gillespy2.Parameter(name="ba", expression=0.003)
    d = gillespy2.Parameter(name="d", expression=0.1302)
    k = gillespy2.Parameter(name="k", expression=0.2)
    mu = gillespy2.Parameter(name="mu", expression=0.025)
    n = gillespy2.Parameter(name="n", expression=65.8494)
    
    # Add Parameters to Model
    model.add_parameter([
        Lambda, alphai, alphas, ba, d, k, mu, n
    ])

    # Define Reactions
    v1 = gillespy2.Reaction(name="v1", reactants={}, products={'Sp': 1}, propensity_function="k*Lambda")
    v2 = gillespy2.Reaction(name="v2", reactants={'I': 1}, products={}, propensity_function="d*I")
    v3 = gillespy2.Reaction(name="v3", reactants={'Sp': 1}, products={'Ip': 1}, propensity_function="((1-alphas)*(1-(1-ba)**n)*I/(Sp + S + Ip + I)+(1-alphas)*(1-alphai)*(1-(1-ba)**n)*Ip/(Sp + S + Ip + I))*Sp")
    v4 = gillespy2.Reaction(name="v4", reactants={'Sp': 1}, products={}, propensity_function="mu*Sp")
    v5 = gillespy2.Reaction(name="v5", reactants={}, products={'S': 1}, propensity_function="(1-k)*Lambda")
    v6 = gillespy2.Reaction(name="v6", reactants={'S': 1}, products={'I': 1}, propensity_function="((1-(1-ba)**n)*I/(Sp + S + Ip + I)+(1-alphai)*(1-(1-ba)**n)*Ip/(Sp + S + Ip + I))*S")
    v7 = gillespy2.Reaction(name="v7", reactants={'S': 1}, products={}, propensity_function="mu*S")
    v8 = gillespy2.Reaction(name="v8", reactants={'Ip': 1}, products={}, propensity_function="mu*Ip")
    v9 = gillespy2.Reaction(name="v9", reactants={'Ip': 1}, products={}, propensity_function="d*Ip")
    v10 = gillespy2.Reaction(name="v10", reactants={'I': 1}, products={}, propensity_function="mu*I")

    # Add Reactions to Model
    model.add_reaction([v1, v2, v3, v4, v5, v6, v7, v8, v9, v10])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=1000)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_zhao1()

start_time = time.time()
results = model.run(algorithm="Tau-Hybrid")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

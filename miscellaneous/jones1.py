import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))

# Create jones1 model
def create_jones1(parameter_values=None):
    model = gillespy2.Model(name="jones1")
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    A = gillespy2.Species(name="A", initial_value=1.0, mode="discrete")
    CC = gillespy2.Species(name="CC", initial_value=18800.0, mode="discrete")
    T = gillespy2.Species(name="T", initial_value=323000.0, mode="discrete")
    Tstar = gillespy2.Species(name="Tstar", initial_value=7800.0, mode="discrete")
    V = gillespy2.Species(name="V", initial_value=120000.0, mode="discrete")
    # Add Vairables to Model
    model.add_species([A, CC, T, Tstar, V])

    # Define Parameters
    KK = gillespy2.Parameter(name="KK", expression=3.36)
    NC = gillespy2.Parameter(name="NC", expression=4.11)
    NN = gillespy2.Parameter(name="NN", expression=285.0)
    a = gillespy2.Parameter(name="a", expression=1.6)
    alpha = gillespy2.Parameter(name="alpha", expression=0.195)
    c = gillespy2.Parameter(name="c", expression=13.0)
    d = gillespy2.Parameter(name="d", expression=0.01)
    delta = gillespy2.Parameter(name="delta", expression=0.7)
    gamma = gillespy2.Parameter(name="gamma", expression=0.0000017475)
    k = gillespy2.Parameter(name="k", expression=0.00000017437)
    Lambda = gillespy2.Parameter(name="Lambda", expression=10000.0)
    mu = gillespy2.Parameter(name="mu", expression=0.07)
    
    # Add Parameters to Model
    model.add_parameter([
        KK, NC, NN, a, alpha, c, d, delta, gamma, k, Lambda, mu
    ])

    # Define Reactions
    v1 = gillespy2.Reaction(name="v1", reactants={'A':1}, products={}, propensity_function="gamma*A*T")
    v2 = gillespy2.Reaction(name="v2", reactants={}, products={'T':1}, propensity_function="Lambda")
    v3 = gillespy2.Reaction(name="v3", reactants={}, products={'T':1}, propensity_function="(a*A*T)/(KK + A)")
    v4 = gillespy2.Reaction(name="v4", reactants={'T':1}, products={'Tstar':1}, propensity_function="k*T*V")
    v5 = gillespy2.Reaction(name="v5", reactants={'T':1}, products={}, propensity_function="d*T")
    v6 = gillespy2.Reaction(name="v6", reactants={'Tstar':1}, products={}, propensity_function="delta*Tstar")
    v7 = gillespy2.Reaction(name="v7", reactants={'CC':1}, products={}, propensity_function="mu*CC")
    v8 = gillespy2.Reaction(name="v8", reactants={'V':1}, products={}, propensity_function="c*V")
    v9 = gillespy2.Reaction(name="v9", reactants={'Tstar':1}, products={'CC':1}, propensity_function="alpha*k*T*V")
    v10 = gillespy2.Reaction(name="v10", reactants={}, products={'V':1}, propensity_function="delta*NN*Tstar")
    v11 = gillespy2.Reaction(name="v11", reactants={}, products={'V':1}, propensity_function="mu*NC*CC")

    # Add Reactions to Model
    model.add_reaction([v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=250)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_jones1()

start_time = time.time()
results = model.run(algorithm="ODE")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()

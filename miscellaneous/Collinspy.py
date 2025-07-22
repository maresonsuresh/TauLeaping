import gillespy2
from matplotlib import pyplot as plt
import time
import os
import sys

sys.path.insert(1, os.path.abspath(os.path.join(os.getcwd(), '../../..')))


# Create the Vilar Oscillator Model
# .............................................................................
def create_collins_toggle(parameter_values=None):
    model = gillespy2.Model(name="CollinsToggle")
    
    # Set System Volume
    model.volume = 1

    # Define Variables (GillesPy2.Species)
    x1 = gillespy2.Species(name="x1", initial_value=76, mode="discrete")
    x2 = gillespy2.Species(name="x2", initial_value=75, mode="discrete")
    x3 = gillespy2.Species(name="x3", initial_value=60, mode="discrete")
    x4 = gillespy2.Species(name="x4", initial_value=60, mode="discrete")

    # Add Vairables to Model
    model.add_species([x1, x2, x3, x4])

    # Define Parameters
    alpha = gillespy2.Parameter(name="alpha", expression=28.98)
    beta = gillespy2.Parameter(name="beta", expression=4)
    c = gillespy2.Parameter(name="c", expression=0.23)
    k = gillespy2.Parameter(name="k", expression=1000)

    # Add Parameters to Model
    model.add_parameter([
        alpha, beta, c, k
    ])

    # Define Reactions
    r1 = gillespy2.Reaction(name="r1", reactants={}, products={'C': 1}, propensity_function="alpha1/(1+pow(V,beta))")
    r2 = gillespy2.Reaction(name="r2", reactants={'A': 1}, products={}, propensity_function="alpha1/(1+pow(V,beta))")
    r3 = gillespy2.Reaction(name="r3", reactants={'C': 1}, products={'R': 1}, propensity_function="alpha1/(1+pow(V,beta))")
    r4 = gillespy2.Reaction(name="r4", reactants={'R': 1}, products={}, propensity_function="alpha1/(1+pow(V,beta))")
    r5 = gillespy2.Reaction(name="r5", reactants={'A': 1, 'Da': 1}, products={'Da_prime': 1}, propensity_function="alpha1/(1+pow(V,beta))")
    r6 = gillespy2.Reaction(name="r6", reactants={'Da_prime': 1}, products={'A': 1, 'Da': 1}, propensity_function="alpha1/(1+pow(V,beta))")
    r7 = gillespy2.Reaction(name="r7", reactants={'Da': 1}, products={'Da': 1, 'Ma': 1}, propensity_function="alpha1/(1+pow(V,beta))")
    r8 = gillespy2.Reaction(name="r8", reactants={'Da_prime': 1}, products={'Da_prime': 1, 'Ma': 1}, propensity_function="alpha1/(1+pow(V,beta))")

    # Add Reactions to Model
    model.add_reaction([r1, r2, r3, r4, r5, r6, r7, r8])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=8000)
    
    # Set Model Timespan
    model.timespan(tspan)
    return model

model = create_collins_toggle()

start_time = time.time()
results = model.run(algorithm="Tau-Hybrid")
end_time = time.time()
simulation_time = end_time - start_time
print(f"Simulation completed in {simulation_time:.4f} seconds.")

# Default plot
results.plot()
plt.show()


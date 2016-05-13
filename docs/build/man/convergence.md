
<a id='Convergence-Simulations-1'></a>

# Convergence Simulations


The convergence simulation


<a id='The-ConvergenceSimulation-Type-1'></a>

## The ConvergenceSimulation Type

<a id='DifferentialEquations.ConvergenceSimulation' href='#DifferentialEquations.ConvergenceSimulation'>#</a>
**`DifferentialEquations.ConvergenceSimulation`** &mdash; *Type*.



ConvergenceSimulation

A type which holds the data from a convergence simulation.


<a id='Related-Functions-1'></a>

## Related Functions

<a id='Base.length-Tuple{DifferentialEquations.ConvergenceSimulation}' href='#Base.length-Tuple{DifferentialEquations.ConvergenceSimulation}'>#</a>
**`Base.length`** &mdash; *Method*.



length(simres::ConvergenceSimulation)

Returns the number of simultations in the Convergence Simulation

<a id='DifferentialEquations.conv_ests' href='#DifferentialEquations.conv_ests'>#</a>
**`DifferentialEquations.conv_ests`** &mdash; *Function*.



conv_ests(error::Vector{Number})

Computes the pairwise convergence estimate for a convergence test done by halving/doubling stepsizes via

log2(error[i+1]/error[i])

Returns the mean of the convergence estimates


<a id='Plot-Functions-1'></a>

## Plot Functions


# Convergence Simulations

The convergence simulation type is useful for deriving order of convergence estimates
from a group of simulations. This object will automatically assemble error vectors
into a more useful manner and provide plotting functionality. Convergence estimates
are also given by pair-wise estimates.

`ConvergenceSimulation`s can either be created by passing the constructor the
appropriate solution array or by using one of the provided `test_convergence` functions.

## The ConvergenceSimulation Type

```@docs
DifferentialEquations.ConvergenceSimulation
```

## Plot Functions

The plot functionality is provided by a Plots.jl recipe. What is plotted is a
line series for each calculated error along the convergence axis. To plot a
convergence simulation, simply use:

```julia
plot(sim::ConvergenceSimulation)
```

All of the functionality (keyword arguments) provided by Plots.jl are able to
be used in this command. Please see the Plots.jl documentation for more information.

## Related Functions

```@docs
Base.length(::ConvergenceSimulation)
DifferentialEquations.conv_ests
DifferentialEquations.calc𝒪estimates
DifferentialEquations.test_convergence
```

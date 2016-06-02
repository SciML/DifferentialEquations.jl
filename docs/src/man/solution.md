# The Solution Type

Each solver has an appropriate solution type. The solution type holds all of the
information about the problem which was solved an its solution. If you enabled
`fullSave=true`, then the solver also includes a time-course of the solution
captured at every `saveSteps` steps.

The solution type has a lot of built in functionality to help analysis. Plotting
functionality is provided for each solution type. To plot the solution, simply use

```julia
plot(sol)
```

The plotting function is implemented as a recipe to Plots.jl and as such receives
all of the features of a Plots.jl plot.

Another feature is the `ConvergenceSimulation`s.One can automatically have
DifferentialEquations.jl perform the error analysis by
passing a `ConvergenceSimulation` a vector of solutions, or using one of the provided
`testConvergence` functions. These will give order of convergence estimates and
provide plotting functionality.

## Solution Types

```@docs
DifferentialEquations.FEMSolution
DifferentialEquations.DESolution
DifferentialEquations.SDESolution
DifferentialEquations.ODESolution
DifferentialEquations.StokesSolution
```

## Related Functions

```@docs
DifferentialEquations.appxTrue!
DifferentialEquations.FEMSolutionTS
```

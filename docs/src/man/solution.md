# The Solution Type

Each solver has an appropriate solution type. The solution type holds all of the
information about the problem which was solved an its solution. If you enabled
`save_timeseries=true`, then the solver also includes a time-course of the solution
captured at every `timeseries_steps` steps.

The solution type has a lot of built in functionality to help analysis. For example,
it has an array interface for accessing the values. We can use

```julia
sol[i]
```

to access the value at timestep `i` (if the timeseres was saved), and

```julia
sol.t[i]
```

to access the value of `t` at timestep `i`. The final value of the simulation,
which is always saved, is saved to

```julia
sol.u
```

If the analytical solution, we also have

```julia
sol.u_analytic # final value
sol.timeseries_analytic # timeseries of analytical solution, saved if save_timesseries == true
```


Plotting functionality is provided for each solution type. To plot the solution, simply use

```julia
plot(sol)
```

The plotting function is implemented as a recipe to Plots.jl and as such receives
all of the features of a Plots.jl plot.

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

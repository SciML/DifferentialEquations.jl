# Defining an ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE

```math
du/dt = f(t,u)
```

`f` should be specified as `f(t,u)` and `u₀` should be an AbstractArray whose
geometry matches the desired geometry of `u`. Note that we are not limited to
numbers or vectors for `u₀`, one is allowed to provide `u₀` as arbitrary
matrices / higher dimension tensors as well.

## Problem Type

```@docs
DifferentialEquations.ODEProblem
```

## Example Problems

Examples problems can be found in <a href="https://github.com/ChrisRackauckas/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl">src/premades/premade_problems.jl</a>

```@docs
DifferentialEquations.twoDimlinearODEExample
DifferentialEquations.linearODEExample
```

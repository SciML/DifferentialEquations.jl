# Defining an ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE

```math
du/dt = f(u,t)dt
```

`f` should be specified as `f(u,t)` and `u₀` should be an AbstractArray whose
geometry matches the desired geometry of `u`. Note that we are not limited to
numbers or vectors for `u₀`, one is allowed to provide `u₀` as arbitrary
matrices / higher dimension tensors as well.

## Problem Type

```@docs
DifferentialEquations.ODEProblem
```

## Example Problems

```@docs
DifferentialEquations.twoDimlinearODEExample
DifferentialEquations.linearODEExample
```

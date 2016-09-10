# Defining an ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE

```math
\frac{du}{dt} = f(t,u)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`),and `u₀` should be an AbstractArray
(or number) whose geometry matches the desired geometry of `u`. Note that we are
not limited to numbers or vectors for `u₀`, one is allowed to provide `u₀` as
arbitrary matrices / higher dimension tensors as well.

## Problem Type

```@docs
DifferentialEquations.ODEProblem
```

## Example Problems

Examples problems can be found in [src/premades/premade_problems.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
prob = prob_ode_linear
sol = solve(prob,[0;1])
```

```@docs
DifferentialEquations.prob_ode_linear
DifferentialEquations.prob_ode_2Dlinear
DifferentialEquations.prob_ode_bigfloatlinear
DifferentialEquations.prob_ode_bigfloat2Dlinear
DifferentialEquations.prob_ode_large2Dlinear
DifferentialEquations.prob_ode_2Dlinear_notinplace
DifferentialEquations.prob_ode_threebody
DifferentialEquations.prob_ode_pleides
DifferentialEquations.prob_ode_vanderpol
DifferentialEquations.prob_ode_vanderpol_stiff
DifferentialEquations.prob_ode_rober
DifferentialEquations.prob_ode_lorenz
DifferentialEquations.prob_ode_rigidbody
```

# Defining a SDE Problem

To define an SDE Problem, you simply need to give the forcing function ``f``,
the noise function `σ`, and the initial condition ``u₀`` which define an SDE

```math
du = f(t,u)dt + Σσᵢ(t,u)dWⁱ
```

`f` and `σ` should be specified as `f(t,u)` and  `σ(t,u)` respectively, and `u₀`
should be an AbstractArray whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`, one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well. A vector
of `σ`s can also be defined to determine an SDE of higher Ito dimension.

## Problem Type

```@docs
DifferentialEquations.SDEProblem
```

## Example Problems

Examples problems can be found in [src/premades/premade_problems.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl)

```@docs
DifferentialEquations.prob_sde_linear
DifferentialEquations.prob_sde_2Dlinear
DifferentialEquations.prob_sde_wave
DifferentialEquations.prob_sde_lorenz
DifferentialEquations.prob_sde_cubic
DifferentialEquations.prob_sde_additive
DifferentialEquations.prob_sde_additivesystem
```

# Defining a FEM Problem

Below are the definitions of the types which specify problems. Some general notes are:

* (t,x) vs (t,x,y): Mathematically one normally specifies equations in 2D as ``f(t,x,y)``.
  However, in this code we use `x` as a vector. Thus you can think of ``x``=`x[:,1]` and
  ``y``=`x[:,2]`. Thus input equations are of the form `f(x,t)` no matter the dimension.
  If time is not included in the problem (for example, a Poisson equation problem),
  then we use `f(x)`. An example is the equation ``u(x,y)= sin(2πx)cos(2πy)/(8π^2)``
  would be specified as `sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)`.
* Linearity: If the equation has linear term, they are specified with functions
  `f(t,x)`. If it is nonlinear, it is specified with functions `f(t,x,u)`. The boundary
  conditions are always `(t,x)`
* Stochastic: By default the equation is deterministic. For each equation, one can
  specify a σ term which adds a stochastic ``σ(t,x,u)dW_t`` term to the equation
  (or with ``σ(t,x)dW_t`` if linear, must match `f`). ``dW_t`` corresponds to the type
  of noise which is chosen. By default this is space-time Gaussian white noise.

## Poisson Equation Problem

```@docs
PoissonProblem
```

## Heat Equation Problem

```@docs
HeatProblem
```

## Example Problems

Examples problems can be found in [src/premades/premade_problems.jl](https://github.com/ChrisRackauckas/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl)


### Poisson Equation

```@docs
DifferentialEquations.prob_poisson_birthdeathinteractingsystem
DifferentialEquations.prob_poisson_noisywave
DifferentialEquations.prob_poisson_birthdeathsystem
DifferentialEquations.prob_poisson_wave
DifferentialEquations.prob_poisson_birthdeath
```

### Heat Equation

```@docs
DifferentialEquations.prob_femheat_birthdeathsystem
DifferentialEquations.prob_femheat_birthdeathinteractingsystem
DifferentialEquations.prob_femheat_diffuse
DifferentialEquations.prob_femheat_stochasticbirthdeath
DifferentialEquations.prob_femheat_moving
DifferentialEquations.heatProblemExample_gierermeinhardt
DifferentialEquations.heatProblemExample_grayscott
DifferentialEquations.prob_femheat_pure
DifferentialEquations.prob_femheat_diffusionconstants
DifferentialEquations.prob_femheat_birthdeath
```

## Related Functions

```@docs
DifferentialEquations.DEProblem
```

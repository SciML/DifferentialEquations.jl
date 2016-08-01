# Defining a FEM Problem

Below are the definitions of the types which specify problems. Some general notes are:

* (x,t) vs (x,y,t): Mathematically one normally specifies equations in 2D as ``f(x,y,t)``.
  However, in this code we use `x` as a vector. Thus you can think of ``x``=`x[:,1]` and
  ``y``=`x[:,2]`. Thus input equations are of the form `f(x,t)` no matter the dimension.
  If time is not included in the problem (for example, a Poisson equation problem),
  then we use `f(x)`. An example is the equation ``u(x,y)= sin(2πx)cos(2πy)/(8π^2)``
  would be specified as `sol(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)`.
* Linearity: If the equation has linear term, they are specified with functions
  `f(x,t)`. If it is nonlinear, it is specified with functions `f(u,x,t)`. The boundary
  conditions are always `(x,t)`
* Stochastic: By default the equation is deterministic. For each equation, one can
  specify a σ term which adds a stochastic ``σ(u,x,t)dW_t`` term to the equation
  (or with ``σ(x,t)dW_t`` if linear, must match `f`). ``dW_t`` corresponds to the type
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

Examples problems can be found in <a href="https://github.com/ChrisRackauckas/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl">src/premades/premade_problems.jl</a>


### Poisson Equation

```@docs
DifferentialEquations.poissonProblemExample_wave
DifferentialEquations.poissonProblemExample_noisyWave
DifferentialEquations.poissonProblemExample_birthdeath
DifferentialEquations.poissonProblemExample_birthdeathinteractingsystem
DifferentialEquations.poissonProblemExample_birthdeathsystem
```

### Heat Equation

```@docs
DifferentialEquations.heatProblemExample_diffuse
DifferentialEquations.heatProblemExample_pure
DifferentialEquations.heatProblemExample_moving
DifferentialEquations.heatProblemExample_birthdeath
DifferentialEquations.heatProblemExample_stochasticbirthdeath
DifferentialEquations.heatProblemExample_birthdeathinteractingsystem
DifferentialEquations.heatProblemExample_gierermeinhardt
DifferentialEquations.heatProblemExample_birthdeathsystem
DifferentialEquations.heatProblemExample_diffusionconstants
DifferentialEquations.heatProblemExample_grayscott
```

## Related Functions

```@docs
DifferentialEquations.DEProblem
```

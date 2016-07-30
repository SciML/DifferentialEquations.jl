# Solver Helpers

This package includes the documentation for the helper functions for the various
solvers.

## ODE Solver Extras

```@docs
Base.length(::ExplicitRK)
DifferentialEquations.ExplicitRK
DifferentialEquations.ODE_DEFAULT_TABLEAU
DifferentialEquations.constructCashKarp
DifferentialEquations.constructRalston
DifferentialEquations.constructDormandPrince
DifferentialEquations.constructDormandPrince8
DifferentialEquations.constructBogakiShampine
DifferentialEquations.constructHuen
DifferentialEquations.constructRKF
DifferentialEquations.constructRKF8
```

## SDE Solver Extras

```@docs
DifferentialEquations.monteCarloSim
DifferentialEquations.RosslerSRI
DifferentialEquations.RosslerSRA
DifferentialEquations.constructSRA1
DifferentialEquations.constructSRIW1
DifferentialEquations.checkSRAOrder
DifferentialEquations.checkSRIOrder
```

## Stationary Stokes

```@docs
DifferentialEquations.GSδq!
DifferentialEquations.GSu!
DifferentialEquations.calc_rp!
DifferentialEquations.update_p!
DifferentialEquations.update_v!
DifferentialEquations.uzawa_p!
DifferentialEquations.stokes_restriction
DifferentialEquations.stokes_prolongation
DifferentialEquations.update_u!
DifferentialEquations.GSv!
```

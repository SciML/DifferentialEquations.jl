# Information on Solvers

## Ordinary Differential Equation Solvers

### ODE

```@docs
solve(::ODEProblem,::Number,::Number)
```

### SDE

```@docs
solve(::SDEProblem,::Number,::Number)
```

## Finite Difference Method Solvers

### Stokes Equation

```@docs
solve(::StokesProblem,::FDMMesh)
```

## Finite Element Method Solvers

```@docs
solve(::FEMmesh,::PoissonProblem)
```

```@docs
solve(::FEMmesh,::HeatProblem)
```

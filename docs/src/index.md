# DifferentialEquations.jl Documentation

This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), and differential delay equations. For ordinary differential equation solvers, see [ODE.jl](https://github.com/JuliaLang/ODE.jl)

This package is for efficient and parallel implementations of research-level algorithms, many of which are quite recent. These algorithms aim to be optimized for HPC applications, including the use of GPUs, Xeon Phis, and multi-node parallelism.

## Functions

{docs}
DifferentialEquations

{docs}
findboundary(elem,bdFlag=[])

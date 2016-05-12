# DifferentialEquations.jl

[![Join the chat at https://gitter.im/ChrisRackauckas/DiffEq.jl](https://badges.gitter.im/ChrisRackauckas/DiffEq.jl.svg)](https://gitter.im/ChrisRackauckas/DiffEq.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), and differential delay equations. For ordinary differential equation solvers, see [ODE.jl](https://github.com/JuliaLang/ODE.jl)

This package is for efficient and parallel implementations of research-level algorithms, many of which are quite recent. These algorithms aim to be optimized for HPC applications, including the use of GPUs, Xeon Phis, and multi-node parallelism.

Currently, finite element solvers for the Poisson and Heat Equations are supplied. Mesh generation tools currently only work for squares, though the solvers will work on general meshes if they are provided. The mesh layout follows the format of [iFEM](http://www.math.uci.edu/~chenlong/programming.html) and many subroutines in the finite element solver are based off of iFEM algorithms.

# Using the package

To install the package, use `Pkg.clone("https://github.com/ChrisRackauckas/DifferentialEquations.jl")` inside the Julia REPL. To load the package, use the command `using DifferentialEquations`. To understand the package in more detail, check out the examples in [test/](test/).

# Current Supported Equations
* PDE Solvers
  * Finite Element Solvers
    * Linear Poisson Equation
    * Semi-linear Poisson Equation
    * Linear Heat Equation
    * Semi-linear Heat Equation (aka Reaction-Diffusion Equation)

# Roadmap
* SODE Solvers
  * Euler-Maruyama
  * Milstein
  * Rossler-SRK
  * Adaptive-SRK
* PDE Solvers
  * Finite difference solvers:
    * Semi-linear Heat Equation (Reaction-Diffusion Equation)
    * Semi-linear Poisson Equation
    * Wave Equation
    * Transport Equation
    * Stokes Equation
* SPDE Solvers
  * Euler-Maruyama
  * IIF-Maruyama
  * IIF-Milstein

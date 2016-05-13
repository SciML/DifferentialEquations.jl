
<a id='DifferentialEquations.jl-Documentation-1'></a>

# DifferentialEquations.jl Documentation


This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), and differential delay equations. For ordinary differential equation solvers, see [ODE.jl](https://github.com/JuliaLang/ODE.jl)


This package is for efficient and parallel implementations of research-level algorithms, many of which are quite recent. These algorithms aim to be optimized for HPC applications, including the use of GPUs, Xeon Phis, and multi-node parallelism. With the easy to use plot/convergence testing algorithms, this package also provides a good sandbox for developing novel numerical schemes.


If at any time you run into documentation that is incomplete/confusing, please contact me via the Gitter channel and I will clear it up!


<a id='Using-the-Package-1'></a>

## Using the Package


To install the package, use the following command inside the Julia REPL:


```julia
Pkg.clone("https://github.com/ChrisRackauckas/DifferentialEquations.jl")
```


To load the package, use the command


```julia
using DifferentialEquations
```


To understand the package in more detail, check out the following tutorials. Example codes for the latest features can be found in [test/](test/). Due to the fast development speed of the package, it is recommended you checkout the main branch:


```
Pkg.checkout("DifferentialEquations")
```


This will keep your local repository up to date with the latest changes.


<a id='Tutorials-1'></a>

## Tutorials


The following tutorials will introduce you to the functionality of DifferentialEquations.jl

- [Poisson Equation Finite Element Method Example](tutorials/femPoisson.md#Poisson-Equation-Finite-Element-Method-Example-1)
- [Heat Equation Finite Element Method Example](tutorials/femHeat.md#Heat-Equation-Finite-Element-Method-Example-1)
- [Stochastic Finite Element Examples](tutorials/femStochastic.md#Stochastic-Finite-Element-Examples-1)
    - [Finite Element Stochastic Poisson Equation](tutorials/femStochastic.md#Finite-Element-Stochastic-Poisson-Equation-1)
    - [Finite Element Stochastic Heat Equation](tutorials/femStochastic.md#Finite-Element-Stochastic-Heat-Equation-1)

<a id='Manual-1'></a>

## Manual

- [Overview of DifferentialEquations.jl Usage](man/overview.md#Overview-of-DifferentialEquations.jl-Usage-1)
- [Defining a Problem](man/problem.md#Defining-a-Problem-1)
    - [Poisson Equation Problem](man/problem.md#Poisson-Equation-Problem-1)
    - [Heat Equation Problem](man/problem.md#Heat-Equation-Problem-1)
    - [Example Problems](man/problem.md#Example-Problems-1)
    - [Related Functions](man/problem.md#Related-Functions-1)
- [Information on Solvers](man/solvers.md#Information-on-Solvers-1)
    - [Finite Element Method Solvers](man/solvers.md#Finite-Element-Method-Solvers-1)
- [The Solution Type](man/solution.md#The-Solution-Type-1)
    - [Related Functions](man/solution.md#Related-Functions-1)
- [Plot Functions](man/plot.md#Plot-Functions-1)
- [Convergence Simulations](man/convergence.md#Convergence-Simulations-1)

<a id='Internal-Documentation-1'></a>

## Internal Documentation

- [Internal Finite Element Tools](internals/femTools.md#Internal-Finite-Element-Tools-1)
    - [General](internals/femTools.md#General-1)
    - [Mesh Tools](internals/femTools.md#Mesh-Tools-1)
    - [Solver Tools](internals/femTools.md#Solver-Tools-1)
    - [Error Tools](internals/femTools.md#Error-Tools-1)
- [Extra Functions](internals/extras.md#Extra-Functions-1)

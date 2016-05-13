
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


To understand the package in more detail, check out the following tutorials. Example codes for the latest features can be found in [test/](test/). Note that for many of the examples in the test folder, you may wish to run them at lower Δx or Δt. These values were taken to be large in order make unit tests run faster! Due to the fast development speed of the package, it is recommended you checkout the main branch:


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
- [Common Terms](man/commonTerms.md#Common-Terms-1)
- [Defining a Problem](man/problem.md#Defining-a-Problem-1)
    - [Poisson Equation Problem](man/problem.md#Poisson-Equation-Problem-1)
    - [Heat Equation Problem](man/problem.md#Heat-Equation-Problem-1)
    - [Example Problems](man/problem.md#Example-Problems-1)
    - [Related Functions](man/problem.md#Related-Functions-1)
- [Mesh Generation](man/mesh.md#Mesh-Generation-1)
    - [Mesh Specification](man/mesh.md#Mesh-Specification-1)
    - [Mesh Type](man/mesh.md#Mesh-Type-1)
    - [Mesh Generation Functions](man/mesh.md#Mesh-Generation-Functions-1)
- [Information on Solvers](man/solvers.md#Information-on-Solvers-1)
    - [Finite Element Method Solvers](man/solvers.md#Finite-Element-Method-Solvers-1)
- [The Solution Type](man/solution.md#The-Solution-Type-1)
    - [Related Functions](man/solution.md#Related-Functions-1)
- [Plot Functions](man/plot.md#Plot-Functions-1)
    - [Related Functions](man/plot.md#Related-Functions-1)
- [Convergence Simulations](man/convergence.md#Convergence-Simulations-1)
    - [The ConvergenceSimulation Type](man/convergence.md#The-ConvergenceSimulation-Type-1)
    - [Related Functions](man/convergence.md#Related-Functions-1)
    - [Plot Functions](man/convergence.md#Plot-Functions-1)

<a id='Internal-Documentation-1'></a>

## Internal Documentation

- [Internal Finite Element Tools](internals/femTools.md#Internal-Finite-Element-Tools-1)
    - [General](internals/femTools.md#General-1)
    - [Mesh Tools](internals/femTools.md#Mesh-Tools-1)
    - [Solver Tools](internals/femTools.md#Solver-Tools-1)
    - [Error Tools](internals/femTools.md#Error-Tools-1)
- [Extra Functions](internals/extras.md#Extra-Functions-1)

<a id='Index-1'></a>

## Index

- [`DifferentialEquations.CG2`](internals/extras.md#DifferentialEquations.CG2)
- [`DifferentialEquations.quadfbasis2`](internals/extras.md#DifferentialEquations.quadfbasis2)
- [`DifferentialEquations`](internals/femTools.md#DifferentialEquations)
- [`DifferentialEquations.CFLμ`](internals/femTools.md#DifferentialEquations.CFLμ)
- [`DifferentialEquations.CFLν`](internals/femTools.md#DifferentialEquations.CFLν)
- [`DifferentialEquations.accumarray`](internals/femTools.md#DifferentialEquations.accumarray)
- [`DifferentialEquations.assemblematrix`](internals/femTools.md#DifferentialEquations.assemblematrix)
- [`DifferentialEquations.getH1error`](internals/femTools.md#DifferentialEquations.getH1error)
- [`DifferentialEquations.getL2error`](internals/femTools.md#DifferentialEquations.getL2error)
- [`DifferentialEquations.gradbasis`](internals/femTools.md#DifferentialEquations.gradbasis)
- [`DifferentialEquations.gradu`](internals/femTools.md#DifferentialEquations.gradu)
- [`DifferentialEquations.meshgrid`](internals/femTools.md#DifferentialEquations.meshgrid)
- [`DifferentialEquations.quadfbasis`](internals/femTools.md#DifferentialEquations.quadfbasis)
- [`DifferentialEquations.quadpts`](internals/femTools.md#DifferentialEquations.quadpts)
- [`Base.length`](man/convergence.md#Base.length-Tuple{DifferentialEquations.ConvergenceSimulation})
- [`DifferentialEquations.ConvergenceSimulation`](man/convergence.md#DifferentialEquations.ConvergenceSimulation)
- [`DifferentialEquations.conv_ests`](man/convergence.md#DifferentialEquations.conv_ests)
- [`DifferentialEquations.convplot_fullΔt`](man/convergence.md#DifferentialEquations.convplot_fullΔt)
- [`DifferentialEquations.convplot_fullΔx`](man/convergence.md#DifferentialEquations.convplot_fullΔx)
- [`DifferentialEquations.convplot_h1vsΔt`](man/convergence.md#DifferentialEquations.convplot_h1vsΔt)
- [`DifferentialEquations.convplot_h1vsΔx`](man/convergence.md#DifferentialEquations.convplot_h1vsΔx)
- [`DifferentialEquations.convplot_l2vsΔt`](man/convergence.md#DifferentialEquations.convplot_l2vsΔt)
- [`DifferentialEquations.convplot_l2vsΔx`](man/convergence.md#DifferentialEquations.convplot_l2vsΔx)
- [`DifferentialEquations.convplot_maxvsΔt`](man/convergence.md#DifferentialEquations.convplot_maxvsΔt)
- [`DifferentialEquations.convplot_maxvsΔx`](man/convergence.md#DifferentialEquations.convplot_maxvsΔx)
- [`DifferentialEquations.convplot_node2vsΔt`](man/convergence.md#DifferentialEquations.convplot_node2vsΔt)
- [`DifferentialEquations.convplot_node2vsΔx`](man/convergence.md#DifferentialEquations.convplot_node2vsΔx)
- [`DifferentialEquations.FEMmesh`](man/mesh.md#DifferentialEquations.FEMmesh)
- [`DifferentialEquations.fem_squaremesh`](man/mesh.md#DifferentialEquations.fem_squaremesh)
- [`DifferentialEquations.findboundary`](man/mesh.md#DifferentialEquations.findboundary)
- [`DifferentialEquations.notime_squaremesh`](man/mesh.md#DifferentialEquations.notime_squaremesh)
- [`DifferentialEquations.parabolic_squaremesh`](man/mesh.md#DifferentialEquations.parabolic_squaremesh)
- [`DifferentialEquations.setboundary`](man/mesh.md#DifferentialEquations.setboundary)
- [`DifferentialEquations.convplot`](man/plot.md#DifferentialEquations.convplot)
- [`DifferentialEquations.showmesh`](man/plot.md#DifferentialEquations.showmesh)
- [`DifferentialEquations.solplot`](man/plot.md#DifferentialEquations.solplot)
- [`DifferentialEquations.solplot_animation`](man/plot.md#DifferentialEquations.solplot_animation)
- [`DifferentialEquations.solplot_appx`](man/plot.md#DifferentialEquations.solplot_appx)
- [`DifferentialEquations.solplot_appxvstrue`](man/plot.md#DifferentialEquations.solplot_appxvstrue)
- [`DifferentialEquations.HeatProblem`](man/problem.md#DifferentialEquations.HeatProblem)
- [`DifferentialEquations.PdeProblem`](man/problem.md#DifferentialEquations.PdeProblem)
- [`DifferentialEquations.PoissonProblem`](man/problem.md#DifferentialEquations.PoissonProblem)
- [`DifferentialEquations.heatProblemExample_birthdeath`](man/problem.md#DifferentialEquations.heatProblemExample_birthdeath)
- [`DifferentialEquations.heatProblemExample_diffuse`](man/problem.md#DifferentialEquations.heatProblemExample_diffuse)
- [`DifferentialEquations.heatProblemExample_moving`](man/problem.md#DifferentialEquations.heatProblemExample_moving)
- [`DifferentialEquations.heatProblemExample_pure`](man/problem.md#DifferentialEquations.heatProblemExample_pure)
- [`DifferentialEquations.heatProblemExample_stochasticbirthdeath`](man/problem.md#DifferentialEquations.heatProblemExample_stochasticbirthdeath)
- [`DifferentialEquations.poissonProblemExample_noisyWave`](man/problem.md#DifferentialEquations.poissonProblemExample_noisyWave)
- [`DifferentialEquations.poissonProblemExample_wave`](man/problem.md#DifferentialEquations.poissonProblemExample_wave)
- [`DifferentialEquations.FEMSolution`](man/solution.md#DifferentialEquations.FEMSolution)
- [`DifferentialEquations.PdeSolution`](man/solution.md#DifferentialEquations.PdeSolution)
- [`DifferentialEquations.appxTrue!`](man/solution.md#DifferentialEquations.appxTrue!)
- [`DifferentialEquations.fem_solveheat`](man/solvers.md#DifferentialEquations.fem_solveheat)
- [`DifferentialEquations.fem_solvepoisson`](man/solvers.md#DifferentialEquations.fem_solvepoisson)

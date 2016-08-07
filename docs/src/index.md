# DifferentialEquations.jl Documentation

This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include ordinary differential equations (ODEs), stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), differential algebraic equations, and differential delay equations. It includes algorithms from very recent research, as well as algorithms optimized for HPC applications. It integrates with the Julia package sphere, for example using Juno's progress meter, and wraps other differential equation solvers so that many different methods for solving the equations can be accessed by simply switching a keyword argument.

All of the algorithms are thoroughly tested to ensure accuracy. Convergence tests are included in the [test/](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/test) folder. The algorithms were also tested to show correctness with nontrivial behavior such as Turing morphogenesis. Example IJulia notebooks
[can be found in the examples folder](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/examples). If you find any example where there seems
to be an error, please open an issue.

If you have any questions, or just want to chat about solvers/using the package, please feel free to message me in the [Gitter channel](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge). For bug reports, feature requests, etc., please submit an issue. If you're interested in contributing, please see the [Contributor's Guide](/internals/contributors_guide).

## Using the Package

To install the package, use the following command inside the Julia REPL:
```julia
Pkg.add("DifferentialEquations")
```

For all of the latest features, switch to the master branch via:

```julia
Pkg.checkout("DifferentialEquations")
```

To load the package, use the command:

```julia
using DifferentialEquations
```

To understand the package in more detail, check out the following tutorials in the manual. Examples
IJulia notebooks using DifferentialEquations can be found [in the examples folder](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/examples).
Codes for the latest features can be found in [test/](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/test).

For the most up to date on using the package information, please contact me [via the repository Gitter](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl)
or [read the latest documentation](http://chrisrackauckas.github.io/DifferentialEquations.jl/latest/)

## Supported Equations

For PDEs, one can optionally specify a noise equation. The solvers currently have
stochastic variants for handling Gaussian Space-time white noise SPDEs.

* ODEs
* SODEs
* (Stochastic) PDEs

  * Linear Poisson Equation
  * Semi-linear Poisson Equation
  * Linear Heat Equation
  * Semi-linear Heat Equation (aka Reaction-Diffusion Equation)
  * Stationary Stokes Equation

## Implemented Solvers

For PDEs, [method] denotes an additional version for handling stochastic partial
differential equations. SDE solvers and ODE solvers take in general sized inputs.
For example, if u₀ is a matrix (and your problem functions are designed to work
with matrices), then the solver will use the matrices without error.

* ODEs

  * Optimized Explicit Solvers

    * Euler
    * Midpoint Method
    * RK4
    * Feagin's Order 10/8 Method
    * Feagin's Order 12/10 Method
    * Feagin's Order 14/12 Method

  * General Explicit (Adaptive) Runge-Kutta Methods

    * Huen's Method
    * Cash-Karp
    * Runge-Kutta-Fehlberg (RKF) 4/5
    * Ralston's Method
    * Bogaki-Shampine
    * Dormand-Prince 4/5
    * Runge-Kutta-Fehlberg (RKF) 7/8
    * Dormand-Prince 7/8

  * Stiff Solvers. Requires [NLsolve.jl](https://github.com/EconForge/NLsolve.jl) and optionally [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). See [Conditional Dependencies](/man/conditional_dependencies).

    * Implicit Euler
    * Trapezoidal
    * Rosenbrock32

  * Wrappers for ODEInterface.jl. See [Conditional Dependencies](/man/conditional_dependencies).

    * dorpi5 - Hairer's DP5(4)
    * dop853 - Hairer's DP8(5,3)
    * odex - Extrapolation algorithm based on explicit midpoint rule
    * radau5 - Implicit Runge-Kutta order 5
    * radau - Implicit Runge-Kutta variable order 5-13
    * seulex - Extrapolation based on linear implicit Euler

  * Wrappers for ODE.jl. See [Conditional Dependencies](/man/conditional_dependencies).

    * ode23 - Bogacki-Shampine's method
    * ode45 - Dormand-Prince  4/5
    * ode78 - Runge-Kutta-Fehlberg  7/8
    * ode23s - Rosenbrock method 2/3
    * ode1 - Forward Euler
    * midpoint - Midpoint method
    * ode2_heun - Huen's method
    * ode4 - RK4
    * ode45_fe - Runge-Kutta-Fehlberg 4/5


* SODEs

  * Euler-Maruyama
  * Milstein
  * Rossler-SRK


* Finite Element Solvers (Stochastic) PDEs

  * Semilinear Poisson Equation

    * See implicit solvers

  * Semilinear Heat Equation (Reaction-Diffusion)

    * Forward Euler [Maruyama]
    * Backward Euler [Maruyama]
    * Semi-implicit Crank-Nicholson [Maruyama]
    * Semi-implicit Backward Euler [Maruyama]

  * Linear Heat Equation

    * Forward Euler [Maruyama]
    * Backward Euler [Maruyama]
    * Crank-Nicholson [Maruyama]


* Implicit Solvers

  * Direct
  * Factorizations (LU, Cholesky, QR, SVD)
  * Conjugate-Gradient (CG)
  * GMRES

## Roadmap

* ODE Solvers

  * Stabilized stiff - ROCK2 and ROCK4

* SODE Solvers

  * Adaptive-SRK

* Finite difference solvers

  * Semi-linear Heat Equation (Reaction-Diffusion Equation)
  * Semi-linear Poisson Equation
  * Wave Equation
  * Transport Equation

* Stochastic PDE Solvers

  * Implicit Integration Factor (IIF) Maruyama
  * Implicit Integration Factor (IIF) Milstein

* DDE Solvers

  * Wrap RETARD and RADAR5
  * Implement standard Runge-Kutta DDE solvers

* Algebraic differential equations

  * Implement standard solvers and add to ODEProblem type

* Linear Solvers

  * Finite Difference Geometric Multigrids
  * Algebraic Multigrids via pyAMG

* Performance

  * Improve FEM performance
  * Implement threaded versions
  * Test ParallelAccelerator.jl on solvers
  * Add Xeon Phi / GPU variants

* Misc

  * Davie-Gaines convergence analysis
  * Add benchmarking tools
  * Improve MonteCarloSimulation


## IJulia Notebook Tutorials

If you have [IJulia](https://github.com/JuliaLang/IJulia.jl) installed, you can access
extra tutorials in the supplied IJulia notebooks via:

```julia
using IJulia
cd(Pkg.dir("DifferentialEquations")*"/examples")
notebook()
```

Otherwise, these notebooks can be viewed [via the Github repository](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/examples)
(note that Github renders them slightly incorrectly, so it will look better in IJulia!).

## Tutorials

The following tutorials will introduce you to the functionality of DifferentialEquations.jl
More examples can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/ChrisRackauckas/DifferentialEquations.jl/tree/master/examples).

```@contents
Pages = [
    "tutorials/ode_example.md",
    "tutorials/sde_example.md",
    "tutorials/fempoisson_example.md",
    "tutorials/femheat_example.md",
    "tutorials/femstochastic_example.md"
    ]
Depth = 2
```

## Solver Options

These pages describe the options available in the solvers.

```@contents
Pages = [
  "solvers/ode_solve.md",
  "solvers/sde_solve.md",
  "solvers/fempoisson_solve.md",
  "solvers/femheat_solve.md",
  "solvers/fdmstokes_solve.md"
]
Depth = 2
```

## Manual

```@contents
Pages = [
    "man/overview.md",
    "man/ODEProblem.md",
    "man/SDEProblem.md",
    "man/FEMProblem.md",
    "man/StokesProblem.md",
    "man/mesh.md",
    "man/solution.md",
    "man/plot.md",
    "man/convergence.md",
    "man/conditional_dependencies.md"
    "man/progress_bar.md"
]
Depth = 2
```

## Internal Documentation

```@contents
Pages = [
  "internals/contributors_guide.md",
  "internals/fem_tools.md",
  "internals/extras.md",
  "internals/solver_helpers.md"
]
Depth = 2
```

## Index

```@index
```

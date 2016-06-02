# DifferentialEquations.jl Documentation

This is a package for solving numerically solving differential equations in Julia by Chris Rackauckas. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include ordinary differential equations (ODEs), stochastic ordinary differential equations (SODEs or SDEs), stochastic partial differential equations (SPDEs), partial differential equations (with both finite difference and finite element methods), and differential delay equations.

This package is for efficient and parallel implementations of research-level algorithms, many of which are quite recent. These algorithms aim to be optimized for HPC applications, including the use of GPUs, Xeon Phis, and multi-node parallelism. With the easy to use plot/convergence testing algorithms, this package also provides a good sandbox for developing novel numerical schemes. Since this package is designed for long computations, one of the features of this package is the existence of tools for inspecting a long calculation. These include optional printing and, if the user is using Juno, a progress meter (with time estimates once implemented on Juno's end). All of the methods have associated tests to ensure that the theoretical order of accuracy is attained.

If you have any questions, or just want to chat about solvers/using the package, please feel free to message me in the Gitter channel. For bug reports, feature requests, etc., please submit an issue.

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

To understand the package in more detail, check out the following tutorials. Example
codes for the latest features can be found in [test/](https://github.com/ChrisRackauckas/DifferentialEquations.jl/test/). Note that for many of
the examples in the test folder, you may wish to run them at lower Δx or Δt.
These values were taken to be large in order make unit tests run faster!

For the most up to date on using the package information, please contact me [via the repository Gitter](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl)
or [read the latest documentation](http://chrisrackauckas.github.io/DifferentialEquations.jl/latest/)

## Tutorials

The following tutorials will introduce you to the functionality of DifferentialEquations.jl

```@contents
Pages = [
    "tutorials/ode.md",
    "tutorials/sde.md",
    "tutorials/femPoisson.md",
    "tutorials/femHeat.md",
    "tutorials/femStochastic.md"
    ]
Depth = 2
```

## Manual

```@contents
Pages = [
    "man/overview.md",
    "man/odeProblem.md",
    "man/sdeProblem.md",
    "man/femProblem.md",
    "man/stokesProblem.md",
    "man/mesh.md",
    "man/solvers.md",
    "man/solution.md",
    "man/plot.md",
    "man/convergence.md",
]
Depth = 2
```

## Internal Documentation

```@contents
Pages = [
  "internals/femTools.md",
  "internals/extras.md",
  "internals/solverHelpers.md"
]
Depth = 2
```

## Index

```@index
```

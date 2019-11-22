# DifferentialEquations.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Travis](https://travis-ci.org/JuliaDiffEq/DifferentialEquations.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DifferentialEquations.jl)
[![AppVoyer](https://ci.appveyor.com/api/projects/status/1smlr9ryfqfx1ear?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/differentialequations-jl-1sx90)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://docs.juliadiffeq.org/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://docs.juliadiffeq.org/dev/)
[![DOI](https://zenodo.org/badge/58516043.svg)](https://zenodo.org/badge/latestdoi/58516043)

This is a suite for numerically solving differential equations written in Julia 
and available for use in Julia, Python, and R. The
purpose of this package is to supply efficient Julia implementations of solvers
for various differential equations. Equations within the realm of this package
include:

- Discrete equations (function maps, discrete stochastic (Gillespie/Markov)
  simulations)
- Ordinary differential equations (ODEs)
- Split and Partitioned ODEs (Symplectic integrators, IMEX Methods)
- Stochastic ordinary differential equations (SODEs or SDEs)
- Random differential equations (RODEs or RDEs)
- Differential algebraic equations (DAEs)
- Delay differential equations (DDEs)
- Mixed discrete and continuous equations (Hybrid Equations, Jump Diffusions)
- (Stochastic) partial differential equations ((S)PDEs) (with both finite
  difference and finite element methods)

The well-optimized DifferentialEquations solvers benchmark as the some of the fastest
implementations, using classic algorithms and ones from recent research which
routinely outperform the "standard" C/Fortran methods, and include algorithms
optimized for high-precision and HPC applications. At the same time, it wraps
the classic C/Fortran methods, making it easy to switch over to them whenever
necessary. Solving differential equations with different methods from 
different languages and packages can be done by changing one line of code, 
allowing for easy benchmarking to ensure you are using the fastest method possible.

DifferentialEquations.jl integrates with the Julia package sphere with: 

- GPU accleration through CUDAnative.jl and CuArrays.jl
- Automated sparsity detection with [SparsityDetection.jl](https://github.com/JuliaDiffEq/SparsityDetection.jl)
- Automatic Jacobian coloring with [SparseDiffTools.jl](https://github.com/JuliaDiffEq/SparseDiffTools.jl), allowing for fast solutions
  to problems with sparse or structured (Tridiagonal, Banded, BlockBanded, etc.) Jacobians
- Allowing the specification of linear solvers for maximal efficiency
- Progress meter integration with the Juno IDE for estimated time to solution 
- Automatic plotting of time series and phase plots 
- Built-in interpolations
- Wraps for common C/Fortran methods like Sundials and Hairer's radau
- Arbitrary precision with BigFloats and Arbfloats
- Arbitrary array types, allowing the definition of differential equations on 
  matrices and distributed arrays
- Unit checked arithmetic with Unitful

Additionally, DifferentialEquations.jl comes with built-in analysis features, including:

- [Forward and adjoint local sensitivity analysis](http://docs.juliadiffeq.org/latest/analysis/sensitivity.html) for fast gradient computations
- [Optimization-based and Bayesian parameter estimation](http://docs.juliadiffeq.org/latest/analysis/parameter_estimation.html)
- Neural differential equations with [DiffEqFlux.jl](https://github.com/JuliaDiffEq/DiffEqFlux.jl)
  for efficient scientific machine learning (scientific ML) and scientific AI.
- [Automatic distributed, multithreaded, and GPU parallelism of ensemble trajectories](http://docs.juliadiffeq.org/latest/features/ensemble.html)
- [Global sensitivity analysis](http://docs.juliadiffeq.org/latest/analysis/global_sensitivity.html)
- [Uncertainty quantification](http://docs.juliadiffeq.org/latest/analysis/uncertainty_quantification.html)

This gives a powerful mixture of speed and productivity features to help you
solve and analyze your differential equations faster.

For information on using the package,
[see the stable documentation](http://docs.juliadiffeq.org/stable/). Use the
[latest documentation](http://docs.juliadiffeq.org/latest/) for the version of
the documentation which contains the un-released features.

All of the algorithms are thoroughly tested to ensure accuracy via convergence
tests. The algorithms are continuously tested to show correctness.
IJulia tutorial notebooks
[can be found at DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl).
Benchmarks
[can be found at DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl).
If you find any equation where there seems to be an error, please open an issue.

If you have any questions, or just want to chat about solvers/using the package,
please feel free to chat in the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge).
For bug reports, feature requests, etc., please submit an issue. If you're
interested in contributing, please see the
[Developer Documentation](https://juliadiffeq.github.io/DiffEqDevDocs.jl/latest/).

## Supporting and Citing

The software in this ecosystem was developed as part of academic research. If you
would like to help support it, please star the repository as such metrics may
help us secure funding in the future. If you use JuliaDiffEq software as part
of your research, teaching, or other activities, we would be grateful if you
could cite our work.
[Please see our citation page for guidelines](http://juliadiffeq.org/citing.html).

--------------------------------

## Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](https://youtu.be/KPEqYtEd-zY)

## Video Introduction

[![Video Introduction to DifferentialEquations.jl](https://user-images.githubusercontent.com/1814174/27973992-e236a9a4-6310-11e7-84af-2b66097cecf9.PNG)](https://youtu.be/75SCMIRlNXM)

## Comparison to MATLAB, R, Julia, Python, C, Mathematica, Maple, and Fortran

<a href="http://www.stochasticlifestyle.com/wp-content/uploads/2019/08/de_solver_software_comparsion.pdf"><img src="http://www.stochasticlifestyle.com/wp-content/uploads/2019/08/de_solver_software_comparsion-1.png" alt="Comparison Of Differential Equation Solver Software" align="middle"/></a>

[See the corresponding blog post](http://www.stochasticlifestyle.com/comparison-differential-equation-solver-suites-matlab-r-julia-python-c-fortran/)

## Example Images

<img src="https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/assets/DifferentialEquations_Example.png" align="middle"  />

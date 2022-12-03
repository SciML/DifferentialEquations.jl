# DifferentialEquations.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/dev/modules/DiffEqDocs/)

[![codecov](https://codecov.io/gh/SciML/DifferentialEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/DifferentialEquations.jl)
[![Build Status](https://github.com/SciML/DifferentialEquations.jl/workflows/CI/badge.svg)](https://github.com/SciML/DifferentialEquations.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

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
- Stochastic differential-algebraic equations (SDAEs)
- Random differential equations (RODEs or RDEs)
- Differential algebraic equations (DAEs)
- Delay differential equations (DDEs)
- Neutral, retarded, and algebraic delay differential equations (NDDEs, RDDEs, and DDAEs)
- Stochastic delay differential equations (SDDEs)
- Experimental support for stochastic neutral, retarded, and algebraic delay differential equations (SNDDEs, SRDDEs, and SDDAEs)
- Mixed discrete and continuous equations (Hybrid Equations, Jump Diffusions)
- (Stochastic) partial differential equations ((S)PDEs) (with both finite
  difference and finite element methods)
  
The well-optimized DifferentialEquations solvers benchmark as some of the fastest
implementations of classic algorithms. It also includes algorithms from recent
research which routinely outperform the "standard" C/Fortran methods, and algorithms
optimized for high-precision and HPC applications. Simultaneously, it wraps
the classic C/Fortran methods, making it easy to switch over to them whenever
necessary. Solving differential equations with different methods from
different languages and packages can be done by changing one line of code,
allowing for easy benchmarking to ensure you are using the fastest method possible.

DifferentialEquations.jl integrates with the Julia package sphere with:

- GPU acceleration through [CUDA.jl](https://cuda.juliagpu.org/stable/) and [DiffEqGPU.jl](https://docs.sciml.ai/DiffEqGPU/dev/)
- Automated sparsity detection with [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)
- Automatic Jacobian coloring with [SparseDiffTools.jl](https://docs.sciml.ai/SparseDiffTools/stable/), allowing for fast solutions
  to problems with sparse or structured (Tridiagonal, Banded, BlockBanded, etc.) Jacobians
- Allowing the specification of linear solvers for maximal efficiency with [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/)
- Progress meter integration with the Visual Studio Code IDE for estimated time to solution
- Automatic plotting of time series and phase plots
- Built-in interpolations
- Wraps for common C/Fortran methods like Sundials and Hairer's radau
- Arbitrary precision with BigFloats and Arbfloats
- Arbitrary array types, allowing the definition of differential equations on
  matrices and distributed arrays
- Unit checked arithmetic with Unitful

Additionally, DifferentialEquations.jl comes with built-in analysis features, including:

Additionally, DifferentialEquations.jl comes with built-in analysis features, including:

- [Forward and Adjoint Sensitivity Analysis (Automatic Differentiation)](https://docs.sciml.ai/SciMLSensitivity/stable/) for fast gradient computations
- [Parameter Estimation and Bayesian Analysis](https://docs.sciml.ai/Overview/stable/highlevels/inverse_problems/)
- Neural differential equations with [DiffEqFlux.jl](https://docs.sciml.ai/DiffEqFlux/stable/)
  for efficient scientific machine learning (scientific ML) and scientific AI.
- Automatic distributed, multithreaded, and GPU [Parallel Ensemble Simulations](https://diffeq.sciml.ai/dev/features/ensemble/)
- [Global Sensitivity Analysis](https://docs.sciml.ai/GlobalSensitivity/stable/)
- [Uncertainty Quantification](https://docs.sciml.ai/Overview/stable/highlevels/uncertainty_quantification/)

This gives a powerful mixture of speed and productivity features to help you
solve and analyze your differential equations faster.

For information on using the package,
[see the stable documentation](https://diffeq.sciml.ai/stable/). Use the
[in-development documentation](https://diffeq.sciml.ai/dev/) for the version of
the documentation which contains the unreleased features.

All of the algorithms are thoroughly tested to ensure accuracy via convergence
tests. The algorithms are continuously tested to show correctness.
IJulia tutorial notebooks
[can be found at DiffEqTutorials.jl](https://github.com/SciML/SciMLTutorials.jl).
Benchmarks
[can be found at DiffEqBenchmarks.jl](https://github.com/SciML/SciMLBenchmarks.jl).
If you find any equation where there seems to be an error, please open an issue.

If you have any questions, or just want to chat about solvers/using the package,
please feel free to chat in the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge).
For bug reports, feature requests, etc., please submit an issue. If you're
interested in contributing, please see the
[Developer Documentation](http://devdocs.sciml.ai/latest/).

## Supporting and Citing

The software in this ecosystem was developed as part of academic research. If you
would like to help support it, please star the repository, as such metrics may
help us secure funding in the future. If you use SciML software as part
of your research, teaching, or other activities, we would be grateful if you
could cite our work.
[Please see our citation page for guidelines](http://sciml.ai/citing.html).

--------------------------------

## Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](https://youtu.be/KPEqYtEd-zY)

## Video Introduction

[![Video Introduction to DifferentialEquations.jl](https://user-images.githubusercontent.com/1814174/27973992-e236a9a4-6310-11e7-84af-2b66097cecf9.PNG)](https://youtu.be/75SCMIRlNXM)

## Comparison with MATLAB, R, Julia, Python, C, Mathematica, Maple, and Fortran

<a href="http://www.stochasticlifestyle.com/wp-content/uploads/2019/08/de_solver_software_comparsion.pdf"><img src="http://www.stochasticlifestyle.com/wp-content/uploads/2019/08/de_solver_software_comparsion-1.png" alt="Comparison Of Differential Equation Solver Software" align="middle"/></a>

[See the corresponding blog post](http://www.stochasticlifestyle.com/comparison-differential-equation-solver-suites-matlab-r-julia-python-c-fortran/)

## Example Images

<img src="https://raw.githubusercontent.com/SciML/DifferentialEquations.jl/master/assets/DifferentialEquations_Example.png" align="middle"  />

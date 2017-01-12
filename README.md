# DifferentialEquations.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Travis](https://travis-ci.org/JuliaDiffEq/DifferentialEquations.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DifferentialEquations.jl)
[![AppVoyer](https://ci.appveyor.com/api/projects/status/1smlr9ryfqfx1ear?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/differentialequations-jl-1sx90)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](http://docs.juliadiffeq.org/latest/)
[![DOI](https://zenodo.org/badge/58516043.svg)](https://zenodo.org/badge/latestdoi/58516043)

This is a suite for numerically solving differential equations in Julia. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations. Equations within the realm of this package include:

- Ordinary differential equations (ODEs)
- Stochastic ordinary differential equations (SODEs or SDEs)
- Algebraic differential equations (DAEs)
- Delay differential equations (DDEs)
- (Stochastic) partial differential equations ((S)PDEs) (with both finite difference and finite element methods)

The well-optimized DifferentialEquations solvers benchmark as the fastest Julia implementations, using classic algorithms and ones from recent research, and include algorithms optimized for high-precision and HPC applications.  It integrates with the Julia package sphere, for example using Juno's progress meter, automatic plotting, built-in interpolations, and wraps other differential equation solvers so that many different methods for solving the equations can be accessed by simply switching a keyword argument. It utilizes Julia's generality to be able to solve problems specified with arbitrary number types (types with units like Unitful, and arbitrary precision numbers like BigFloats and ArbFloats), arbitrary sized arrays (ODEs on matrices), and more. This gives a powerful mixture of speed and productivity features to help you solve and analyze your differential equations faster.

For information on using the package, [see the documentation](http://docs.juliadiffeq.org/latest/).

All of the algorithms are thoroughly tested to ensure accuracy via convergence tests. The algorithms are continuously tested to show correctness. IJulia tutorial notebooks [can be found at DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl). Benchmarks [can be found at DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl). If you find any equation where there seems to be an error, please open an issue.

If you have any questions, or just want to chat about solvers/using the package, please feel free to chat in the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge). For bug reports, feature requests, etc., please submit an issue. If you're interested in contributing, please see the [Developer Documentation](https://juliadiffeq.github.io/DiffEqDevDocs.jl/latest/).

--------------------------------


<img src="https://raw.githubusercontent.com/JuliaDiffEq/DifferentialEquations.jl/master/assets/DifferentialEquations_Example.png" align="middle"  />

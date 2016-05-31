# TODO List

## General
* Change all solvers to work on systems of equations
* Upgrade to v0.5

## Algorithms
* Finite difference solvers:
  * Semi-linear Heat Equation (Reaction-Diffusion Equation)
  * Semi-linear Poisson Equation
  * Wave Equation
  * Transport Equation
  * Implicit Integration Factor (IIF) Maruyama
  * Implicit Integration Factor (IIF) Milstein
* Finite Difference Geometric Multigrids
* Algebraic Multigrids via pyAMG
* Adaptive-SRK

## Performance
* Implement simd/threading in core loops?
* Test ParallelAccelerator.jl on functionalities
* Turn off bounds checking in core loops
* Add an option for fast math
* Add Xeon Phi / GPU variants for PDE solvers

## Misc
* Davie-Gaines convergence analysis
* CG and GMRES don't work with systems, change IterativeSolvers

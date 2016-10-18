# Finite Element Method Heat Equation Solvers

## Avaliable Methods

[method] denotes an additional version for handling stochastic partial
differential equations.

* Finite Element Solvers (Stochastic) PDEs

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

## Solver Documentation

```@docs
solve(::FEMmesh,::HeatProblem)
```

# Stochastic Differential Equation Solvers

## Implemented Solvers

In addition to the standard Euler-Maruyama method, specialized versions of higher
order Runge-Kutta methods are implemented which give increased accuracy and speed.

  * Euler-Maruyama
  * Milstein
  * Rossler-SRK

## Solver Documentation

```@docs
solve(::AbstractSDEProblem,::AbstractArray)
```

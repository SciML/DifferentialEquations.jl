# Breaking Changes in v7.0

- ParameterizedFunctions.jl, along with a few other modeling libraries, are no 
  longer exported by DifferentialEquations.jl. You must do `using ParameterizedFunctions` 
  to use the `@ode_def macro. As it has not been part of the documentation for 3 years,
  this is breaking but should not effect most users. This reduces the dependencies and
  using time of the library by about half.
- OrdinaryDiffEq.jl has a new linear solver system based on the LinearSolve.jl library.
  `linsolve` arguments now take LinearSolve.jl solvers.

#__precompile__()

module DifferentialEquations

  using Reexport

  @reexport using DiffEqBase
  @reexport using StochasticDiffEq
  @reexport using FiniteElementDiffEq
  @reexport using DiffEqDevTools
  @reexport using OrdinaryDiffEq
  @reexport using AlgebraicDiffEq
  @reexport using StokesDiffEq
  @reexport using DiffEqParamEstim
  @reexport using Plots

end # module

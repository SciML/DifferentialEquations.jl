__precompile__()

module DifferentialEquations

  using Reexport

  @reexport using DiffEqBase
  @reexport using DiffEqNoiseProcess
  @reexport using RecursiveArrayTools

  @reexport using SteadyStateDiffEq
  @reexport using StochasticDiffEq
  @reexport using OrdinaryDiffEq
  @reexport using BoundaryValueDiffEq
  using Sundials
  @reexport using DelayDiffEq

  @reexport using DiffEqCallbacks
  @reexport using DiffEqMonteCarlo
  @reexport using DiffEqJump

  @reexport using DiffEqFinancial
  @reexport using MultiScaleArrays
  @reexport using DiffEqPhysics

  @reexport using DimensionalPlotRecipes

  using LinearAlgebra

  import DiffEqBase: solve

  include("default_solve.jl")
  include("default_arg_parsing.jl")
  include("ode_default_alg.jl")
  include("sde_default_alg.jl")
  include("dae_default_alg.jl")
  include("dde_default_alg.jl")
  include("discrete_default_alg.jl")
  include("rode_default_alg.jl")
  include("steady_state_default_alg.jl")

end # module

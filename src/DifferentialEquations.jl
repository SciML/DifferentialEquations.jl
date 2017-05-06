__precompile__(false)

module DifferentialEquations

  using Reexport

  @reexport using DiffEqBase
  @reexport using DiffEqPDEBase
  @reexport using DiffEqNoiseProcess
  @reexport using RecursiveArrayTools

  @reexport using SteadyStateDiffEq
  @reexport using StochasticDiffEq
  @reexport using FiniteElementDiffEq
  @reexport using OrdinaryDiffEq
  @reexport using AlgebraicDiffEq
  @reexport using Sundials
  @reexport using DelayDiffEq

  @reexport using ParameterizedFunctions
  @reexport using DiffEqParamEstim
  @reexport using DiffEqSensitivity
  @reexport using DiffEqCallbacks
  @reexport using DiffEqMonteCarlo
  @reexport using DiffEqDevTools
  @reexport using DiffEqJump
  @reexport using DiffEqUncertainty

  @reexport using DiffEqFinancial
  @reexport using DiffEqBiological
  @reexport using MultiScaleArrays

  @reexport using PyDSTool

  @reexport using DimensionalPlotRecipes

  import DiffEqBase: solve

  include("default_solve.jl")
  include("default_arg_parsing.jl")
  include("ode_default_alg.jl")
  include("sde_default_alg.jl")
  include("dae_default_alg.jl")
  include("dde_default_alg.jl")
  include("fem_default_alg.jl")


end # module

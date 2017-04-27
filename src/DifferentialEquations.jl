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
  @reexport using StokesDiffEq
  @reexport using Sundials
  @reexport using DelayDiffEq

  @reexport using ParameterizedFunctions
  @reexport using DiffEqParamEstim
  @reexport using DiffEqSensitivity
  @reexport using DiffEqCallbacks
  @reexport using DiffEqMonteCarlo
  @reexport using DiffEqDevTools
  @reexport using DiffEqJump

  @reexport using DiffEqFinancial
  @reexport using DiffEqBiological
  @reexport using MultiScaleArrays

  @reexport using PyDSTool

  @reexport using DimensionalPlotRecipes

  import DiffEqBase: solve

  include("ode_default_alg.jl")
  include("sde_default_alg.jl")
  include("dae_default_alg.jl")
  include("dde_default_alg.jl")
  include("fem_default_alg.jl")

  function solve(prob::DEProblem,args...;default_set=false,kwargs...)
    if default_set == true && !isempty(args)
      error("The chosen algorithm, "*string(args[1])*", does not exist.
        Please verify that the appropriate solver package has been installed.")
    elseif default_set == true
      error("No algorithm is chosen but default_set=true.")
    end
    alg,extra_kwargs = default_algorithm(prob;kwargs...)
    solve(prob,alg,args...;default_set=true,kwargs...,extra_kwargs...)
  end

  export default_algorithm


end # module

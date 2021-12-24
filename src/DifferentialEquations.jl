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
  @reexport using DiffEqJump

  using LinearAlgebra

  import DiffEqBase: solve
  import LinearSolve

  include("default_solve.jl")
  include("default_arg_parsing.jl")
  include("ode_default_alg.jl")
  include("sde_default_alg.jl")
  include("dae_default_alg.jl")
  include("dde_default_alg.jl")
  include("discrete_default_alg.jl")
  include("rode_default_alg.jl")
  include("steady_state_default_alg.jl")
  include("bvp_default_alg.jl")

end # module

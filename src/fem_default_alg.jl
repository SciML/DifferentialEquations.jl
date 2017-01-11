function default_algorithm{islinear,isstochastic,MeshType<:FEMMesh}(
          prob::AbstractPoissonProblem{islinear,isstochastic,MeshType};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=FEMDiffEqPoisson() # Standard default
  uEltype = eltype(prob.u0)

  alg, extra_kwargs
end

function default_algorithm{islinear,isstochastic,MeshType<:FEMMesh}(
          prob::AbstractHeatProblem{islinear,isstochastic,MeshType};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=FEMDiffEqHeatEuler() # Standard default
  uEltype = eltype(prob.u0)

  :alg_hints ∈ keys(o) ? alg_hints = o[:alg_hints] : alg_hints = Symbol[:nonstiff]

  if :stiff ∈ alg_hints && :nonstiff ∈ alg_hints
    error("The problem must either be designated as stiff or non-stiff")
  end

  if :stiff ∈ alg_hints
    alg = FEMDiffEqHeatSemiImplicitCrankNicholson()
  end

  alg, extra_kwargs
end

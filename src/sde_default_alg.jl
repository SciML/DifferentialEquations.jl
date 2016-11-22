function default_algorithm{uType,tType,isinplace,NoiseClass,F,F2,F3}(prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=SRIW1() # Standard default
  uEltype = eltype(prob.u0)

  :alg_hints ∈ keys(o) ? alg_hints = o[:alg_hints] : alg_hints = Symbol[:nonstiff]

  if :additive ∈ alg_hints
    alg = SRA1()
  end


  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)

  alg,extra_kwargs
end

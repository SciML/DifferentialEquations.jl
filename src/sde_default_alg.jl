function default_algorithm{uType,tType,isinplace,ND}(prob::AbstractSDEProblem{uType,tType,isinplace,ND};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=SRIW1() # Standard default
  uEltype = eltype(prob.u0)

  alg_hints = get_alg_hints(o)

  if :additive ∈ alg_hints
    alg = SRA1()
  end

  if prob.noise_rate_prototype != nothing
    alg = EM()
  end

  if :stratonovich ∈ alg_hints
    alg = EulerHeun()
  end

  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)

  alg,extra_kwargs
end

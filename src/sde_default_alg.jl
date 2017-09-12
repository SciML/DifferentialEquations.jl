function default_algorithm{uType,tType,isinplace,ND}(prob::AbstractSDEProblem{uType,tType,isinplace,ND};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=SRIW1() # Standard default
  uEltype = eltype(prob.u0)

  alg_hints = get_alg_hints(o)

  if :additive ∈ alg_hints
    alg = SRA1()
  end

  if :commutative ∈ alg_hints
    alg = RKMilCommute()
  end

  if :stiff ∈ alg_hints
    alg = ImplicitRKMil()
  end

  if :stratonovich ∈ alg_hints
    if :stiff ∈ alg_hints
      alg = ImplicitRKMil(interpretation=:stratonovich)
    else
      alg = RKMil(interpretation=:stratonovich)
    end
  end

  if prob.noise_rate_prototype != nothing || prob.noise != nothing
    if :stratonovich ∈ alg_hints
      if :stiff ∈ alg_hints
        alg = ImplicitEulerHeun()
      else
        alg = EulerHeun()
      end
    else
      if :stiff ∈ alg_hints
        alg = ImplicitEM()
      else
        alg = EM()
      end
    end
  end

  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)

  alg,extra_kwargs
end

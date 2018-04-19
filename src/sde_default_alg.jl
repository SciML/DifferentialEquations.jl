function default_algorithm{uType,tType,isinplace,ND}(prob::AbstractSDEProblem{uType,tType,isinplace,ND};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=SOSRI() # Standard default
  uEltype = eltype(prob.u0)

  alg_hints = get_alg_hints(o)

  if :commutative ∈ alg_hints
    alg = RKMilCommute()
  end

  if :stiff ∈ alg_hints
    alg = ImplicitRKMil(autodiff=false)
  end

  if :stratonovich ∈ alg_hints
    if :stiff ∈ alg_hints
      alg = ImplicitRKMil(autodiff=false,interpretation=:stratonovich)
    else
      alg = RKMil(interpretation=:stratonovich)
    end
  end

  if prob.noise_rate_prototype != nothing || prob.noise != nothing
    if :stratonovich ∈ alg_hints
      if :stiff ∈ alg_hints
        alg = ImplicitEulerHeun(autodiff=false)
      else
        alg = LambaEulerHeun()
      end
    else
      if :stiff ∈ alg_hints
        alg = ISSEM(autodiff=false)
      else
        alg = LambaEM()
      end
    end
  end

  if :additive ∈ alg_hints
    if :stiff ∈ alg_hints
      alg = SKenCarp(autodiff=false)
    else
      alg = SOSRA()
    end
  end

  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)

  alg,extra_kwargs
end

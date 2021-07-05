function default_algorithm(prob::DiffEqBase.AbstractSDEProblem{uType,tType,isinplace,ND};kwargs...) where {uType,tType,isinplace,ND}
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=SOSRI() # Standard default
  uEltype = eltype(prob.u0)

  alg_hints = get_alg_hints(o)

  if :commutative ∈ alg_hints
    alg = RKMilCommute()
  end

  is_stiff = :stiff ∈ alg_hints
  is_stratonovich = :stratonovich ∈ alg_hints
  if is_stiff || prob.f.mass_matrix !== I
    alg = ImplicitRKMil(autodiff=false)
  end

  if is_stratonovich
    if is_stiff || prob.f.mass_matrix !== I
      alg = ImplicitRKMil(autodiff=false,interpretation=:stratonovich)
    else
      alg = RKMil(interpretation=:stratonovich)
    end
  end

  if prob.noise_rate_prototype != nothing || prob.noise != nothing
    if is_stratonovich
      if is_stiff || prob.f.mass_matrix !== I
        alg = ImplicitEulerHeun(autodiff=false)
      else
        alg = LambaEulerHeun()
      end
    else
      if is_stiff || prob.f.mass_matrix !== I
        alg = ISSEM(autodiff=false)
      else
        alg = LambaEM()
      end
    end
  end

  if :additive ∈ alg_hints
    if is_stiff || prob.f.mass_matrix !== I
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

function default_algorithm{uType,tType,inplace}(prob::AbstractODEProblem{uType,tType,inplace};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=Tsit5 # Standard default
  uEltype = eltype(prob.u0)

  :alg_hints ∈ keys(o) ? alg_hints = o[:alg_hints] : alg_hints = Symbol[:nonstiff]

  if :stiff ∈ alg_hints && :nonstiff ∈ alg_hints
    error("The problem must either be designated as stiff or non-stiff")
  end

  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)


  if :nonstiff ∈ alg_hints
    if uEltype <: AbstractFloat
      if !(uEltype <: Float64)
        # Most likely higher precision, so use a higher order method
        alg = Vern8
      end
    end
  else # The problem is stiff
    if uEltype <: Float64 # Sundials only works on Float64!
      alg = CVODE_BDF
    else
      alg = Rosenbrock23
    end
  end
  alg,extra_kwargs
end

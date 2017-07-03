function default_algorithm{uType,tType,inplace}(prob::AbstractODEProblem{uType,tType,inplace};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=Tsit5() # Standard default
  uEltype = eltype(prob.u0)

  alg_hints = get_alg_hints(o)
  tol_level = get_tolerance_level(o)
  callbacks = callbacks_exists(o)
  mm = mass_matrix_exists(prob)

  if :stiff ∈ alg_hints && :nonstiff ∈ alg_hints
    error("The problem must either be designated as stiff or non-stiff")
  end

  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)


  if :nonstiff ∈ alg_hints
    if uEltype <: AbstractFloat
      if !(uEltype <: Float64) || tol_level == :extreme_tol
        # Most likely higher precision, so use a higher order method
        alg = Vern9()
      elseif tol_level == :low_tol
        alg = Vern7()
      elseif tol_level == :med_tol
        alg = Tsit5()
      else # tol_level == :high_tol
        alg = BS3()
      end
    end
  else # The problem is stiff
    if uEltype <: Float64 # Sundials only works on Float64!
      if tol_level == :high_tol || callbacks || mm # But not in these cases
        alg = Rosenbrock23()
      else
        alg = CVODE_BDF()
      end
    else
      alg = Rosenbrock23()
    end
  end
  alg,extra_kwargs
end

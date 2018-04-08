function default_algorithm{uType,tType,inplace}(prob::AbstractODEProblem{uType,tType,inplace};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=AutoTsit5(Rosenbrock23()) # Standard default
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
  elseif :stiff ∈ alg_hints # The problem is stiff
    if uType <: Array{Float64} && !mm && length(prob.u0) > 10000
      # Use Krylov method when huge!
      alg = CVODE_BDF(linear_solver=:BCG)
    elseif uType <: Array{Float64} && !mm && length(prob.u0) > 1000
      # Sundials only works on Float64!
      # Sundials is fast when problems are large enough
      alg = CVODE_BDF()
    else
      alg = Rosenbrock23()
    end
    if tol_level == :high_tol
      alg = Rosenbrock23()
    else
      if eltype(prob.u0) <: Float64
        alg = Rodas4()
      else # This is only for Julia v0.6 lufact! bug
        alg = Rodas4(linsolve=LinSolveFactorize(qrfact!))
      end
    end
  else # :auto ∈ alg_hints
    if uEltype <: AbstractFloat
      if !(uEltype <: Float64) || tol_level == :extreme_tol
        # Most likely higher precision, so use a higher order method
        alg = AutoVern9(Rodas5())
      elseif tol_level == :low_tol
        alg = AutoVern7(Rodas4())
      else # :med or low
        alg = AutoTsit5(Rosenbrock23())
      end
    end
  end
  alg,extra_kwargs
end

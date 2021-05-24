function default_algorithm(prob::DiffEqBase.AbstractODEProblem{uType,tType,inplace};kwargs...) where {uType,tType,inplace}
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg=AutoTsit5(Rosenbrock23(autodiff=false)) # Standard default
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

  if typeof(prob.f) <: SplitFunction
    alg = KenCarp4(autodiff=false)
  elseif typeof(prob.f) <: DynamicalODEFunction
    if tol_level == :low_tol || tol_level == :med_tol
      alg = DPRKN6()
    elseif tol_level == :high_tol
      alg = DPRKN8()
    else
      alg = DPRKN12()
    end
  else # Standard ODE
    if :nonstiff ∈ alg_hints || length(prob.u0) > 10000
      # Don't default to implicit here because of memory requirements
      # And because the linear system gets unruly
      if (!(uEltype <: Float64) && !(uEltype <: Float32) && !(uEltype <: Complex)) || tol_level == :extreme_tol
        # Most likely higher precision, so use a higher order method
        alg = Vern9(lazy=!callbacks)
      elseif tol_level == :low_tol
        alg = Vern7(lazy=!callbacks)
      elseif tol_level == :med_tol
        alg = Tsit5()
      else # tol_level == :high_tol
        alg = BS3()
      end
    elseif :stiff ∈ alg_hints # The problem is stiff
      if uType <: Array{Float64} && !mm && length(prob.u0) > 1000
        # Use Krylov method when huge!
        alg = CVODE_BDF(linear_solver=:GMRES)
      elseif length(prob.u0) > 200
        # Sundials only works on Float64!
        # Sundials is fast when problems are large enough
        alg = QNDF()
      elseif tol_level == :high_tol
        alg = Rosenbrock23(autodiff=false)
      else
        alg = Rodas4(autodiff=false)
      end
    else # :auto ∈ alg_hints
      if (!(uEltype <: Float64) && !(uEltype <: Float32)) || tol_level == :extreme_tol
        # Most likely higher precision, so use a higher order method
        alg = AutoVern9(Rodas5(autodiff=false),lazy=!callbacks)
      elseif tol_level == :low_tol
        alg = AutoVern7(Rodas4(autodiff=false),lazy=!callbacks)
      else # :med or low
        alg = AutoTsit5(Rosenbrock23(autodiff=false))
      end
    end
  end
  alg,extra_kwargs
end

function default_algorithm{uType,tType,lType,isinplace}(prob::AbstractDDEProblem{uType,tType,lType,isinplace};kwargs...)
  o = Dict{Symbol,Any}(kwargs)
  extra_kwargs = Any[]; alg = MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false))) # Standard default
  uEltype = eltype(prob.u0)

  alg_hints = get_alg_hints(o)
  tol_level = get_tolerance_level(o)

  if tol_level == :extreme_tol || tol_level == :low_tol
    stiff_alg = Rodas5(autodiff=false)
  else
    stiff_alg = Rosenbrock23(autodiff=false)
  end

  if :stiff ∈ alg_hints
    alg = MethodOfSteps(stiff_alg)
  elseif :nonstiff ∈ alg_hints || length(prob.u0) > 10000
    alg = MethodOfSteps(Tsit5())
  else # :auto
    if tol_level == :extreme_tol || tol_level == :low_tol
      alg = MethodOfSteps(AutoVern7(stiff_alg))
    else
      alg = MethodOfSteps(AutoTsit5(stiff_alg))
    end
  end

  # If adaptivity is not set and the tType is not a float, turn off adaptivity
  # Bad interaction with ForwardDiff
  #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)

  alg,extra_kwargs
end

function default_algorithm(prob::DiffEqBase.AbstractODEProblem{uType, tType, inplace};
    kwargs...) where {uType, tType, inplace}
    o = Dict{Symbol, Any}(kwargs)
    extra_kwargs = Any[]
    alg = AutoTsit5(Rosenbrock23(autodiff = false)) # Standard default
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
        alg = KenCarp4(autodiff = false)
    elseif typeof(prob.f) <: DynamicalODEFunction
        if tol_level == :low_tol || tol_level == :med_tol
            alg = Tsit5()
        else
            alg = Vern7(lazy = !callbacks)
        end
    else # Standard ODE
        if :nonstiff ∈ alg_hints
            # Don't default to implicit here because of memory requirements
            # And because the linear system gets unruly
            if (!(uEltype <: Float64) && !(uEltype <: Float32) && !(uEltype <: Complex)) ||
               tol_level == :extreme_tol || tol_level == :low_tol
                # Most likely higher precision, so use a higher order method
                alg = Vern7(lazy = !callbacks)
            else
                alg = Tsit5()
            end
        elseif :stiff ∈ alg_hints || mm # The problem is stiff
            if length(prob.u0) > 500
                # Use Krylov method when huge!
                alg = FBDF(autodiff = false, linsolve = LinearSolve.KrylovJL_GMRES())
            elseif length(prob.u0) > 500
                alg = FBDF(autodiff = false)
            elseif tol_level == :high_tol
                alg = Rosenbrock23(autodiff = false)
            else
                alg = Rodas5P(autodiff = false)
            end
        else # :auto ∈ alg_hints
            if (!(uEltype <: Float64) && !(uEltype <: Float32)) || tol_level == :extreme_tol  ||
                tol_level == :low_tol
                # Most likely higher precision, so use a higher order method
                if length(prob.u0) > 500
                    alg = AutoVern7(KenCarp47(autodiff = false, linsolve = LinearSolve.KrylovJL_GMRES()), lazy = !callbacks)
                elseif length(prob.u0) > 50
                    alg = AutoVern7(KenCarp47(autodiff = false), lazy = !callbacks)
                else
                    alg = AutoVern7(Rodas5P(autodiff = false), lazy = !callbacks)
                end
            else # :med or low
                if length(prob.u0) > 500
                    alg = AutoTsit5(TRBDF2(autodiff = false, linsolve = LinearSolve.KrylovJL_GMRES()))
                elseif length(prob.u0) > 50
                    alg = AutoTsit5(TRBDF2(autodiff = false))
                else
                    alg = AutoTsit5(Rosenbrock23(autodiff = false))
                end
            end
        end
    end
    alg, extra_kwargs
end

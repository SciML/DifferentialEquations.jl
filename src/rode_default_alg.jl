function default_algorithm(prob::DiffEqBase.AbstractRODEProblem; kwargs...)
    o = Dict{Symbol, Any}(kwargs)
    extra_kwargs = Any[]
    alg = RandomEM() # Standard default

    alg_hints = get_alg_hints(o)

    # If adaptivity is not set and the tType is not a float, turn off adaptivity
    # Bad interaction with ForwardDiff
    #!(tType <: AbstractFloat) && (:adaptive ∉ keys(o)) && push!(extra_kwargs,:adaptive=>false)

    alg, extra_kwargs
end

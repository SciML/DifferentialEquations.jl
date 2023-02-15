function DiffEqBase.__solve(prob::DiffEqBase.DEProblem,
                            alg::Union{Nothing, DiffEqBase.DEAlgorithm},
                            args...; default_set = false, kwargs...)
    if default_set == true
        error("The chosen algorithm, $alg, does not exist. Please verify that the appropriate solver package has been installed.")
    end
    alg, extra_kwargs = default_algorithm(prob; kwargs...)
    alg = DiffEqBase.prepare_alg(alg, prob.u0, prob.p, prob)

    if alg isa SciMLBase.AbstractNonlinearAlgorithm && prob isa SciMLBase.AbstractSteadyStateProblem
      prob = NonlinearProblem(prob)
    end
    DiffEqBase.__solve(prob, alg, args...; default_set = true, kwargs..., extra_kwargs...)
end

function DiffEqBase.__init(prob::DiffEqBase.DEProblem,
                           alg::Union{Nothing, DiffEqBase.DEAlgorithm},
                           args...; default_set = false, kwargs...)
    if default_set == true
        error("The chosen algorithm, $alg, does not exist. Please verify that the appropriate solver package has been installed.")
    end
    alg, extra_kwargs = default_algorithm(prob; kwargs...)
    alg = DiffEqBase.prepare_alg(alg, prob.u0, prob.p, prob)
    DiffEqBase.__init(prob, alg, args...; default_set = true, kwargs..., extra_kwargs...)
end

#=
function DiffEqBase.__solve(prob::JumpProcesses.JumpProblem,alg::Nothing,
               args...;default_set=false,kwargs...)
  if default_set == true
    error("The chosen algorithm, $alg, does not exist. Please verify that the appropriate solver package has been installed.")
  end
  alg,extra_kwargs = default_algorithm(prob.prob;kwargs...)
  DiffEqBase.__solve(prob,alg,args...;default_set=true,kwargs...,extra_kwargs...)
end
=#

export default_algorithm, DefaultAlgorithm

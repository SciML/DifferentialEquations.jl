immutable DefaultAlgorithm <: DEAlgorithm end

function solve(prob::DEProblem,args...;default_set=false,kwargs...)
  if length(args)>0 && !(typeof(args[1]) <: DEAlgorithm)
    error("Inappropiate solve command. The arguments do not make sense. Likely, you gave an algorithm which does not actually exist (or does not <:DEAlgorithm)")
  else
    solve(prob::DEProblem,nothing,args...;default_set=false,kwargs...)
  end
end

function solve(prob::DEProblem,alg::DEAlgorithm,args...;default_set=false,kwargs...)
  error("The chosen algorithm, $alg, does not exist. Please verify that the appropriate solver package has been installed.")
end

function solve(prob::DEProblem,::Union{Void,DefaultAlgorithm},args...;default_set=false,kwargs...)
  if default_set == true
    error("No algorithm is chosen but default_set=true.")
  end
  alg,extra_kwargs = default_algorithm(prob;kwargs...)
  solve(prob,alg,args...;default_set=true,kwargs...,extra_kwargs...)
end

export default_algorithm, DefaultAlgorithm

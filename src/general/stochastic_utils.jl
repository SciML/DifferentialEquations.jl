"""
`getNoise(N,node,elem;noisetype=:White)`

Returns a random vector corresponding to the noise type which was chosen.
"""
function getNoise(u,node,elem;noisetype=:White)
  if noisetype==:White
    return(ChunkedArray(randn,u))
  end
end

"""
`monteCarloSim(Δt::Number,prob::SDEProblem)`

Performs a parallel Monte-Carlo simulation to solve the SDE problem with Δt numMonte times.
Returns a vector of solution objects.

### Keyword Arguments
* `T` - Final time. Default is 1.
* `numMonte` - Number of Monte-Carlo simulations to run. Default is 10000
* `save_timeseries` - Denotes whether save_timeseries should be turned on in each run. Default is false.
"""
function monteCarloSim(prob::AbstractSDEProblem,tspan=[0,1];numMonte=10000,save_timeseries=false,kwargs...)
  elapsedTime = @elapsed solutions = pmap((i)->solve(prob,tspan;save_timeseries=save_timeseries,kwargs...),1:numMonte)
  solutions = convert(Array{SDESolution},solutions)
  if prob.knownanalytic
    N = size(solutions,1)
    errors = Dict() #Should add type information
    error_means  = Dict()
    error_medians= Dict()
    for k in keys(solutions[1].errors)
      errors[k] = reshape(Float64[sol.errors[k] for sol in solutions],size(solutions)...)
      error_means[k] = mean(errors[k])
      error_medians[k]=median(errors[k])
    end
  end
  return(MonteCarloSimulation(solutions,errors,error_means,error_medians,N,elapsedTime))
end

type MonteCarloSimulation
  solutions::Array{SDESolution}
  errors
  error_means
  error_medians
  N
  elapsedTime
end

Base.length(sim::MonteCarloSimulation) = sim.N
Base.endof( sim::MonteCarloSimulation) = length(sim)
Base.getindex(sim::MonteCarloSimulation,i::Int) = sim.solutions[i]
Base.getindex(sim::MonteCarloSimulation,i::Int,I::Int...) = sim.solutions[i][I]

function print(io::IO, sim::MonteCarloSimulation)
  println(io,"$(typeof(sim)) of length $(length(sim)).")
  println(io,"\n-----------Errors-----------")
  for (k,v) in sim.errors
    println(io,"$k: $v")
  end
end

function show(io::IO,sim::MonteCarloSimulation)
  println(io,"$(typeof(sim)) of length $(length(sim)).")
end

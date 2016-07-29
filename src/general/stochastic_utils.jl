"""
getNoise(N,node,elem;noisetype=:White)

Returns a random vector corresponding to the noise type which was chosen.
"""
function getNoise(u,node,elem;noisetype=:White)
  if noisetype==:White
    return(ChunkedArray(randn,u))
  end
end

"""
monteCarloSim(Δt::Number,prob::SDEProblem)

Performs a parallel Monte-Carlo simulation to solve the SDE problem with Δt numMonte times.
Returns a vector of solution objects.

### Keyword Arguments
* T - Final time. Default is 1.
* numMonte - Number of Monte-Carlo simulations to run. Default is 10000
* save_timeseries - Denotes whether save_timeseries should be turned on in each run. Default is true.
* alg - Algorithm for solving the SDEs. Default is "EM"
"""
function monteCarloSim(prob::SDEProblem;Δt::Number=0,tspan=[0,1],numMonte=10000,save_timeseries=false,alg=:SRIW1Optimized,adaptive=false,abstol=1e-3,reltol=1e-2,adaptivealg=:RSwM3,qmax=4)
  elapsedTime = @elapsed solutions = pmap((i)->solve(prob,tspan,Δt=Δt,save_timeseries=save_timeseries,alg=alg,adaptive=adaptive,abstol=abstol,reltol=reltol,adaptivealg=adaptivealg,qmax=qmax),1:numMonte)
  solutions = convert(Array{SDESolution},solutions)
  if prob.knownanalytic
    N = size(solutions,1)
    errors = Dict() #Should add type information
    means  = Dict()
    medians= Dict()
    for k in keys(solutions[1].errors)
      errors[k] = reshape(Float64[sol.errors[k] for sol in solutions],size(solutions)...)
      means[k] = mean(errors[k])
      medians[k]=median(errors[k])
    end
  end
  return(MonteCarloSimulation(solutions,errors,means,medians,N,elapsedTime))
end

type MonteCarloSimulation
  solutions::Array{SDESolution}
  errors
  means
  medians
  N
  elapsedTime
end

Base.length(sim::MonteCarloSimulation) = sim.N

"""
getNoise(N,node,elem;noiseType="White")

Returns a random vector corresponding to the noise type which was chosen.
"""
function getNoise(u,node,elem;noiseType="White")
  if noiseType=="White"
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
* fullSave - Denotes whether fullSave should be turned on in each run. Default is true.
* alg - Algorithm for solving the SDEs. Default is "EM"
"""
function monteCarloSim(prob::SDEProblem;Δt::Number=0,tspan=[0,1],numMonte=10000,fullSave=true,alg="EM")
  solutions = pmap((i)->solve(prob,tspan,Δt=Δt,fullSave=fullSave,alg=alg),1:numMonte)
  solutions = convert(Array{SDESolution},solutions)
  return(solutions)
end

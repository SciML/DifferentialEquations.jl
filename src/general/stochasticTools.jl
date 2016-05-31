"""
getNoise(N,node,elem;noiseType="White")

Returns a random vector corresponding to the noise type which was chosen.
"""
function getNoise(u,node,elem;noiseType="White")
  if noiseType=="White"
    return(ChunkedArray(randn,u))
  end
end

function monteCarloSim(Δt::Number,prob::SDEProblem;T=1,numMonte=10000,fullSave=true,alg="EM")
  solutions = pmap((i)->solve(prob,Δt,T,fullSave=fullSave,alg=alg),1:numMonte)
  solutions = convert(Array{SDESolution},solutions)
  return(solutions)
end

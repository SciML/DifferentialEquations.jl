using DifferentialEquations
prob = linearSDEExample(α=1/10,β=1/20,u₀=1/2)

## Solve and plot

msims = Vector{MonteCarloSimulation}(7)
@progress for i=0:6
  msims[i+1] = monteCarloSim(prob::SDEProblem,Δt=1//2^(4),adaptive=true,numMonte=1000,abstol=10.0^(-i),reltol=0)
end


prob = waveSDEExample()

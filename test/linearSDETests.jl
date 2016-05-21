using DifferentialEquations
prob = linearSDEExample()

## Solve and plot
sol =solve(prob::SDEProblem,1//2^(4),1,fullSave=true,alg="SRI")

PyPlot.plot(sol.tFull,sol.uFull)
PyPlot.plot(sol.tFull,sol.solFull)

## Convergence Testing
Δts = 1.//2.^(10:-1:4) #14->7 good plot

convsim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")

convsim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")

convsim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")

abs(convsim.𝒪est["l2"]-.5) + abs(convsim2.𝒪est["l∞"]-1) + abs(convsim3.𝒪est["final"]-1.5)<.5 #High tolerance since low Δts for testing!

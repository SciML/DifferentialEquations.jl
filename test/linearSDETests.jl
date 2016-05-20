using DifferentialEquations

prob = linearSDEExample()

sol =solve(prob::SDEProblem,1//2^(5),1,fullSave=true)

Δts = 1.//2.^(12:-1:5)

convsim = testConvergence(Δts,prob,numMonte=Int(1e2))

convsim.𝒪est["final"]-.5 < 0.1

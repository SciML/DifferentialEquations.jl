@everywhere using DifferentialEquations
srand(100)
Δts = 1.//2.^(10:-1:4) #14->7 good plot

prob = waveSDEExample()
sim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")
sim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")
sim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")
bool1 = abs(sim.𝒪est["l2"]-.5) + abs(sim2.𝒪est["l∞"]-1) + abs(sim3.𝒪est["final"]-1.5) <.5 #High tolerance since low Δts for testing!

prob = cubicSDEExample()
sim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")
sim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")
sim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")
bool2 = abs(sim.𝒪est["l2"]-.5) + abs(sim2.𝒪est["l∞"]-1) + abs(sim3.𝒪est["final"]-1.5) <.5 #High tolerance since low Δts for testing!

## Convergence Testing
prob = additiveSDEExample()
sim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")
sim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")
sim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")
sim4 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRA")
bool3 = abs(sim.𝒪est["l2"]-1) + abs(sim2.𝒪est["l∞"]-1) + abs(sim3.𝒪est["final"]-2) + abs(sim4.𝒪est["final"]-2) <.4 #High tolerance since low Δts for testing!

bool1 && bool2 && bool3

using DifferentialEquations
srand(100)
Δts = 1.//2.^(10:-1:4) #14->7 good plot

prob = waveSDEExample()
convsim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")
convsim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")
convsim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")
bool1 = abs(convsim.𝒪est["l2"]-.5) + abs(convsim2.𝒪est["l∞"]-1) + abs(convsim3.𝒪est["final"]-1.5) <.5 #High tolerance since low Δts for testing!

prob = cubicSDEExample()
convsim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")
convsim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")
convsim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")
bool2 = abs(convsim.𝒪est["l2"]-.5) + abs(convsim2.𝒪est["l∞"]-1) + abs(convsim3.𝒪est["final"]-1.5) <.5 #High tolerance since low Δts for testing!

## Convergence Testing
prob = additiveSDEExample()
convsim = testConvergence(Δts,prob,numMonte=Int(1e1),alg="EM")
convsim2 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="RKMil")
convsim3 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRI")
convsim4 = testConvergence(Δts,prob,numMonte=Int(1e1),alg="SRA")
bool3 = abs(convsim.𝒪est["l2"]-1) + abs(convsim2.𝒪est["l∞"]-1) + abs(convsim3.𝒪est["final"]-2) + abs(convsim4.𝒪est["final"]-2) <.2 #High tolerance since low Δts for testing!

bool1 && bool2 && bool3

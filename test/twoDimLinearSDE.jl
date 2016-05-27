using DifferentialEquations
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,1//2^(4),1,fullSave=true,alg="RKMil")

fig = PyPlot.figure("pyplot_appx_vs_true",figsize=(10,10))
PyPlot.plot(sol.tFull,sol.uFull[:,1])
PyPlot.plot(sol.tFull,sol.solFull[:,1])
PyPlot.plot(sol.tFull,sol.uFull[:,2])
PyPlot.plot(sol.tFull,sol.solFull[:,2])

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(10:-1:4) #14->7 good plot

convsim = testConvergence(Δts,prob,numMonte=Int(5e1),alg="EM")

convsim2 = testConvergence(Δts,prob,numMonte=Int(5e1),alg="RKMil")

convsim3 = testConvergence(Δts,prob,numMonte=Int(5e1),alg="SRI")

abs(convsim.𝒪est["l2"]-.5) + abs(convsim2.𝒪est["l∞"]-1) + abs(convsim3.𝒪est["final"]-1.5)<.2 #High tolerance since low Δts for testing!

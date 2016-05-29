using DifferentialEquations, Plots
srand(100)
prob = linearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,1//2^(4),1,fullSave=true,alg="SRI")

plot(sol)

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(11:-1:4) #14->7 good plot with higher num Monte

sim = testConvergence(Δts,prob,numMonte=Int(1e2),alg="EM")

sim2 = testConvergence(Δts,prob,numMonte=Int(1e2),alg="RKMil")

sim3 = testConvergence(Δts,prob,numMonte=Int(1e2),alg="SRI")

plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,600))

abs(sim.𝒪est["l2"]-.5) + abs(sim2.𝒪est["l∞"]-1) + abs(sim3.𝒪est["final"]-1.5)<.2 #High tolerance since low Δts for testing!

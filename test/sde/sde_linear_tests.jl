using DifferentialEquations, Plots
srand(100)
prob = linearSDEExample()

## Solve and plot
println("Solve and Plot")
#@code_warntype solve(prob::SDEProblem,[0,1],Δt=1//2^(4),alg=:SRI)

TEST_PLOT && plot(sol)
NUM_MONTE = 10
## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(9:-1:4) #14->7 good plot with higher num Monte

sim = test_convergence(Δts,prob,numMonte=NUM_MONTE,alg=:EM)

sim2 = test_convergence(Δts,prob,numMonte=NUM_MONTE,alg=:RKMil)

sim3 = test_convergence(Δts,prob,numMonte=NUM_MONTE,alg=:SRI)

#TEST_PLOT && plot(plot(sim),plot(sim2),plot(sim3),layout=@layout([a b c]),size=(1200,600))

abs(sim.𝒪est[:l2]-.5) + abs(sim2.𝒪est[:l∞]-1) + abs(sim3.𝒪est[:final]-1.5)<.441  #High tolerance since low Δts for testing!

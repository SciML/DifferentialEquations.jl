using DifferentialEquations
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,Δt=1//2^(4),save_timeseries=true,alg=:SRI)
sol =solve(prob::SDEProblem,[0,1],Δt=0,save_timeseries=true,alg=:SRIW1Optimized)

#Now do the simulation 10000 times in parallel. Return an array
solArr = monteCarloSim(prob::SDEProblem,Δt=1//2^(4))

#First index is the sime, so sol.timeseries[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.timeseries[..,1] returns the time series for the
#first row, and sol.timeseries[..,2] returns the time series for the second.
TEST_PLOT && plot(sol,plottrue=true)

## Convergence Testing
println("Convergence Test on 2D Linear")
Δts = 1.//2.^(8:-1:4) #14->7 good plot

sim = test_convergence(Δts,prob,numMonte=5,alg=:EM)

sim2 = test_convergence(Δts,prob,numMonte=5,alg=:RKMil)

sim3 = test_convergence(Δts,prob,numMonte=5,alg=:SRI)

abs(sim.𝒪est[:l2]-.5) + abs(sim2.𝒪est[:l∞]-1) + abs(sim3.𝒪est[:final]-1.5)<.6 #High tolerance since low Δts for testing!

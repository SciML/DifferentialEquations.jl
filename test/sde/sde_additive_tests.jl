using DifferentialEquations
srand(100)

prob = prob_sde_additive
sol =solve(prob::SDEProblem,[0,1],Δt=1/2^(3),save_timeseries=true,alg=:SRA)
sol =solve(prob::SDEProblem,[0,1],Δt=1/2^(3),save_timeseries=true,alg=:SRA1Optimized)

prob = prob_sde_additivesystem

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,[0,1],Δt=1/2^(3),save_timeseries=true,alg=:SRA)
sol =solve(prob::SDEProblem,[0,1],Δt=1/2^(3),save_timeseries=true,alg=:SRA1Optimized)

#Now do the simulation 10000 times in parallel. Return an array
solArr = monteCarloSim(prob::SDEProblem,Δt=1//2^(3),numMonte=5)

#First index is the sime, so sol.timeseries[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.timeseries[..,1] returns the time series for the
#first row, and sol.timeseries[..,2] returns the time series for the second.
TEST_PLOT && plot(sol,plot_analytic=true)

## Convergence Testing
println("Convergence Test on MultiDimAdditive")
Δts = 1./2.^(7:-1:4) #14->7 good plot

sim = test_convergence(Δts,prob,numMonte=5,alg=:SRA)

sim2 = test_convergence(Δts,prob,numMonte=5,alg=:SRA1Optimized)

abs(sim.𝒪est[:l2]-2) + abs(sim2.𝒪est[:l∞]-2) <.1 #High tolerance since low Δts for testing!

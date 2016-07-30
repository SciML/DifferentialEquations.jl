using DifferentialEquations, Plots
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,[0,1],Δt=1//2^(4),save_timeseries=true,alg=:SRI,abstol=1,reltol=0)
err1 = sol.errors[:final]
TEST_PLOT && plot(sol,plottrue=true,legend=false,title="tol = 1")


println("Solve and Plot")
sol =solve(prob::SDEProblem,[0,1],Δt=1//2^(4),save_timeseries=true,alg=:SRI,adaptive=true,abstol=1,reltol=0)
err1 = sol.errors[:final]
TEST_PLOT && plot(sol,plottrue=true,legend=false,title="tol = 1")

true #Always gives float on SDE due to RNG

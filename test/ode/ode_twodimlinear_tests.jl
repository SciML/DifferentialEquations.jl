using DifferentialEquations
prob = twoDimlinearODEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:Euler)

#First index is the sime, so sol.timeseries[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.timeseries[..,1] returns the time series for the
#first row, and sol.timeseries[..,2] returns the time series for the second.
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:ExplicitRK)

true

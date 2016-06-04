using DifferentialEquations, Plots, EllipsisNotation
prob = lorenzAttractorODEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::ODEProblem,[0,100];Δt=1//2^(4),fullSave=true,alg="ExplicitRK",adaptive=true)

#First index is the sime, so sol.uFull[1,..] is the initial condition
#Last indices are the indexes of the variables. Since our initial condition
#Has 4 rows and two columns, sol.uFull[..,1] returns the time series for the
#first row, and sol.uFull[..,2] returns the time series for the second.
plot(sol.uFull[..,1],sol.uFull[..,2],sol.uFull[..,3])
#gui()

true

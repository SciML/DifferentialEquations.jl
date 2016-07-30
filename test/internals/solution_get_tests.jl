using DifferentialEquations,Plots
prob = twoDimlinearODEExample()

## Solve and plot
println("Test getindex")
tab = constructBogakiShampine()
sol =solve(prob::ODEProblem,save_timeseries=true,alg=:ExplicitRK,adaptive=true,tableau=tab)
sol[1] == prob.u₀ && sol[end] == sol.u

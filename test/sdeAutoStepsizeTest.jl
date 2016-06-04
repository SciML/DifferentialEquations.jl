using DifferentialEquations, Plots
srand(100)
prob = twoDimlinearSDEExample()

## Solve and plot
println("Solve and Plot")
sol =solve(prob::SDEProblem,fullSave=true,alg="SRI",adaptive=true)

#Now do the simulation 10000 times in parallel. Return an array
solArr = monteCarloSim(prob::SDEProblem,Δt=1//2^(4))

plot(sol,plottrue=true)
gui()

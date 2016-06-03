using DifferentialEquations,Plots
prob = twoDimlinearODEExample()

## Solve and plot
println("Solve and Plot")
tab = constructBogakiShampine()
sol =solve(prob::ODEProblem,fullSave=true,alg="ExplicitRK",adaptive=true,tableau=tab)
plot(sol,plottrue=true)
Δt₀ = sol.tFull[2]

1e-7 < Δt₀ < .01

sol =solve(prob::ODEProblem,fullSave=true,alg="Euler")
plot(sol,plottrue=true)
Δt₀ = sol.tFull[2]

1e-7 < Δt₀ < .01

tab = constructDormandPrince8()
sol3 =solve(prob::ODEProblem,fullSave=true,alg="ExplicitRK",adaptive=true,tableau=tab)
plot(sol3,plottrue=true)
Δt₀ = sol3.tFull[2]

1e-7 < Δt₀ < .01

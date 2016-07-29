using DifferentialEquations,Plots
prob = twoDimlinearODEExample()

## Solve and plot
println("Solve and Plot")
tab = constructBogakiShampine()
sol =solve(prob::ODEProblem,save_timeseries=true,alg=:ExplicitRK,adaptive=true,tableau=tab)
plot(sol,plottrue=true)
Δt₀ = sol.ts[2]

bool1 = 1e-7 < Δt₀ < .1

sol =solve(prob::ODEProblem,save_timeseries=true,alg=:Euler)
plot(sol,plottrue=true)
Δt₀ = sol.ts[2]

bool2 = 1e-7 < Δt₀ < .01

tab = constructDormandPrince8()
sol3 =solve(prob::ODEProblem,save_timeseries=true,alg=:ExplicitRK,adaptive=true,tableau=tab)
plot(sol3,plottrue=true)
Δt₀ = sol3.ts[2]

bool3 = 1e-7 < Δt₀ < .3

bool1 && bool2 && bool3

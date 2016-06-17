using DifferentialEquations,Plots
prob = twoDimlinearODEExample()

## Solve and plot
println("Solve and Plot")
tab = constructBogakiShampine()
sol =solve(prob::ODEProblem,[0,1],Δt=1/2^4,fullSave=true,alg="ExplicitRK",adaptive=true,tableau=tab)
plot(sol,plottrue=true)
#gui()

tab = constructDormandPrince()
sol2 =solve(prob::ODEProblem,[0,1],Δt=1/2^4,fullSave=true,alg="ExplicitRK",adaptive=true,tableau=tab)
plot(sol2,plottrue=true)
#gui()

tab = constructRKF8()
sol3 =solve(prob::ODEProblem,[0,1],Δt=1/2^4,fullSave=true,alg="ExplicitRK",adaptive=true,tableau=tab)
plot(sol3,plottrue=true)
#gui()

length(sol.tFull)>length(sol2.tFull)>length(sol3.tFull)

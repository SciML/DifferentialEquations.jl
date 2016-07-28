# Introduction to the ODE Solvers

using DifferentialEquations

prob = linearODEExample()
Δt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.
sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:dopri5)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:dop853)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:odex)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:seulex)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:radau)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:radau5)
plot(sol,plottrue=true)

prob = twoDimlinearODEExample()

sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:dopri5)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:dop853)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:odex)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:seulex)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:radau)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),fullSave=true,alg=:radau5)
plot(sol,plottrue=true)

true

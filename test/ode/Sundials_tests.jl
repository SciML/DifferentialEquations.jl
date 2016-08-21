using DifferentialEquations

prob = linearODEExample()
Δt = 1//2^(4)
tspan = 0:1/2^(4):1
@time sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:cvode)

prob = twoDimlinearODEExample()
@time sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:cvode)

prob = twoDimlinearODEExample!()
@time sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:cvode)

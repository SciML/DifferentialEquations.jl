# Introduction to the ODE Solvers

using DifferentialEquations

prob = prob_ode_linear
Δt = 1//2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.
sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:dopri5)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:dop853)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:odex)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:seulex)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:radau)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:radau5)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:dopri5)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:dop853)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:odex)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:seulex)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:radau)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:radau5)

#=
prob = prob_ode_bigfloat2Dlinear

sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:dopri5)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:dop853)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:odex)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:seulex)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:radau)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:radau5)
=#

true

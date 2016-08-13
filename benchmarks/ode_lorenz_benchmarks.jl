using DifferentialEquations, Plots, EllipsisNotation, ODE

prob = linearODEExample()
using BenchmarkTools
sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:Rosenbrock32,save_timeseries=true,timeseries_steps=1)
@benchmark sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)
@benchmark sol2 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,save_timeseries=true,timeseries_steps=1)

prob = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
prob = twoDimlinearODEExample!()
sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:Rosenbrock32,save_timeseries=true,timeseries_steps=1)
@benchmark sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)

prob2 = twoDimlinearODEExample(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
@benchmark sol2 =solve(prob2::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,abstol=1e-6,reltol=1e-6)

@benchmark sol1 =solve(prob2::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)

sol =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:ExplicitRK,abstol=1e-2,reltol=1e-2)
sol =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:ode45,abstol=1e-2,reltol=1e-2)

#prob = lorenzAttractorODEExample()

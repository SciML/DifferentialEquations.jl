using DifferentialEquations, Plots, EllipsisNotation, ODE
#=
prob = linearODEExample()
using BenchmarkTools
b1 = @benchmark sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)
b2 = @benchmark sol2 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,save_timeseries=true,timeseries_steps=1)
=#

prob = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)

#b3 = @benchmark sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)

prob2 = twoDimlinearODEExample(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
#b4 = @benchmark sol2 =solve(prob2::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,abstol=1e-6,reltol=1e-6)

#b5 = @benchmark sol1 =solve(prob2::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)
sol =solve(prob::ODEProblem,[0,10],alg=:ExplicitRK,abstol=1e-6,reltol=1e-6,progressbar=true,progress_steps=1)

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:ExplicitRK,γ=0.8)
elapsed2 = @elapsed sol2 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:ExplicitRKVectorized,γ=0.8)

elapsed3 = @elapsed sol3 =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:ode45)
elapsed4 = @elapsed sol4 =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:dopri5)

eff1 = sol1.errors[:final]/elapsed1
eff2 = sol2.errors[:final]/elapsed2
eff3 = sol3.errors[:final]/elapsed3
eff4 = sol4.errors[:final]/elapsed4

eff1/eff3
eff4/eff1
eff4/eff3

#prob = lorenzAttractorODEExample()

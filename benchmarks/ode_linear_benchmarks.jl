using DifferentialEquations, Plots, EllipsisNotation, ODE
probnum = linearODEExample()
prob = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
prob2 = twoDimlinearODEExample(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
using BenchmarkTools



b1 = @benchmark sol1 =solve(probnum::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)
b2 = @benchmark sol2 =solve(probnum::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,save_timeseries=true,timeseries_steps=1)
b3 = @benchmark sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)
b4 = @benchmark sol2 =solve(prob2::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,abstol=1e-6,reltol=1e-6)
b5 = @benchmark sol1 =solve(prob2::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=true,timeseries_steps=1)

# Test Progressbar

sol =solve(prob::ODEProblem,[0,10],Δt=1/2^(4),alg=:DP5,abstol=1e-6,reltol=1e-6,progressbar=true,progress_steps=1)

# Precompile

sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:DP5,timechoicealg=:Simple)
sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:DP5,timechoicealg=:Lund)
sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:ExplicitRK,timechoicealg=:Lund)
sol3 =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:ode45)
sol4 =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:dopri5)

## Standard Tolerance

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];alg=:DP5)
elapsed3 = @elapsed sol3 =solve(prob2::ODEProblem,[0,10];abstol=1e-3,roltol=1e-6,alg=:ode45) # Fix ODE.jl to be normal...
elapsed4 = @elapsed sol4 =solve(prob2::ODEProblem,[0,10];alg=:dopri5)

## Test how much choices matter
elapsed2 = @elapsed sol2 =solve(prob::ODEProblem,[0,10];alg=:DP5Vectorized)
elapsed5 = @elapsed sol5 =solve(prob::ODEProblem,[0,10];alg=:DP5,timechoicealg=:Simple)
elapsed6 = @elapsed sol6 =solve(prob::ODEProblem,[0,10];fullnormalize=true,alg=:DP5)
elapsed7 = @elapsed sol7 =solve(prob::ODEProblem,[0,10];alg=:ExplicitRK)

eff1 = 1/(sol1.errors[:final]*elapsed1)
eff2 = 1/(sol2.errors[:final]*elapsed2)
eff3 = 1/(sol3.errors[:final]*elapsed3)
eff4 = 1/(sol4.errors[:final]*elapsed4)
eff5 = 1/(sol5.errors[:final]*elapsed5)
eff6 = 1/(sol6.errors[:final]*elapsed6)
eff7 = 1/(sol7.errors[:final]*elapsed7)

rat2 = eff1/eff3
rat3 = eff1/eff4
rat5 = eff4/eff3

#Choices
rat1 = eff1/eff2
rat4 = eff1/eff5
rat5 = eff1/eff6
rat6 = eff1/eff7

## Test at ODE.jl defaults

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];abstol=1e-8,reltol=1e-5,alg=:DP5)
elapsed3 = @elapsed sol3 =solve(prob2::ODEProblem,[0,10];alg=:ode45) # Fix ODE.jl to be normal...
eff1 = 1/(sol1.errors[:final]*elapsed1); eff3 = 1/(sol3.errors[:final]*elapsed3)
rat2 = eff1/eff3

## Low Tolerance

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];reltol=1e-6,alg=:DP5)
elapsed2 = @elapsed sol2 =solve(prob::ODEProblem,[0,10];reltol=1e-6,alg=:DP5Vectorized)

elapsed3 = @elapsed sol3 =solve(prob2::ODEProblem,[0,10];abstol=1e-6,roltol=1e-6,alg=:ode45)
elapsed4 = @elapsed sol4 =solve(prob2::ODEProblem,[0,10];abstol=1e-6,roltol=1e-6,alg=:dopri5)

elapsed5 = @elapsed sol5 =solve(prob::ODEProblem,[0,10];reltol=1e-6,alg=:ExplicitRK,timechoicealg=:Simple)

eff1 = 1/(sol1.errors[:final]*elapsed1)
eff2 = 1/(sol2.errors[:final]*elapsed2)
eff3 = 1/(sol3.errors[:final]*elapsed3)
eff4 = 1/(sol4.errors[:final]*elapsed4)
eff5 = 1/(sol5.errors[:final]*elapsed5)

rat1 = eff1/eff2
rat2 = eff1/eff3
rat3 = eff1/eff4
rat4 = eff1/eff5
rat5 = eff4/eff3

## Number

elapsed1 = @elapsed sol1 =solve(probnum::ODEProblem,[0,10];reltol=1e-6,alg=:DP5)
elapsed2 = @elapsed sol2 =solve(probnum::ODEProblem,[0,10];reltol=1e-6,alg=:ExplicitRK)

elapsed3 = @elapsed sol3 =solve(probnum::ODEProblem,[0,10];abstol=1e-6,roltol=1e-6,alg=:ode45)
elapsed4 = @elapsed sol4 =solve(probnum::ODEProblem,[0,10];reltol=1e-6,alg=:dopri5)

elapsed5 = @elapsed sol5 =solve(probnum::ODEProblem,[0,10];reltol=1e-6,alg=:DP5,timechoicealg=:Simple)

eff1 = 1/(sol1.errors[:final]*elapsed1)
eff2 = 1/(sol2.errors[:final]*elapsed2)
eff3 = 1/(sol3.errors[:final]*elapsed3)
eff4 = 1/(sol4.errors[:final]*elapsed4)
eff5 = 1/(sol5.errors[:final]*elapsed5)

rat1 = eff1/eff2
rat2 = eff1/eff3
rat3 = eff1/eff4
rat4 = eff1/eff5
rat5 = eff4/eff3

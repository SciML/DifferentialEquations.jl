### Linear

using DifferentialEquations
probnum = linearODEExample()
prob = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
prob2 = twoDimlinearODEExample(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
using BenchmarkTools



b1 = @benchmark sol1 =solve(probnum::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,timeseries_steps=1)
b2 = @benchmark sol2 =solve(probnum::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4,timeseries_steps=1)
b4 = @benchmark sol2 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:ode4)
b5 = @benchmark sol1 =solve(prob::ODEProblem,[0,1];Δt=1/2^(10),alg=:RK4,save_timeseries=false)

# Test Progressbar

sol =solve(prob::ODEProblem,[0,1,10],Δt=1/2^(4),alg=:DP5,abstol=1e-6,reltol=1e-3,progressbar=true,progress_steps=1)

# Precompile

sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:DP5,timechoicealg=:Simple)
sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:DP5,timechoicealg=:Lund)
sol1 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:ExplicitRK,timechoicealg=:Lund)
sol3 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:ode45)
sol3 =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:ode45)
sol4 =solve(prob::ODEProblem,[0,10];Δt=1/2^(4),alg=:dopri5)
sol4 =solve(prob2::ODEProblem,[0,10];Δt=1/2^(4),alg=:dopri5)
## Standard Tolerance

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];alg=:DP5)
elapsed3 = @elapsed sol3 =solve(prob2::ODEProblem,[0,10];abstol=1e-3,roltol=1e-6,alg=:ode45) # Fix ODE.jl to be normal...
elapsed4 = @elapsed sol4 =solve(prob2::ODEProblem,[0,10];alg=:dopri5)

## Test how much choices matter
elapsed2 = @elapsed sol2 =solve(prob::ODEProblem,[0,10];alg=:DP5Vectorized)
elapsed5 = @elapsed sol5 =solve(prob::ODEProblem,[0,10];alg=:DP5,timechoicealg=:Simple)
elapsed6 = @elapsed sol6 =solve(prob::ODEProblem,[0,10];fullnormalize=true,alg=:DP5)
elapsed7 = @elapsed sol7 =solve(prob::ODEProblem,[0,10];alg=:ExplicitRK) # Higher β by default
elapsed8 = @elapsed sol8 =solve(prob::ODEProblem,[0,10];alg=:DP5,β=0.04) # Higher β by default
elapsed9 = @elapsed sol9 =solve(prob::ODEProblem,[0,10];alg=:BS3,β=0.04) # Higher β by default
elapsed10 = @elapsed sol10 =solve(prob::ODEProblem,[0,10];alg=:cvode_Adams)

eff1 = 1/(sol1.errors[:final]*elapsed1)
eff2 = 1/(sol2.errors[:final]*elapsed2)
eff3 = 1/(sol3.errors[:final]*elapsed3)
eff4 = 1/(sol4.errors[:final]*elapsed4)
eff5 = 1/(sol5.errors[:final]*elapsed5)
eff6 = 1/(sol6.errors[:final]*elapsed6)
eff7 = 1/(sol7.errors[:final]*elapsed7)
eff8 = 1/(sol8.errors[:final]*elapsed8)
eff9 = 1/(sol9.errors[:final]*elapsed9)
eff10= 1/(sol10.errors[:final]*elapsed10)

rat2 = eff1/eff3
rat3 = eff1/eff4
rat5 = eff4/eff3

#Choices
rat1 = eff1/eff2
rat4 = eff1/eff5
rat5 = eff1/eff6
rat6 = eff1/eff7
rat7 = eff1/eff8

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

# Other Methods

t5 = @elapsed sol5 =solve(prob::ODEProblem,[0,10],alg=:BS3,reltol=1e-6)
e5 = sol5.errors[:final]
eff5 = 1/(t5*e5)

t6 = @elapsed sol5 =solve(prob::ODEProblem,[0,10],alg=:BS5,reltol=1e-6)
e6 = sol5.errors[:final]
eff6 = 1/(t6*e6)

t7 = @elapsed sol5 =solve(prob::ODEProblem,[0,10],alg=:Tsit5,reltol=1e-6)
e7 = sol5.errors[:final]
eff7 = 1/(t7*e7)

eff1/eff5
eff1/eff6
eff1/eff7

# Higher Order

t8 = @elapsed sol8 =solve(prob::ODEProblem,[0,10],alg=:DP8,reltol=1e-6,β=0.07)
e8 = sol8.errors[:final]
eff8 = 1/(t8*e8)

t9 = @elapsed sol9 =solve(prob::ODEProblem,[0,10],alg=:dop853,reltol=1e-6)
e9 = sol9.errors[:final]
eff9 = 1/(t9*e9)

t9 = @elapsed sol9 =solve(prob::ODEProblem,[0,10],alg=:ode78,reltol=1e-6)
e9 = sol9.errors[:final]
eff9 = 1/(t9*e9)

## Number

elapsed1 = @elapsed sol1 =solve(probnum::ODEProblem,[0,10];reltol=1e-6,alg=:DP5)
elapsed2 = @elapsed sol2 =solve(probnum::ODEProblem,[0,10];reltol=1e-6,alg=:ExplicitRK,β=0.06)

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

##### Three-Body

using DifferentialEquations
## Define the ThreeBody Problem
const threebody_μ = parse(BigFloat,"0.012277471"); const threebody_μ′ = 1 - threebody_μ

f = (t,u,du) -> begin
  # 1 = y₁
  # 2 = y₂
  # 3 = y₁'
  # 4 = y₂'
  D₁ = ((u[1]+threebody_μ)^2 + u[2]^2)^(3/2)
  D₂ = ((u[1]-threebody_μ′)^2 + u[2]^2)^(3/2)
  du[1] = u[3]
  du[2] = u[4]
  du[3] = u[1] + 2u[4] - threebody_μ′*(u[1]+threebody_μ)/D₁ - threebody_μ*(u[1]-threebody_μ′)/D₂
  du[4] = u[2] - 2u[3] - threebody_μ′*u[2]/D₁ - threebody_μ*u[2]/D₂
end

prob_ode_threebody = ODEProblem(f,[0.994, 0.0, 0.0, parse(BigFloat,"-2.00158510637908252240537862224")])
prob = prob_ode_threebody

t₀ = 0.0; T = parse(BigFloat,"17.0652165601579625588917206249")

### Progress bar run
sol =solve(prob::ODEProblem,[t₀,T],alg=:DP5,progressbar=true,progress_steps=1)
### Warmups
solve(prob::ODEProblem,[t₀,T],alg=:DP5)
solve(prob::ODEProblem,[t₀,T],alg=:dopri5)
solve(prob::ODEProblem,[t₀,T],alg=:ode45,abstol=1e-6,reltol=1e-3)

t1 = @elapsed sol =solve(prob::ODEProblem,[t₀,T],alg=:DP5,abstol=1e-6,reltol=1e-3)
e1 = norm(sol.u-prob.u₀)
eff1 = 1/(t1*e1)

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,T],alg=:ode45,abstol=1e-6,reltol=1e-3)
e2 = norm(sol2.u-prob.u₀)
eff2 = 1/(t2*e2)

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,T],alg=:ode45,abstol=1e-6,reltol=1e-3,norm=(y)->vecnorm(y,2),minstep=1e-40) # Fails
e2 = norm(sol2.u-prob.u₀)
eff2 = 1/(t2*e2)

t3 = @elapsed sol3 =solve(prob::ODEProblem,[t₀,T],alg=:dopri5,abstol=1e-6,reltol=1e-3)
e3 = norm(sol3.u-prob.u₀)
eff3 = 1/(t3*e3)

t4 = @elapsed sol4 =solve(prob::ODEProblem,[t₀,T],alg=:ExplicitRK,abstol=1e-6,reltol=1e-3)
e4 = norm(sol4.u-prob.u₀)
eff4 = 1/(t4*e4)

t4 = @elapsed sol4 =solve(prob::ODEProblem,[t₀,T],alg=:BS5,abstol=1e-6,reltol=1e-3)
e4 = norm(sol4.u-prob.u₀)
eff4 = 1/(t4*e4)

eff1/eff2
eff1/eff3
eff1/eff4

### Lower Tolerance

t1 = @elapsed sol =solve(prob::ODEProblem,[t₀,T],alg=:DP5,abstol=1e-9,reltol=1e-5)
e1 = norm(sol.u-prob.u₀)
eff1 = 1/(t1*e1)

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,T],alg=:ode45,abstol=1e-9,reltol=1e-5)
e2 = norm(sol2.u-prob.u₀)
eff2 = 1/(t2*e2)

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,T],Δt=sol.t[2],alg=:ode45,abstol=1e-9,reltol=1e-5,norm=(y)->vecnorm(y,2))
e2 = norm(sol2.u-prob.u₀)
eff2 = 1/(t2*e2)

t3 = @elapsed sol3 =solve(prob::ODEProblem,[t₀,T],alg=:dopri5,abstol=1e-9,reltol=1e-5)
e3 = norm(sol3.u-prob.u₀)
eff3 = 1/(t3*e3)

t4 = @elapsed sol4 =solve(prob::ODEProblem,[t₀,T],alg=:DP5,abstol=1e-9,reltol=1e-5,β=0.065)
e4 = norm(sol4.u-prob.u₀)
eff4 = 1/(t4*e4)

t5 = @elapsed sol5 =solve(prob::ODEProblem,[t₀,T],alg=:BS5,abstol=1e-9,reltol=1e-5)
e5 = norm(sol5.u-prob.u₀)
eff5 = 1/(t5*e5)

eff1/eff2
eff1/eff3
eff1/eff4
eff1/eff5
eff4/eff2

### Longer

t1 = @elapsed sol =solve(prob::ODEProblem,[t₀,2T],lg=:DP5,abstol=1e-9,reltol=1e-7)
e1 = norm(sol.u-prob.u₀)
eff1 = 1/(t1*e1)

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,2T],Δt=sol.t[2],alg=:ode45,abstol=1e-9,reltol=1e-7)
e2 = norm(sol2.u-prob.u₀)
eff2 = 1/(t2*e2)

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,2T],Δt=sol.t[2],alg=:ode45,abstol=1e-9,reltol=1e-7,norm=(y)->vecnorm(y,2))
e2 = norm(sol2.u-prob.u₀)
eff2 = 1/(t2*e2)

t3 = @elapsed sol3 =solve(prob::ODEProblem,[t₀,2T],alg=:dopri5,abstol=1e-9,reltol=1e-7)
e3 = norm(sol3.u-prob.u₀)
eff3 = 1/(t3*e3)

t4 = @elapsed sol =solve(prob::ODEProblem,[t₀,2T],alg=:BS5,abstol=1e-9,reltol=1e-7)
e4 = norm(sol.u-prob.u₀)
eff4 = 1/(t4*e4)

t5 = @elapsed sol =solve(prob::ODEProblem,[t₀,2T],alg=:Tsit5,abstol=1e-9,reltol=1e-7)
e5 = norm(sol.u-prob.u₀)
eff5 = 1/(t5*e5)

eff1/eff2
eff1/eff3
eff1/eff4

t5 = @elapsed sol5 =solve(prob::ODEProblem,[t₀,2T],alg=:Feagin14,abstol=1e-9,reltol=1e-7)
e5 = norm(sol5.u-prob.u₀)
eff5 = 1/(t5*e5)

t6 = @elapsed sol6 =solve(prob::ODEProblem,[t₀,2T],alg=:ode78,abstol=1e-9,reltol=1e-7)
e6 = norm(sol6.u-prob.u₀)
eff6 = 1/(t6*e6)

t5 = @elapsed sol5 =solve(prob::ODEProblem,[t₀,2T],alg=:ExplicitRK,abstol=1e-9,reltol=1e-7)
e5 = norm(sol5.u-prob.u₀)
eff5 = 1/(t5*e5)

tab = constructBogakiShampine5()
t5 = @elapsed sol5 =solve(prob::ODEProblem,[t₀,2T],alg=:ExplicitRK,abstol=1e-9,reltol=1e-7,tableau=tab)
e5 = norm(sol5.u-prob.u₀)
eff5 = 1/(t5*e5)

tab = constructTsitouras5()
t5 = @elapsed sol5 =solve(prob::ODEProblem,[t₀,2T],alg=:ExplicitRK,abstol=1e-9,reltol=1e-7,tableau=tab)
e5 = norm(sol5.u-prob.u₀)
eff5 = 1/(t5*e5)

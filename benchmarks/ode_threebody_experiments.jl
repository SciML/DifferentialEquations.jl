using DifferentialEquations
prob = threebodyODEExample()

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

t2 = @elapsed sol2 =solve(prob::ODEProblem,[t₀,T],Δt=sol.t[2],alg=:ode45,abstol=1e-6,reltol=1e-3,norm=(y)->vecnorm(y,2),minstep=1e-40) # Fails
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

using DifferentialEquations

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u) -> (linear_bigα*u)
analytic = (t,u₀) -> u₀*exp(linear_bigα*t)
"""Linear ODE on Float64"""
prob_ode_bigfloatlinear = ODEProblem(f,parse(BigFloat,"0.5"),analytic=analytic)

f = (t,u,du) -> begin
for i in 1:length(u)
  du[i] = linear_bigα*u[i]
end
end
"""2D Linear ODE, bigfloats"""
prob_ode_bigfloat2Dlinear = ODEProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,analytic=analytic)

probnum = prob_ode_linear
probnumbig = prob_ode_bigfloatlinear
#prob    = prob_ode_large2Dlinear
prob = prob_ode_2Dlinear
probbig = prob_ode_bigfloat2Dlinear
Δts = 1.//2.^(7:-1:4)
testTol = .2
bools = Vector{Bool}(0)

## DP5

sim = test_convergence(Δts,probnum,alg=:DP5)
push!(bools,abs(sim.𝒪est[:l2]-5) < testTol)
sim = test_convergence(Δts,prob,alg=:DP5Vectorized)
push!(bools,abs(sim.𝒪est[:l2]-5) < testTol)
sim = test_convergence(Δts,prob,alg=:DP5)
push!(bools,abs(sim.𝒪est[:l2]-5) < testTol)

sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,adaptive=false,save_timeseries=false)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5Vectorized,adaptive=false,save_timeseries=false)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRKVectorized,adaptive=false,save_timeseries=false)

push!(bools,minimum(sol1.u - sol2.u .< 3e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,adaptive=false,save_timeseries=false)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false)

push!(bools,minimum(sol1.u - sol2.u .< 3e-10))

sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,β=0.04)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,β=0.04)


sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,β=0.04)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,β=0.04)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRKVectorized,β=0.04)

### BS3
sim = test_convergence(Δts,probnum,alg=:BS3)
push!(bools,abs(sim.𝒪est[:l2]-3) < testTol)
sim = test_convergence(Δts,prob,alg=:BS3Vectorized)
push!(bools,abs(sim.𝒪est[:l2]-3) < testTol)
sim = test_convergence(Δts,prob,alg=:BS3)
push!(bools,abs(sim.𝒪est[:l2]-3) < testTol)

tab = constructBogakiShampine3()
sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^1,alg=:BS3,adaptive=false,save_timeseries=false)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^1,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:BS3Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRKVectorized,tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^1,alg=:BS3,adaptive=false,save_timeseries=false)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^1,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,2],Δt=1/2^6,alg=:BS3Vectorized)
sol2 =solve(prob::ODEProblem,[0,2],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,2],Δt=1/2^6,alg=:BS3)

### BS5
Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:BS5)

sim = test_convergence(Δts,prob,alg=:BS5)

tab = constructBogakiShampine5()
sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5,adaptive=false,save_timeseries=false)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:BS5,adaptive=false,save_timeseries=false)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5)

### Tsit5

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnum,alg=:Tsit5)
#sim = test_convergence(Δts,prob,alg=:Tsit5)

tab = constructTsitouras5()
sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:Tsit5,adaptive=false,save_timeseries=false)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:Tsit5)
sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:Tsit5,adaptive=false,save_timeseries=false)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Tsit5Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Tsit5)

### Vern6

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnumbig,alg=:Vern6)
#sim = test_convergence(Δts,probbig,alg=:Vern6)

tab = constructVernerEfficient6(BigFloat)
sol1 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern6,adaptive=false,save_timeseries=false)
sol2 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:Vern6)
sol2 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:Vern6,adaptive=false,save_timeseries=false)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern6Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern6)

### Vern7

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnumbig,alg=:Vern7)
#sim = test_convergence(Δts,probbig,alg=:Vern7)

tab = constructVerner7(BigFloat)
sol1 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern7,adaptive=false,save_timeseries=false)
sol2 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:Vern7)
sol2 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:Vern7,adaptive=false,save_timeseries=false)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern7Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern7)

### TanYam7

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnumbig,alg=:TanYam7)
#sim = test_convergence(Δts,probbig,alg=:TanYam7)

tab = constructTanakaYamashitaEfficient7(BigFloat)
sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:TanYam7,adaptive=false,save_timeseries=false)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:TanYam7)
sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:TanYam7,adaptive=false,save_timeseries=false)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TanYam7Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TanYam7)

### Vern8

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnumbig,alg=:Vern8)
#sim = test_convergence(Δts,probbig,alg=:Vern8)

tab = constructVerner8(BigFloat)
sol1 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern8,adaptive=false,save_timeseries=false)
sol2 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:Vern8)
sol2 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:Vern8,adaptive=false,save_timeseries=false)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern8Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern8)

### DP8

Δts = 1.//2.^(3:-1:1)
#sim = test_convergence(Δts,probnumbig,alg=:DP8)
#sim = test_convergence(Δts,probbig,alg=:DP8)

sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8,adaptive=false,save_timeseries=false)

sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:DP8,adaptive=false,save_timeseries=false)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8Vectorized)
sol3 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8)
#sol4 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:dop853)

### TsitPap8

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnumbig,alg=:TsitPap8)
#sim = test_convergence(Δts,probbig,alg=:TsitPap8)

tab = constructTsitourasPapakostas8(BigFloat)
sol1 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:TsitPap8,adaptive=false,save_timeseries=false)
sol2 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,sol1.u - sol2.u < 1e-10)

sol1 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:TsitPap8)
sol2 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:TsitPap8,adaptive=false,save_timeseries=false)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(sol1.u - sol2.u .< 1e-10))

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TsitPap8Vectorized)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TsitPap8)

### Vern9

Δts = 1.//2.^(6:-1:3)
#sim = test_convergence(Δts,probnumbig,alg=:Vern9)
#sim = test_convergence(Δts,probbig,alg=:Vern9)

tab = constructVernerEfficient9(BigFloat)
sol1 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern9,adaptive=false,save_timeseries=false)
sol2 =solve(probnumbig::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,abs(sol1.u - sol2.u) < 1e-15)

sol1 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:Vern9)
sol2 =solve(probnumbig::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:Vern9,adaptive=false,save_timeseries=false)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

push!(bools,minimum(abs(sol1.u - sol2.u) .< 1e-15))

sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern9Vectorized)
sol2 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern9)

minimum(bools)

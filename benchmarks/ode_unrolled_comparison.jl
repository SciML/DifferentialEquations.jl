using DifferentialEquations
using BenchmarkTools
probnum = linearODEExample(); probnumbig = linearODEExample(u₀=BigFloat(1.0))
prob    = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
probbig = twoDimlinearODEExample!(α=ones(BigFloat,4,2),u₀=map(BigFloat,rand(4,2)).*ones(4,2)/2)
prob2   = twoDimlinearODEExample(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
Δts = 1.//2.^(7:-1:4)

## DP5

sim = test_convergence(Δts,probnum,alg=:DP5)
sim = test_convergence(Δts,prob,alg=:DP5)

@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5Vectorized,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRKVectorized,adaptive=false,save_timeseries=false)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false)

sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5)
sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,β=0.04)


sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK)
sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRKVectorized)

### BS3
sim = test_convergence(Δts,probnum,alg=:BS3)
sim = test_convergence(Δts,prob,alg=:BS3)

tab = constructBogakiShampine3()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^1,alg=:BS3,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^1,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:BS3)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^1,alg=:BS3,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^1,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

sol1 =solve(prob::ODEProblem,[0,2],Δt=1/2^6,alg=:BS3Vectorized)
sol2 =solve(prob::ODEProblem,[0,2],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
sol3 =solve(prob::ODEProblem,[0,2],Δt=1/2^6,alg=:BS3)

### BS5
Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:BS5)
sim = test_convergence(Δts,prob,alg=:BS5)

tab = constructBogakiShampine5()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:BS5,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5Vectorized)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
@time sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:BS5)

### Tsit5

Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:Tsit5)
sim = test_convergence(Δts,prob,alg=:Tsit5)

tab = constructTsitouras5()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:Tsit5,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:Tsit5)
@time sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:Tsit5,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Tsit5Vectorized)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
@time sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Tsit5)

### Vern6

Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:Vern6)
sim = test_convergence(Δts,prob,alg=:Vern6)

tab = constructTsitouras5()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern6,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:Vern6)
@time sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:Vern6,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern6Vectorized)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
@time sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern6)

### TanYam7

Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:TanYam7)
sim = test_convergence(Δts,prob,alg=:TanYam7)

tab = constructTsitouras5()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:TanYam7,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:TanYam7)
@time sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:TanYam7,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TanYam7Vectorized)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
@time sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TanYam7)

### DP8

Δts = 1.//2.^(3:-1:1)
sim = test_convergence(Δts,probnumbig,alg=:DP8)
sim = test_convergence(Δts,probbig,alg=:DP8)

@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8,adaptive=false,save_timeseries=false)

@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8)

@time sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^3,alg=:DP8,adaptive=false,save_timeseries=false)

@time sol1 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8Vectorized)
@time sol3 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:DP8)
@time sol4 =solve(probbig::ODEProblem,[0,10],Δt=1/2^6,alg=:dop853)

### TsitPap8

Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:TsitPap8)
sim = test_convergence(Δts,prob,alg=:TsitPap8)

tab = constructTsitouras5()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:TsitPap8,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:TsitPap8)
@time sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:TsitPap8,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TsitPap8Vectorized)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
@time sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:TsitPap8)

### Vern9

Δts = 1.//2.^(6:-1:3)
sim = test_convergence(Δts,probnum,alg=:Vern9)
sim = test_convergence(Δts,prob,alg=:Vern9)

tab = constructTsitouras5()
@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern9,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

sol1.u - sol2.u < 1e-10

@time sol1 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:Vern9)
@time sol2 =solve(probnum::ODEProblem,[0,70],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:Vern9,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^3,alg=:ExplicitRK,adaptive=false,save_timeseries=false, tableau=tab)

minimum(sol1.u - sol2.u .< 1e-10)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern9Vectorized)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,tableau=tab)
@time sol3 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:Vern9)

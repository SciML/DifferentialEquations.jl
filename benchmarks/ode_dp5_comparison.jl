using DifferentialEquations
using BenchmarkTools
probnum = linearODEExample()
prob = twoDimlinearODEExample!(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)
prob2 = twoDimlinearODEExample(;α=ones(100,100),u₀=rand(100,100).*ones(100,100)/2)

@time sol1 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,adaptive=false,save_timeseries=false)
@time sol2 =solve(probnum::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5Vectorized,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRKVectorized,adaptive=false,save_timeseries=false)

@time sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5,adaptive=false,save_timeseries=false)
@time sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK,adaptive=false,save_timeseries=false)

sol1 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:DP5)
sol2 =solve(prob::ODEProblem,[0,10],Δt=1/2^6,alg=:ExplicitRK)

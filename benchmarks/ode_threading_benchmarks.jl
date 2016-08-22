using DifferentialEquations
probnum = linearODEExample()
prob = twoDimlinearODEExample!(;α=ones(500,500),u₀=rand(500,500).*ones(500,500)/2)
prob2 = twoDimlinearODEExample(;α=ones(500,500),u₀=rand(500,500).*ones(500,500)/2)
using BenchmarkTools

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];alg=:DP5)
elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];alg=:DP5Threaded)

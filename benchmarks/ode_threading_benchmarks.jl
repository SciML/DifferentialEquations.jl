using DifferentialEquations
probnum = prob_ode_linear
prob = prob_ode_2Dlinear
using BenchmarkTools

elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];alg=:DP5)
elapsed1 = @elapsed sol1 =solve(prob::ODEProblem,[0,10];alg=:DP5Threaded)

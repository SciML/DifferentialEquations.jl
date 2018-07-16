using DifferentialEquations, Test

lags = [.2]
f = function (u,h,p,t)
  out = -h(p,t-.2) + u
end
h = (p,t) -> 0.0


prob = DDEProblem(f,1.0,h,(0.0,10.0),constant_lags = lags)

sol = solve(prob)

@test typeof(sol.alg) <: CompositeAlgorithm
@test typeof(sol.alg.algs[1]) <: Tsit5
@test typeof(sol.alg.algs[2]) <: Rosenbrock23

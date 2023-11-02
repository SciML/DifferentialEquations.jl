using DifferentialEquations, Test

lags = [0.2]
f = function (u, h, p, t)
    out = -h(p, t - 0.2) + u
end
h = (p, t) -> 0.0

prob = DDEProblem(f, 1.0, h, (0.0, 10.0), constant_lags = lags)

sol = solve(prob)

@test sol.alg isa CompositeAlgorithm
@test sol.alg.algs[1] isa Tsit5
@test sol.alg.algs[2] isa Rosenbrock23

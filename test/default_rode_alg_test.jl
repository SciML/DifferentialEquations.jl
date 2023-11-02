using DifferentialEquations, Test

f = (u, p, t, W) -> 1.01u .+ 0.87u .* W
u0 = 1.00
tspan = (0.0, 1.0)
prob = RODEProblem(f, u0, tspan)
sol = solve(prob, dt = 1 / 100)

@test sol.alg isa RandomEM

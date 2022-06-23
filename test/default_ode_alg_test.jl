using DifferentialEquations, Test

f_2dlinear = (du, u, p, t) -> (@. du = p * u)
f_2dlinear_analytic = (u0, p, t) -> @. u0 * exp(p * t)
prob_ode_2Dlinear = ODEProblem(ODEFunction(f_2dlinear, analytic = f_2dlinear_analytic),
                               rand(4, 2), (0.0, 1.0), 1.01)

alg, kwargs = default_algorithm(prob_ode_2Dlinear; dt = 1 // 2^(4))
sol = solve(prob_ode_2Dlinear; dt = 1 // 2^(4))

@test typeof(sol.alg.algs[1]) <: Tsit5
@test typeof(sol.alg.algs[2]) <: Rosenbrock23

sol = solve(prob_ode_2Dlinear; reltol = 1e-1)

@test typeof(sol.alg.algs[1]) <: Tsit5
@test typeof(sol.alg.algs[2]) <: Rosenbrock23

sol = solve(prob_ode_2Dlinear; reltol = 1e-7)

@test typeof(sol.alg.algs[1]) <: Vern7
@test typeof(sol.alg.algs[2]) <: Rodas4

sol = solve(prob_ode_2Dlinear; reltol = 1e-10)

@test typeof(sol.alg.algs[1]) <: Vern9
@test typeof(sol.alg.algs[2]) <: Rodas5

sol = solve(prob_ode_2Dlinear; alg_hints = [:stiff])

@test typeof(sol.alg) <: Rodas4

sol = solve(prob_ode_2Dlinear; alg_hints = [:stiff], reltol = 1e-1)

@test typeof(sol.alg) <: Rosenbrock23

const linear_bigα = parse(BigFloat, "1.01")
f = (du, u, p, t) -> begin for i in 1:length(u)
    du[i] = linear_bigα * u[i]
end end
(::typeof(f))(::Type{Val{:analytic}}, u0, p, t) = u0 * exp(linear_bigα * t)
prob_ode_bigfloat2Dlinear = ODEProblem(f, map(BigFloat, rand(4, 2)) .* ones(4, 2) / 2,
                                       (0.0, 1.0))

sol = solve(prob_ode_bigfloat2Dlinear; dt = 1 // 2^(4))
@test typeof(sol.alg.algs[1]) <: Vern9
@test typeof(sol.alg.algs[2]) <: Rodas5

default_algorithm(prob_ode_bigfloat2Dlinear; alg_hints = [:stiff])

sol = solve(prob_ode_bigfloat2Dlinear; alg_hints = [:stiff])

@test typeof(sol.alg) <: Rodas4

sol = solve(prob_ode_bigfloat2Dlinear, nothing; alg_hints = [:stiff])

@test typeof(sol.alg) <: Rodas4

struct FooAlg end

@test_throws ErrorException solve(prob_ode_bigfloat2Dlinear, FooAlg(); default_set = true)

struct FooAlg2 <: DiffEqBase.DEAlgorithm end

@test_throws ErrorException solve(prob_ode_bigfloat2Dlinear, FooAlg2(); default_set = true)

prob = ODEProblem(f, rand(4, 2) .* ones(4, 2) / 2, (0.0, 1.0))

sol = solve(prob; alg_hints = [:stiff])

@test typeof(sol.alg) <: Rodas4

sol = solve(prob; alg_hints = [:stiff], reltol = 1e-1)

@test typeof(sol.alg) <: Rosenbrock23

sol = solve(prob; alg_hints = [:stiff], callback = CallbackSet())

@test typeof(sol.alg) <: Rodas4

prob = ODEProblem(f, rand(4, 2) .* ones(4, 2) / 2, (0.0, 1.0))

alg, kwargs = default_algorithm(prob; alg_hints = [:stiff])

@test typeof(alg) <: Rodas4

m = 1.0
ω = 1.0

function mass_system!(du, u, p, t)
    # a(t) = (1/m) w^2 x
    (1 / m) * (ω^2) * u[1]
end

v0 = 0.0
u0 = 1.0
tspan = (0.0, 10.0)

prob = SecondOrderODEProblem(mass_system!, v0, u0, tspan)
sol = solve(prob)

@test typeof(sol.alg) <: Tsit5

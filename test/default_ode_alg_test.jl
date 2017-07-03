using DifferentialEquations, DiffEqProblemLibrary, Base.Test

alg, kwargs = default_algorithm(prob_ode_2Dlinear;dt=1//2^(4))
sol =solve(prob_ode_2Dlinear;dt=1//2^(4))

@test typeof(sol.alg) == typeof(alg)
@test typeof(sol.alg) <: Tsit5

sol =solve(prob_ode_2Dlinear;reltol=1e-1)

@test typeof(sol.alg) <: BS3

sol =solve(prob_ode_2Dlinear;reltol=1e-7)

@test typeof(sol.alg) <: Vern7

sol =solve(prob_ode_2Dlinear;reltol=1e-10)

@test typeof(sol.alg) <: Vern9

sol =solve(prob_ode_2Dlinear;alg_hints=[:stiff])

@test typeof(sol.alg) <: CVODE_BDF

sol =solve(prob_ode_2Dlinear;alg_hints=[:stiff],reltol=1e-1)

@test typeof(sol.alg) <: Rosenbrock23

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
(p::typeof(f))(::Type{Val{:analytic}},t,u0) = u0*exp(linear_bigα*t)
prob_ode_bigfloat2Dlinear = ODEProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,(0.0,1.0))

sol =solve(prob_ode_bigfloat2Dlinear;dt=1//2^(4))
@test typeof(sol.alg) <: Vern8

default_algorithm(prob_ode_bigfloat2Dlinear;alg_hints=[:stiff])

sol =solve(prob_ode_bigfloat2Dlinear;alg_hints=[:stiff])

@test typeof(sol.alg) <: Rosenbrock23

sol =solve(prob_ode_bigfloat2Dlinear,nothing;alg_hints=[:stiff])

@test typeof(sol.alg) <: Rosenbrock23

immutable FooAlg end

@test_throws ErrorException solve(prob_ode_bigfloat2Dlinear,FooAlg();default_set=true)

immutable FooAlg2 <: DEAlgorithm end

@test_throws ErrorException solve(prob_ode_bigfloat2Dlinear,FooAlg2();default_set=true)

prob = ODEProblem(f,rand(4,2).*ones(4,2)/2,(0.0,1.0))

sol =solve(prob;alg_hints=[:stiff])

@test typeof(sol.alg) <: CVODE_BDF

sol =solve(prob;alg_hints=[:stiff],reltol=1e-1)

@test typeof(sol.alg) <: Rosenbrock23

sol =solve(prob;alg_hints=[:stiff],callback=CallbackSet())

@test typeof(sol.alg) <: Rosenbrock23

prob_mm = ODEProblem(f,rand(4,2).*ones(4,2)/2,(0.0,1.0),mass_matrix=nothing)

alg, kwargs = default_algorithm(prob_mm;alg_hints=[:stiff])

@test typeof(alg) <: Rosenbrock23

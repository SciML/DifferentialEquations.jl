using DifferentialEquations, DiffEqProblemLibrary, Base.Test

alg, kwargs = default_algorithm(prob_ode_2Dlinear;dt=1//2^(4))
sol =solve(prob_ode_2Dlinear;dt=1//2^(4))

@test typeof(sol.alg) == typeof(alg)
@test typeof(sol.alg) == Tsit5

sol =solve(prob_ode_2Dlinear;alg_hints=[:stiff])

@test typeof(sol.alg) == CVODE_BDF{:Newton,:Dense}

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
(p::typeof(f))(::Type{Val{:analytic}},t,u0) = u0*exp(linear_bigα*t)
prob_ode_bigfloat2Dlinear = ODEProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,(0.0,1.0))

sol =solve(prob_ode_bigfloat2Dlinear;dt=1//2^(4))
@test typeof(sol.alg) == Vern7

default_algorithm(prob_ode_bigfloat2Dlinear;alg_hints=[:stiff])

sol =solve(prob_ode_bigfloat2Dlinear;alg_hints=[:stiff])

@test typeof(sol.alg) <: Rosenbrock23

immutable FooAlg end

@test_throws ErrorException solve(prob_ode_bigfloat2Dlinear,FooAlg;default_set=true)

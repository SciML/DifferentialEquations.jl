## Breakout these since no other test of their adaptivity

using DifferentialEquations

Δts = 1.//2.^(8:-1:4)
testTol = 0.2
bools = Vector{Bool}(0)
### Rosenbrock23

prob = prob_ode_linear

sim = test_convergence(Δts,prob,alg=:Rosenbrock23)
push!(bools,abs(sim.𝒪est[:final]-2) < testTol)

sol = solve(prob,[0;1],alg=:Rosenbrock23)
push!(bools,length(sol) < 20)

prob = prob_ode_2Dlinear

sim = test_convergence(Δts,prob,alg=:Rosenbrock23)
push!(bools,abs(sim.𝒪est[:final]-2) < testTol)

sol = solve(prob,[0;1],alg=:Rosenbrock23)
push!(bools,length(sol) < 20)

### Rosenbrock32

prob = prob_ode_linear

sim = test_convergence(Δts,prob,alg=:Rosenbrock32)
push!(bools,abs(sim.𝒪est[:final]-3) < testTol)

sol = solve(prob,[0;1],alg=:Rosenbrock32)
push!(bools,length(sol) < 20)

prob = prob_ode_2Dlinear

sim = test_convergence(Δts,prob,alg=:Rosenbrock32)
push!(bools,abs(sim.𝒪est[:final]-3) < testTol)

sol = solve(prob,[0;1],alg=:Rosenbrock32)
push!(bools,length(sol) < 20)

minimum(bools)

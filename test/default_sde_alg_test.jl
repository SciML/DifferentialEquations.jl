using DifferentialEquations, DiffEqProblemLibrary, Base.Test

using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
import DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_additive

srand(100)

prob = prob_sde_additive
sol =solve(prob,dt=1/2^(3))
@test typeof(sol.alg) <: SOSRI

sol =solve(prob,dt=1/2^(3),alg_hints=[:additive])
@test typeof(sol.alg) <: SOSRA

sol =solve(prob,dt=1/2^(3),alg_hints=[:stratonovich])
@test StochasticDiffEq.alg_interpretation(sol.alg) == :stratonovich
@test typeof(sol.alg) <: RKMil

f = (du,u,p,t) -> du.=1.01u
g = function (du,u,p,t)
  du[1,1] = 0.3u[1]
  du[1,2] = 0.6u[1]
  du[1,3] = 0.9u[1]
  du[1,4] = 0.12u[2]
  du[2,1] = 1.2u[1]
  du[2,2] = 0.2u[2]
  du[2,3] = 0.3u[2]
  du[2,4] = 1.8u[2]
end
prob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=zeros(2,4))

sol =solve(prob,dt=1/2^(3))
@test typeof(sol.alg) <: LambaEM

sol =solve(prob,dt=1/2^(3),alg_hints=[:stiff])
@test typeof(sol.alg) <: ISSEM

sol =solve(prob,dt=1/2^(3),alg_hints=[:additive])
@test typeof(sol.alg) <: SOSRA

sol =solve(prob,dt=1/2^(3),alg_hints=[:stratonovich])
@test typeof(sol.alg) <: LambaEulerHeun

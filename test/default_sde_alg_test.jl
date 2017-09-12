using DifferentialEquations, DiffEqProblemLibrary, Base.Test

srand(100)

prob = prob_sde_additive
sol =solve(prob,dt=1/2^(3))
@test typeof(sol.alg) <: SRIW1

sol =solve(prob,dt=1/2^(3),alg_hints=[:additive])
@test typeof(sol.alg) <: SRA1

sol =solve(prob,dt=1/2^(3),alg_hints=[:stratonovich])
@test StochasticDiffEq.alg_interpretation(sol.alg) == :stratonovich
@test typeof(sol.alg) <: RKMil

f = (t,u,du) -> du.=1.01u
g = function (t,u,du)
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

sol =solve(prob,dt=1/2^(3),alg_hints=[:additive])
@test typeof(sol.alg) <: EM

sol =solve(prob,dt=1/2^(3),alg_hints=[:stratonovich])
@test typeof(sol.alg) <: EulerHeun

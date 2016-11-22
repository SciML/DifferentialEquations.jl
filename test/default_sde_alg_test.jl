using DifferentialEquations

srand(100)

prob = prob_sde_additive
sol =solve(prob,dt=1/2^(3))
@test typeof(sol.alg) == SRIW1

sol =solve(prob,dt=1/2^(3),alg_hints=[:additive])
@test typeof(sol.alg) == SRA1

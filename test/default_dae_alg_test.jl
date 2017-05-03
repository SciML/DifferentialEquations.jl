using DifferentialEquations

prob = prob_dae_resrob
sol =solve(prob)

@test typeof(sol.alg) <: IDA

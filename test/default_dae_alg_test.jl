using DifferentialEquations, Test, Sundials

using DiffEqProblemLibrary.DAEProblemLibrary: importdaeproblems; importdaeproblems()
import DiffEqProblemLibrary.DAEProblemLibrary: prob_dae_resrob

prob = prob_dae_resrob
sol = solve(prob)

@test typeof(sol.alg) <: Sundials.IDA

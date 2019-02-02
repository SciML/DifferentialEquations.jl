using DifferentialEquations, Test

using DiffEqProblemLibrary.DAEProblemLibrary: importdaeproblems; importdaeproblems()
import DiffEqProblemLibrary.DAEProblemLibrary: prob_dae_resrob

prob = prob_dae_resrob
sol = solve(prob)

@test typeof(sol.alg) <: IDA

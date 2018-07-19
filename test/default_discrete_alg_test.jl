using DifferentialEquations, Test

prob = DiscreteProblem(zeros(2),(0.0,1.0))
sol = solve(prob)
@test typeof(sol.alg) <: FunctionMap

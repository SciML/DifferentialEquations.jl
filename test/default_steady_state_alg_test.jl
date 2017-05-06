using DifferentialEquations, Base.Test

f = function (t,u,du)
  du[1] = 2 - 2u[1]
  du[2] = u[1] - 4u[2]
end
u0 = zeros(2)
prob = SteadyStateProblem(f,u0)

sol = solve(prob)

@test typeof(sol.alg)<: SSRootfind

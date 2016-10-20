using DifferentialEquations

f = function (t,u,du)
  du[1][1] = 0.85u[1][1]
end

u = [[1]]

prob = ODEProblem(f,u)

sol = solve(prob,Δt=1/2^4)

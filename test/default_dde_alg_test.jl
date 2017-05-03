using DifferentialEquations


lags = [.2]
f = function (t,u,h)
  out = -h(t-.2) + u
end
h = (t) -> 0.0


prob = ConstantLagDDEProblem(f,h,1.0,lags,(0.0,10.0))

sol = solve(prob)

@test typeof(sol.alg) <: Tsit5

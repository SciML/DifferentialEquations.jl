using DifferentialEquations


lags = [.2]
f = function (t,u,h)
  out = -h(t-.2) + u
end
h = (t) -> 0.0


prob = DDEProblem(f,h,1.0,(0.0,10.0),lags)

sol = solve(prob)

@test typeof(sol.alg) <: Tsit5

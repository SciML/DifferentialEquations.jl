using DifferentialEquations, ParameterizedFunctions

f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1

prob = ODEProblem(f,[1.0;1.0])
sol = solve(prob,[0;10];sensitivity_params=[:a,:b])
plot(sol,sensitivity=true)

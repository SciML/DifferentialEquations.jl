using DifferentialEquations, NLsolve, ParameterizedFunctions

f = @ode_def BallBounce begin
  dy =  v
  dv = -g
end g=9.81

function event_f(t,u) # Event when event_f(t,u,k) == 0
  u[1]
end

function apply_event!(u)
  u[2] = -u[2]
end

const Δt_safety = 1
const interp_points = 10

callback = @ode_callback begin
  @ode_event event_f apply_event! interp_points Δt_safety
end

u0 = [50.0,0.0]
prob = ODEProblem(f,u0)
tspan = [0;15]

sol = solve(prob,tspan,callback=callback)
Plots.plotly()
plot(sol,denseplot=true)

sol = solve(prob,tspan,callback=callback,alg=:Vern6)
#plot(sol,denseplot=true)

bounced = ODEProblem(f,sol[8])
sol_bounced = solve(bounced,tspan,callback=callback,alg=:Vern6,Δt=sol.t[9]-sol.t[8])
#plot(sol_bounced,denseplot=true)
sol_bounced(0.04) # Complete density
bool1 = maximum(maximum.(map((i)->sol.k[9][i]-sol_bounced.k[2][i],1:length(sol.k[9])))) == 0


sol2= solve(prob,tspan,callback=callback,alg=:Vern6,adaptive=false,Δt=1/2^4)
#plot(sol2)

sol2= solve(prob,tspan,alg=:Vern6)

sol3= solve(prob,tspan,alg=:Vern6,saveat=[.5])

bool2 = sol2(3) ≈ sol(3)

bool1 && bool2

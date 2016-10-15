using DifferentialEquations, NLsolve

const α = 0.3
f = function (t,u,du)
  for i in 1:length(u)
    du[i] = α*u[i]
  end
end

function event_f(t,u) # Event when event_f(t,u,k) == 0
  1-maximum(u)
end

function apply_event!(u,cache)
  @ode_change_cachesize cache length+1
  maxidx = findmax(u)[2]
  Θ = rand()
  u[maxidx] = Θ
  u[end] = 1-Θ
end

const Δt_safety = 1
const interp_points = 10
callback = @ode_callback begin
  @ode_event event_f apply_event! interp_points Δt_safety
end
u0 = [0.2]
prob = ODEProblem(f,u0)
tspan = [0;10]
sol = solve(prob,tspan,callback=callback)

#=
Plots.plotlyjs()
plot(sol)

plot(sol.t,map((x)->length(x),sol[:]),lw=3,
     ylabel="Number of Cells",xlabel="Time")
ts = linspace(0,10,100)
plot(ts,map((x)->x[1],sol.(ts)),lw=3,
     ylabel="Amount of X in Cell 1",xlabel="Time")
=#

true

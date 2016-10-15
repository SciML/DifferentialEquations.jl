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

@noinline function apply_event!(u,cache)
  #@ode_change_cachesize cache length+1
  maxidx = findmax(u)[2]
  Θ = rand()
  u[maxidx] = Θ
  u[end] = 1-Θ
end

const Δt_safety = 1
const interp_points = 10
const rootfind_event_loc = true
callback = @ode_callback begin
  @ode_event event_f apply_event! rootfind_event_loc interp_points Δt_safety
end
u0 = [0.2]
prob = ODEProblem(f,u0)
tspan = [0;10]
sol = solve(prob,tspan,callback=callback)

for alg in DifferentialEquations.DIFFERENTIALEQUATIONSJL_ALGORITHMS
  if !contains(string(alg),"Vectorized") && !contains(string(alg),"Threaded") && alg ∉ DifferentialEquations.DIFFERENTIALEQUATIONSJL_IMPLICITALGS
    println(alg)
    sol = solve(prob,tspan,callback=callback,alg=alg)
  end
end

callback_no_interp = @ode_callback begin
  @ode_event event_f apply_event! false 0
end

for alg in DifferentialEquations.DIFFERENTIALEQUATIONSJL_ALGORITHMS
  if !contains(string(alg),"Vectorized") && !contains(string(alg),"Threaded") && alg ∉ DifferentialEquations.DIFFERENTIALEQUATIONSJL_IMPLICITALGS
    println(alg)
    sol = solve(prob,tspan,callback=callback_no_interp,alg=alg,dense=false)
  end
end

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

using DifferentialEquations

## Setup Tests

prob = prob_ode_linear

probs = Vector{ODEProblem}(2)
probs[1] = prob
probs[2] = prob_ode_2Dlinear

tspan = [0,1]
tspans = Vector{Vector{Int64}}(2)
tspans[1] = tspan; tspans[2] = tspan
abstols = 1./10.^(3:10)
reltols = 1./10.^(3:10)


setups = [Dict(:alg=>:RK4);Dict(:alg=>:Euler);Dict(:alg=>:BS3);
          Dict(:alg=>:Midpoint);Dict(:alg=>:BS5);Dict(:alg=>:DP5)]


t1 = @elapsed sol = solve(prob::ODEProblem,tspan,alg=:RK4,Δt=1/2^(4))
t2 = @elapsed sol2 = solve(prob::ODEProblem,tspan;Δt=1/2^(4),setups[1]...)

bool = (sol2[end] == sol[end])

## Shootout Tests

println("Shootout Tests")

shoot = ode_shootout(prob,tspan,setups,Δt=1/2^(4))

show(shoot)
println(shoot)
shoot[end]

set = ode_shootoutset(probs,tspans,setups;Δt=1/2^(4))

println(set[1])
println(set[:])
set[end]
set[1][:]

## WorkPrecision Tests
println("WorkPrecision Tests")
wp = ode_workprecision(prob,tspan,abstols,reltols;alg=:DP5,name="Dormand-Prince 4/5")

wp[1]
wp[:]
wp[end]
show(wp)

wp_set = ode_workprecision_set(prob,tspan,abstols,reltols,setups;Δt=1/2^4,numruns=2)

wp_set[1]
wp_set[:]
wp_set[end]
println(wp_set)
show(wp_set)
bool2 = (minimum(diff(wp_set[1].errors).==0)) # The errors for a fixed timestep method should be constant

bool && bool2

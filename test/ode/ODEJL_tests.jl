using DifferentialEquations

prob = linearODEExample()
Δt = 1/2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.

sol =solve(prob::ODEProblem,tspan;Δt=Δt,save_timeseries=true,alg=:ode1)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:ode23)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode45)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode78)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode23s)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode2_midpoint)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode2_heun)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode4)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode45_fe)

prob = twoDimlinearODEExample()

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode1)
TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,tspan,Δt=Δt,save_timeseries=true,alg=:ode23)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode45)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode78)

#sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode23s) #ODE.jl issues
#TEST_PLOT && plot(sol,plot_analytic=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode2_midpoint)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode2_heun)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode4)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,save_timeseries=true,alg=:ode45_fe)

true

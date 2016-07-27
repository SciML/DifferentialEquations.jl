using DifferentialEquations

prob = linearODEExample()
Δt = 1/2^(4) #The initial timestepping size. It will automatically assigned if not given.
tspan = [0,1] # The timespan. This is the default if not given.

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode1)
plot(sol,plottrue=true)
Plots.gui()

sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:ode23)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode45)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode78)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode23s)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode2_midpoint)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode2_heun)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode4)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode45_fe)
plot(sol,plottrue=true)

prob = twoDimlinearODEExample()

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode1)
plot(sol,plottrue=true)
Plots.gui()

sol =solve(prob::ODEProblem,tspan,Δt=Δt,fullSave=true,alg=:ode23)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode45)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode78)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode23s) #Issues
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode2_midpoint)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode2_heun)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode4)
plot(sol,plottrue=true)

sol =solve(prob::ODEProblem,[0,1];Δt=Δt,fullSave=true,alg=:ode45_fe)
plot(sol,plottrue=true)

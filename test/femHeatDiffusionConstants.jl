######
##FEM Heat Diffusion Constants Tests
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = .5
Δx = 1//2^(4)
Δt = 1//2^(4)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = heatProblemExample_diffusionconstants()

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,save_timeseries=true,timeseries_steps=1)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler,save_timeseries=true,timeseries_steps=1)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()

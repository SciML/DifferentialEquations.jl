######
##FEM Heat Diffusion Constants Tests
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = 100
Δx = 1//2^(3)
Δt = 1//2^(9)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = heatProblemExample_diffusionconstants(max=100)

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler,save_timeseries=true,timeseries_steps=100)

plot(sol,plottrue=false,zlim=(0,3),cbar=false)
gui()

@linux? animate(sol,zlim=(0,3),cbar=false) : println("Animation only works with ImageMagick installation, disabled on osx for testing")

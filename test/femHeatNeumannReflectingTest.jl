######
##FEM Heat Diffusion Constants Tests
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = 1000
Δx = 1//2^(4)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_diffusionconstants(max=100)

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=1)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()

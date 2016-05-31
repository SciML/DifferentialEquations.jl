######
##FEM Heat Diffusion Constants Tests
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = .5
Δx = 1//2^(4)
Δt = 1//2^(4)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_diffusionconstants()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="ImplicitEuler",fullSave=true,saveSteps=1)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=1)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()

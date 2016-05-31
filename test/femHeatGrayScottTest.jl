######
##FEM Heat Nonlinear Test
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = 100
Δx = 1//2^(4)
Δt = 1//2^(14)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_grayscott(D = [1 .0001])

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,saveSteps=1)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()
CFLμ(Δt,Δx)

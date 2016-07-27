######
##FEM Heat Nonlinear Test
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = 10
Δx = 1//2^(4)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_grayscott(;ρ=.3,k=.62,D=[1e-2 .5e-2])

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg=:Euler,fullSave=true,saveSteps=1000)

plot(sol,plottrue=false,zlim=(0,1),cbar=false)
gui()
CFLμ(Δt,Δx)
#animate(sol::FEMSolution;zlim=(0,3),cbar=false)

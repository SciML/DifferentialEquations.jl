######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 5
Δx = 1//2^(1)
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_birthdeathinteractingsystem()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="SemiImplicitEuler")

plot(sol,plottrue=false,zlim=(0,2),cbar=false)

prob = heatProblemExample_birthdeathsystem()
sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="ImplicitEuler",iterations=1000,autodiff=true)

plot(sol,plottrue=false,zlim=(0,2),cbar=false)

true

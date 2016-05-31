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

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler")

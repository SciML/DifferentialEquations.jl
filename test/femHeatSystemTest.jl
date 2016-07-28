######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 5
Δx = 1//2^(1)
Δt = 1//2^(7)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = heatProblemExample_birthdeathinteractingsystem()

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitEuler)

plot(sol,plottrue=false,zlim=(0,2),cbar=false)

prob = heatProblemExample_birthdeathsystem()
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,iterations=1000,autodiff=true)

plot(sol,plottrue=false,zlim=(0,2),cbar=false)

true

######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 5
Δx = 1/2^(1)
Δt = 1/2^(7)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = prob_femheat_birthdeathinteractingsystem
#@code_warntype solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler)
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler)

TEST_PLOT && plot(sol,plot_analytic=false,zlim=(0,2),cbar=false)

prob = prob_femheat_birthdeathsystem
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,iterations=1000,autodiff=true)

TEST_PLOT && plot(sol,plot_analytic=false,zlim=(0,2),cbar=false)

true

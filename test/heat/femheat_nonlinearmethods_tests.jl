######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = prob_femheat_birthdeath


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler)

println("Semi-implicit Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitEuler)

println("Semi-implicit Crank Nicholson")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitCrankNicholson)

Δx = 1//2^(2)
Δt = 1//2^(4)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
println("Implicit Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,autodiff=true)
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,autodiff=false)
TEST_PLOT && plot(sol)

#Returns true if nonlinear solver is correct
bool1 = maximum(abs.(sol.u - .777))<.01

### Stochastic Tests

#Define a parabolic problem
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = prob_femheat_stochasticbirthdeath


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler)

println("Semi-implicit Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitEuler)

#=
# CG and GMRES require size 1 vector, breaks with numvars change
println("Semi-implicit Crank Nicholson")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitCrankNicholson,solver=:CG)

println("Semi-implicit Crank Nicholson GMRES")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitCrankNicholson,solver=:GMRES)
=#

#Define a quicker problem
Δx = 1//2^(1)
Δt = 1//2^(1)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
println("Implicit Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,autodiff=true)
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:ImplicitEuler,autodiff=false)
TEST_PLOT && plot(sol)

bool1

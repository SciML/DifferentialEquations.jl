#######
##FEM Stochastic Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
using DifferentialEquations, Plots
T = 5
Δx = 1//2^(3)
Δt = 1//2^(11)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_stochasticbirthdeath()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,solver="LU")

println("Generating Animation")
animate(sol::FEMSolution;zlims=(0,3),cbar=false)

# Returns true if nothing error'd
true

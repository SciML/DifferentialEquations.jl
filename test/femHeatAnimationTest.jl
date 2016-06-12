#######
##FEM Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
using DifferentialEquations
T = 2
Δx = 1//2^(3)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
prob = heatProblemExample_moving()

sim = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true)

println("Generating Animation")
animate(sim::FEMSolution;zlims=(0,.1),cbar=false)

## Should have moved off the frame.
maximum(sim.u) .< 1e-6

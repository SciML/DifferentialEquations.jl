#######
##FEM Stochastic Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
using DifferentialEquations, Plots#, ImageMagick
T = 5
Δx = 1//2^(3)
Δt = 1//2^(5)
fem_mesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,:neumann)
prob = heatProblemExample_stochasticbirthdeath()

sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:SemiImplicitCrankNicholson,save_timeseries=true,solver=:LU)

println("Generating Animation")
TEST_PLOT && animate(sol::FEMSolution;zlims=(0,3),cbar=false)

# Returns true if nothing error'd
true

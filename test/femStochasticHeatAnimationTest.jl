#######
##FEM Stochastic Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
#Will work on Windows, but will give blurry output
using DifferentialEquations
T = 5
Δx = 1//2^(3)
Δt = 1//2^(11)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_stochasticbirthdeath()

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg="Euler",fullSave=true,solver="LU")

println("Generating Animation")
@linux? animate(sol::FEMSolution;zlim=(0,3),cbar=false) : println("Animation only works with ImageMagick installation, disabled on osx for testing")

# Returns true if nothing error'd
true

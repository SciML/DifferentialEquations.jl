#######
##FEM Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
#Will work on Windows, but will give blurry output
using DifferentialEquations
T = 2
Δx = 1//2^(3)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_moving()

res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler",fullSave=true)

println("Generating Animation")
@linux? solplot_animation(res::FEMSolution;zlim=(0,.1),cbar=false) : println("Animation only works with ImageMagick installation, disabled on osx for testing")

## Should have moved off the frame.
maximum(res.u) .< 1e-6

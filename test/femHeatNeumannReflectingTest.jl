######
##FEM Heat Diffusion Constants Tests
######
using DifferentialEquations, Plots

#Define a parabolic problem
T = 100
Δx = 1//2^(3)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
prob = heatProblemExample_diffusionconstants(max=100)

sol = solve(femMesh::FEMmesh,prob::HeatProblem,alg=:Euler,fullSave=true,saveSteps=100)

plot(sol,plottrue=false,zlim=(0,3),cbar=false)
gui()

@linux? animate(sol,zlim=(0,3),cbar=false) : println("Animation only works with ImageMagick installation, disabled on osx for testing")

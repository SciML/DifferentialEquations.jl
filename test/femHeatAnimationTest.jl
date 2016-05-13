#######
##FEM Heat Animation Test
#######

#Generates an animation for a solution of the heat equation
#Uses Plots.jl, requires matplotlib >=1.5
#Will work on Windows, but will give blurry output
using DifferentialEquations
T = 2
Δx = 1//2^(5)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_moving()

res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler",fullSave=true)
if !isdefined(:testState) #Don't plot during test
  println("in here!")
  solplot_animation(res::FEMSolution;zlim=(0,.1),vmax=.1,cbar=false) #Make animation, excluded from test
end

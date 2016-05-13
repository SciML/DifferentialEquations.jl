######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
pdeProb = heatProblemExample_birthdeath()


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
if !isdefined(:testState) #Don't plot during test
  Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
  Plots.savefig("plot.pdf")
  Plots.savefig("plot.png")
end


println("Semi-implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitEuler")
if !isdefined(:testState) #Don't plot during test
  Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
end

println("Semi-implicit Crank Nicholson")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitCrankNicholson")
if !isdefined(:testState) #Don't plot during test
  Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
end

Δx = 1//2^(2)
Δt = 1//2^(4)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
println("Implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",autodiff=true)
if !isdefined(:testState) #Don't plot during test
  Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
end

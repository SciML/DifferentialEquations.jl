######
##FEM Heat Nonlinear Test
######
using DifferentialEquations

#Define a parabolic problem
T = 1
Δx = 1//2^(4)
Δt = 1//2^(8)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
pdeProb = heatProblemExample_birthdeath()


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)
Plots.savefig("plot.pdf")
Plots.savefig("plot.png")

println("Semi-implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitEuler")
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)

println("Semi-implicit Crank Nicholson")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitCrankNicholson")
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)

println("Implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",autodiff=true)
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)

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

println("Semi-implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitEuler")

println("Semi-implicit Crank Nicholson")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitCrankNicholson")

Δx = 1//2^(2)
Δt = 1//2^(4)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
println("Implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",autodiff=true)
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)

#Returns true if nonlinear solver is correct
bool1 = maximum(abs(res.u - .777))<.01

### Stochastic Tests

#Define a parabolic problem
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
pdeProb = heatProblemExample_stochasticbirthdeath()


#Solve it with a bunch of different algorithms, plot solution
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")

println("Semi-implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitEuler")

println("Semi-implicit Crank Nicholson")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitCrankNicholson",solver="CG")

println("Semi-implicit Crank Nicholson GMRES")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="SemiImplicitCrankNicholson",solver="GMRES")

#Define a quicker problem
Δx = 1//2^(1)
Δt = 1//2^(1)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Neumann")
println("Implicit Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",autodiff=true)
Plots.surface(femMesh.node[:,1],femMesh.node[:,2],res.u,zlim=(0,2),cbar=false)

bool1

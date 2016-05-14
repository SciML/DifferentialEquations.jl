######
##FEM Heat Method Tests
######
using DifferentialEquations

#Define a parabolic problem
T = 1
Δx = 1//2^(3)
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_moving() #also try heatProblemExample_pure() or heatProblemExample_diffuse()

#Solve it with a bunch of different algorithms, plot solution
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")

Δt = 1//2^(4) #Make faster for tests
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
println("Direct")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="Direct")

println("LU")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="LU")

println("CG")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler")

println("Direct")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson",solver="Direct")

println("Cholesky")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson",solver="Cholesky")

println("CG")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson")

#Define a different parabolic problem
T = 1//2^(5)
Δx = 1//2^(3)
Δt = 1//2^(9)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_pure()

#Solve with Euler Method
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")

#Choose a finer mesh, solve with Euler, and add this result to the previous as
#an approximately true solution.
T = 1//2^(5)
Δt = 1//2^(11)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
res2 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
appxTrue!(res,res2)
solplot_appxvstrue(res)

res.nodeErr2<.005 #Returns true if res solution is near the apprxTrue res2 solution

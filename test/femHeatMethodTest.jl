######
##FEM Heat Method Tests
######
using DifferentialEquations

#Define a parabolic problem
T = 1
Δx = 1//2^(5)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_moving() #also try heatProblemExample_pure() or heatProblemExample_diffuse()

#Solve it with a bunch of different algorithms, plot solution
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
solplot_appxvstrue(res,savefile="plot.svg")

println("Direct")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="Direct")
solplot_appxvstrue(res,savefile="plot.svg")

println("LU")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="LU")
solplot_appxvstrue(res,savefile="plot.svg")

println("CG")
res2 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler")
solplot_appxvstrue(res2,savefile="plot.svg")

println("Direct")
res3 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson",solver="Direct")
solplot_appxvstrue(res3,savefile="plot.svg")

println("Cholesky")
res3 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson",solver="Cholesky")
solplot_appxvstrue(res3,savefile="plot.svg")

println("CG")
res4 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson")
solplot_appxvstrue(res4,savefile="plot.svg")

#Define a different parabolic problem
T = 1//8
Δx = 1//2^(5)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_pure()

#Solve with Euler Method
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
solplot(res,savefile="plot.svg")

#Choose a finer mesh, solve with Euler, and add this result to the previous as
#an approximately true solution.
T = 1//8
Δt = 1//2^(16)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
res2 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
appxTrue!(res,res2)

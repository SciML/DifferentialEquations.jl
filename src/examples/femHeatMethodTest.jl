######
##FEM Heat Method Tests
######
using DiffEq

T = 1
Δx = 1//2^(5)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_moving() #also try heatProblemExample_pure() or heatProblemExample_diffuse()

println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
solplot_appxvstrue(res,savefile="plot.svg")

println("Direct")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="Direct")
solplot_appxvstrue(res,savefile="plot.svg")

println("Cholesky")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="Cholesky")
solplot_appxv strue(res,savefile="plot.svg")

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

T = 1//8
Δx = 1//2^(5)
Δt = 1//2^(12)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_pure()

println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
solplot(res,savefile="plot.svg")

T = 1//8
Δt = 1//2^(16)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
res2 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
appxTrue!(res,res2)

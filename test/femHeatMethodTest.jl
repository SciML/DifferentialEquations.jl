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
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

println("Direct")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="Direct")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

println("LU")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler",solver="LU")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

println("CG")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="ImplicitEuler")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

println("Direct")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson",solver="Direct")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

println("Cholesky")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson",solver="Cholesky")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

println("CG")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="CrankNicholson")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

#Define a different parabolic problem
T = 1//8
Δx = 1//2^(3)
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
pdeProb = heatProblemExample_pure()

#Solve with Euler Method
println("Euler")
res = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
if !isdefined(:testState) #Don't plot during test
  solplot_appxvstrue(res,savefile="plot.svg")
end

#Choose a finer mesh, solve with Euler, and add this result to the previous as
#an approximately true solution.
T = 1//8
Δt = 1//2^(7)
femMesh = parabolic_squaremesh([0 1 0 1],Δx,Δt,T,"Dirichlet")
res2 = fem_solveheat(femMesh::FEMmesh,pdeProb::HeatProblem,alg="Euler")
appxTrue!(res,res2)
